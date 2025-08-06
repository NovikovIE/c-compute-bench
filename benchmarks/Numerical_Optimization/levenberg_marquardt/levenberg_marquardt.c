/*
 * levenberg_marquardt: A benchmark for non-linear least squares optimization.
 *
 * This program implements the Levenberg-Marquardt algorithm (LMA) to solve a
 * non-linear curve-fitting problem. LMA is an iterative method that finds the
 * parameters of a model function that best fit a set of data points.
 * It dynamically interpolates between the Gauss-Newton method and gradient descent,
 * making it robust and efficient for a wide range of problems.
 *
 * The benchmark is configured with:
 * - num_parameters: The number of parameters in the model function to be optimized.
 * - num_iterations: The number of optimization steps to perform.
 * - num_data_points: The number of (x, y) data points to fit.
 *
 * The core computation involves repeatedly:
 * 1. Calculating the Jacobian matrix of the model function (approximated via finite differences).
 * 2. Forming and solving the LMA augmented normal equations: (J^T*J + lambda*I) * delta = J^T*r
 *    This step is computationally intensive, involving matrix multiplication and solving a linear system.
 * 3. Updating the model parameters.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define MT_N 624
#define MT_M 397
#define MT_MATRIX_A 0x9908b0dfUL
#define MT_UPPER_MASK 0x80000000UL
#define MT_LOWER_MASK 0x7fffffffUL

static uint32_t mt[MT_N];
static int mt_index = MT_N + 1;

void mt_seed(uint32_t seed) {
    mt[0] = seed;
    for (mt_index = 1; mt_index < MT_N; mt_index++) {
        mt[mt_index] = (1812433253UL * (mt[mt_index - 1] ^ (mt[mt_index - 1] >> 30)) + mt_index);
    }
}

uint32_t mt_rand(void) {
    uint32_t y;
    static const uint32_t mag01[2] = {0x0UL, MT_MATRIX_A};
    if (mt_index >= MT_N) {
        if (mt_index > MT_N) {
             fprintf(stderr, "FATAL: Mersenne Twister not seeded.");
             exit(1);
        }
        for (int i = 0; i < MT_N - MT_M; i++) {
            y = (mt[i] & MT_UPPER_MASK) | (mt[i + 1] & MT_LOWER_MASK);
            mt[i] = mt[i + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (int i = MT_N - MT_M; i < MT_N - 1; i++) {
            y = (mt[i] & MT_UPPER_MASK) | (mt[i + 1] & MT_LOWER_MASK);
            mt[i] = mt[i + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[MT_N - 1] & MT_UPPER_MASK) | (mt[0] & MT_LOWER_MASK);
        mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
        mt_index = 0;
    }
    y = mt[mt_index++];
    y ^= (y >> 11); y ^= (y << 7) & 0x9d2c5680UL; y ^= (y << 15) & 0xefc60000UL; y ^= (y >> 18);
    return y;
}

// Generate a random double between 0 and 1
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- BENCHMARK DATA AND PARAMETERS ---
typedef struct {
    int num_parameters;
    int num_iterations;
    int num_data_points;

    double *params;         // Current parameters being optimized
    double *data_x;
    double *data_y;

    // Workspace arrays for computation
    double *residuals;
    double *jacobian;
    double *jtj;            // J^T * J
    double *jtr;            // J^T * r
    double *jtj_aug;        // Augmented J^T * J for solving
    double *jtj_inv;        // Inverse of augmented J^T * J
    double *delta;          // Parameter update step
    double *next_params;

    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

// --- MODEL FUNCTION ---
// A sum of Gaussian peaks. The number of parameters must be a multiple of 3.
// p[3*i+0] = amplitude, p[3*i+1] = center, p[3*i+2] = width
double model_function(double x, const double *p, int n_params) {
    double y = 0.0;
    for (int i = 0; i < n_params / 3; ++i) {
        double amp = p[3 * i];
        double cen = p[3 * i + 1];
        double wid = p[3 * i + 2];
        if (wid == 0.0) continue; // Avoid division by zero
        double exp_arg = (x - cen) / wid;
        y += amp * exp(-exp_arg * exp_arg);
    }
    return y;
}

// --- LINEAR ALGEBRA HELPERS ---
// Invert a square matrix using Gauss-Jordan elimination
int invert_matrix(const double* input, double* output, int n) {
    double *m = (double*)malloc(n * n * sizeof(double));
    if (!m) return 0;
    memcpy(m, input, n * n * sizeof(double));

    // Initialize output to identity matrix
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            output[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < n; ++i) {
        // Find pivot
        int pivot = i;
        for (int j = i + 1; j < n; ++j)
            if (fabs(m[j * n + i]) > fabs(m[pivot * n + i]))
                pivot = j;

        // Swap rows
        if (pivot != i) {
            for (int k = 0; k < n; ++k) {
                double temp = m[i * n + k];
                m[i * n + k] = m[pivot * n + k];
                m[pivot * n + k] = temp;
                temp = output[i * n + k];
                output[i * n + k] = output[pivot * n + k];
                output[pivot * n + k] = temp;
            }
        }

        double div = m[i * n + i];
        if (fabs(div) < 1e-12) { free(m); return 0; } // Singular matrix

        // Normalize pivot row
        for (int j = 0; j < n; ++j) {
            m[i * n + j] /= div;
            output[i * n + j] /= div;
        }

        // Eliminate other rows
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                double mult = m[j * n + i];
                for (int k = 0; k < n; ++k) {
                    m[j * n + k] -= mult * m[i * n + k];
                    output[j * n + k] -= mult * output[i * n + k];
                }
            }
        }
    }

    free(m);
    return 1;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_parameters num_iterations num_data_points seed\n", argv[0]);
        exit(1);
    }

    g_data.num_parameters = atoi(argv[1]);
    g_data.num_iterations = atoi(argv[2]);
    g_data.num_data_points = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (g_data.num_parameters % 3 != 0 || g_data.num_parameters <= 0) {
        fprintf(stderr, "Error: num_parameters must be a positive multiple of 3.\n");
        exit(1);
    }

    // Allocate all memory
    const int p = g_data.num_parameters;
    const int d = g_data.num_data_points;
    g_data.params = (double*)malloc(p * sizeof(double));
    g_data.data_x = (double*)malloc(d * sizeof(double));
    g_data.data_y = (double*)malloc(d * sizeof(double));
    g_data.residuals = (double*)malloc(d * sizeof(double));
    g_data.jacobian = (double*)malloc(d * p * sizeof(double));
    g_data.jtj = (double*)malloc(p * p * sizeof(double));
    g_data.jtr = (double*)malloc(p * sizeof(double));
    g_data.jtj_aug = (double*)malloc(p * p * sizeof(double));
    g_data.jtj_inv = (double*)malloc(p * p * sizeof(double));
    g_data.delta = (double*)malloc(p * sizeof(double));
    g_data.next_params = (double*)malloc(p * sizeof(double));
    g_data.final_result = 0.0;

    // Generate true parameters for the model
    double *true_params = (double*)malloc(p * sizeof(double));
    for (int i = 0; i < p / 3; ++i) {
        true_params[3 * i] = rand_double() * 5.0 + 1.0;   // Amplitude
        true_params[3 * i + 1] = rand_double() * 10.0; // Center
        true_params[3 * i + 2] = rand_double() * 2.0 + 0.5; // Width
    }

    // Generate data points based on the true model with some noise
    for (int i = 0; i < d; ++i) {
        g_data.data_x[i] = rand_double() * 10.0;
        g_data.data_y[i] = model_function(g_data.data_x[i], true_params, p) + (rand_double() - 0.5) * 0.1;
    }

    // Set initial guess for parameters
    for (int i = 0; i < p; ++i) {
        g_data.params[i] = true_params[i] + (rand_double() - 0.5) * 0.5;
    }

    free(true_params);
}

void run_computation() {
    const int p = g_data.num_parameters;
    const int d = g_data.num_data_points;
    const double h = 1e-6; // Step for finite difference Jacobian

    double lambda = 0.001;
    double current_error;

    for (int iter = 0; iter < g_data.num_iterations; ++iter) {
        // 1. Calculate residuals and current error (sum of squares)
        current_error = 0.0;
        for (int i = 0; i < d; ++i) {
            double model_y = model_function(g_data.data_x[i], g_data.params, p);
            g_data.residuals[i] = g_data.data_y[i] - model_y;
            current_error += g_data.residuals[i] * g_data.residuals[i];
        }

        // 2. Calculate Jacobian using finite differences
        for (int j = 0; j < p; ++j) {
            g_data.params[j] += h;
            for (int i = 0; i < d; ++i) {
                double perturbed_y = model_function(g_data.data_x[i], g_data.params, p);
                g_data.jacobian[i * p + j] = (g_data.residuals[i] - (g_data.data_y[i] - perturbed_y)) / h;
            }
            g_data.params[j] -= h; // Restore original parameter
        }

        // 3. Calculate J^T * J and J^T * r
        for (int i = 0; i < p; ++i) {
            g_data.jtr[i] = 0.0;
            for (int k = 0; k < d; ++k) {
                 g_data.jtr[i] += g_data.jacobian[k * p + i] * g_data.residuals[k];
            }
            for (int j = 0; j < p; ++j) {
                g_data.jtj[i * p + j] = 0.0;
                for (int k = 0; k < d; ++k) {
                    g_data.jtj[i * p + j] += g_data.jacobian[k * p + i] * g_data.jacobian[k * p + j];
                }
            }
        }

        // 4. Solve the augmented system (J^T*J + lambda*I) * delta = J^T*r
        memcpy(g_data.jtj_aug, g_data.jtj, p * p * sizeof(double));
        for (int i = 0; i < p; ++i) {
            g_data.jtj_aug[i * p + i] += lambda;
        }

        if (!invert_matrix(g_data.jtj_aug, g_data.jtj_inv, p)) {
            // Matrix is singular, increase damping and try again
            lambda *= 10;
            continue;
        }

        // Calculate delta = (J^T*J + lambda*I)^-1 * J^T*r
        for (int i = 0; i < p; ++i) {
            g_data.delta[i] = 0.0;
            for (int j = 0; j < p; ++j) {
                g_data.delta[i] += g_data.jtj_inv[i * p + j] * g_data.jtr[j];
            }
        }

        // 5. Evaluate potential new parameters and error
        for (int i = 0; i < p; ++i) {
            g_data.next_params[i] = g_data.params[i] + g_data.delta[i];
        }

        double next_error = 0.0;
        for (int i = 0; i < d; ++i) {
            double model_y = model_function(g_data.data_x[i], g_data.next_params, p);
            double resid = g_data.data_y[i] - model_y;
            next_error += resid * resid;
        }

        // 6. Update parameters and lambda based on success
        if (next_error < current_error) {
            lambda /= 10.0;
            memcpy(g_data.params, g_data.next_params, p * sizeof(double));
        } else {
            lambda *= 10.0;
        }
    }
    
    // Accumulate a final result to prevent dead code elimination.
    for(int i = 0; i < p; ++i) {
        g_data.final_result += g_data.params[i];
    }
}

void cleanup() {
    free(g_data.params);
    free(g_data.data_x);
    free(g_data.data_y);
    free(g_data.residuals);
    free(g_data.jacobian);
    free(g_data.jtj);
    free(g_data.jtr);
    free(g_data.jtj_aug);
    free(g_data.jtj_inv);
    free(g_data.delta);
    free(g_data.next_params);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
