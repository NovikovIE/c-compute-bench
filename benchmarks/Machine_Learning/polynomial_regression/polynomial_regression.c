#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- START MERSENNE TWISTER --- 
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
             fprintf(stderr, "FATAL: Mersenne Twister not seeded.\n");
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND PARAMETERS ---
int num_samples;
int num_features;
int degree;
int poly_features_count;

double *X;      // Input features: num_samples x num_features
double *y;      // Target values: num_samples
double *X_poly; // Polynomial features: num_samples x poly_features_count
double *XTX;    // X_poly^T * X_poly: poly_features_count x poly_features_count
double *XTy;    // X_poly^T * y: poly_features_count
double *theta;  // Coefficients to be solved: poly_features_count
double final_result = 0.0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_samples num_features degree seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    degree = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    poly_features_count = 1 + num_features * degree;

    // Allocate memory for matrices and vectors
    X = (double*)malloc(num_samples * num_features * sizeof(double));
    y = (double*)malloc(num_samples * sizeof(double));
    X_poly = (double*)malloc(num_samples * poly_features_count * sizeof(double));
    XTX = (double*)malloc(poly_features_count * poly_features_count * sizeof(double));
    XTy = (double*)malloc(poly_features_count * sizeof(double));
    theta = (double*)malloc(poly_features_count * sizeof(double));

    if (!X || !y || !X_poly || !XTX || !XTy || !theta) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Initialize input data with random values
    for (int i = 0; i < num_samples * num_features; ++i) {
        X[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
    for (int i = 0; i < num_samples; ++i) {
        y[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

void run_computation() {
    // Step 1: Construct the polynomial features matrix X_poly
    for (int i = 0; i < num_samples; ++i) {
        X_poly[i * poly_features_count] = 1.0; // Bias term
        for (int j = 0; j < num_features; ++j) {
            double current_power = 1.0;
            double base_val = X[i * num_features + j];
            for (int d = 1; d <= degree; ++d) {
                current_power *= base_val;
                int col_idx = 1 + j * degree + (d - 1);
                X_poly[i * poly_features_count + col_idx] = current_power;
            }
        }
    }

    // Step 2: Compute XTX = X_poly^T * X_poly
    for (int i = 0; i < poly_features_count; ++i) {
        for (int j = i; j < poly_features_count; ++j) { // Exploit symmetry
            double sum = 0.0;
            for (int k = 0; k < num_samples; ++k) {
                sum += X_poly[k * poly_features_count + i] * X_poly[k * poly_features_count + j];
            }
            XTX[i * poly_features_count + j] = sum;
            XTX[j * poly_features_count + i] = sum; // Symmetric matrix
        }
    }

    // Step 3: Compute XTy = X_poly^T * y
    for (int i = 0; i < poly_features_count; ++i) {
        double sum = 0.0;
        for (int k = 0; k < num_samples; ++k) {
            sum += X_poly[k * poly_features_count + i] * y[k];
        }
        XTy[i] = sum;
    }

    // Step 4: Solve XTX * theta = XTy using Gaussian elimination with partial pivoting
    int n = poly_features_count;
    for (int i = 0; i < n; i++) {
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(XTX[k * n + i]) > fabs(XTX[max_row * n + i])) {
                max_row = k;
            }
        }

        for (int k = i; k < n; k++) {
            double temp = XTX[i * n + k];
            XTX[i * n + k] = XTX[max_row * n + k];
            XTX[max_row * n + k] = temp;
        }
        double temp = XTy[i];
        XTy[i] = XTy[max_row];
        XTy[max_row] = temp;

        for (int k = i + 1; k < n; k++) {
            if (XTX[i * n + i] == 0.0) continue; // Matrix is singular
            double factor = XTX[k * n + i] / XTX[i * n + i];
            XTy[k] -= factor * XTy[i];
            for (int j = i; j < n; j++) {
                XTX[k * n + j] -= factor * XTX[i * n + j];
            }
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        if (XTX[i * n + i] == 0.0) { // Singular matrix
            theta[i] = 0.0;
            continue;
        }
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += XTX[i * n + j] * theta[j];
        }
        theta[i] = (XTy[i] - sum) / XTX[i * n + i];
    }

    // Step 5: Accumulate result to prevent dead code elimination
    for (int i = 0; i < poly_features_count; ++i) {
        final_result += theta[i];
    }
}

void cleanup() {
    free(X);
    free(y);
    free(X_poly);
    free(XTX);
    free(XTy);
    free(theta);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
