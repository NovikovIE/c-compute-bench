#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

//
// Mersenne Twister (MT19937) a 32-bit pseudorandom number generator.
//
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

//
// Benchmark-specific data and functions
//

typedef struct {
    int num_dimensions;
    int num_iterations;
    double* x;         // Current solution vector
    double* g;         // Gradient vector
    double* p;         // Trial step vector
    double* x_new;     // Candidate new solution
    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

// Calculates objective function (f) and gradient (g) for the Rosenbrock function.
// g can be NULL if only the function value is needed.
double rosenbrock_fg(const double* x, double* g, int n) {
    double f = 0.0;

    if (g) {
        for (int i = 0; i < n; i++) {
            g[i] = 0.0;
        }
    }

    for (int i = 0; i < n - 1; i++) {
        double term1 = x[i+1] - x[i] * x[i];
        double term2 = 1.0 - x[i];
        f += 100.0 * term1 * term1 + term2 * term2;

        if (g) {
            g[i] += -400.0 * x[i] * term1 - 2.0 * term2;
            g[i+1] += 200.0 * term1;
        }
    }
    return f;
}

// Calculates the L2 norm of a vector.
double vec_norm(const double* v, int n) {
    double sum_sq = 0.0;
    for (int i = 0; i < n; i++) {
        sum_sq += v[i] * v[i];
    }
    return sqrt(sum_sq);
}

// Parses command line arguments, allocates memory, and initializes data.
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_dimensions> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_dimensions = atoi(argv[1]);
    g_data.num_iterations = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    g_data.x = (double*)malloc(g_data.num_dimensions * sizeof(double));
    g_data.g = (double*)malloc(g_data.num_dimensions * sizeof(double));
    g_data.p = (double*)malloc(g_data.num_dimensions * sizeof(double));
    g_data.x_new = (double*)malloc(g_data.num_dimensions * sizeof(double));

    if (!g_data.x || !g_data.g || !g_data.p || !g_data.x_new) {
        fprintf(stderr, "Failed to allocate memory.\n");
        exit(1);
    }

    // Initialize the starting point x with random values in [-2, 2]
    for (int i = 0; i < g_data.num_dimensions; i++) {
        g_data.x[i] = ((double)mt_rand() / (double)UINT32_MAX) * 4.0 - 2.0;
    }
    g_data.final_result = 0.0;
}

// Executes the core trust-region optimization algorithm.
void run_computation() {
    int n = g_data.num_dimensions;
    
    // Trust-region parameters
    double delta = 1.0;          // Trust region radius
    const double max_delta = 10.0;    // Maximum radius
    const double eta = 0.1;           // Step acceptance threshold
    
    for (int i = 0; i < g_data.num_iterations; i++) {
        double f_curr = rosenbrock_fg(g_data.x, g_data.g, n);
        double g_norm = vec_norm(g_data.g, n);

        if (g_norm < 1e-9) {
            break; // Converged
        }

        // Solve the subproblem using the Cauchy point (steepest descent step)
        // constrained by the trust region radius.
        double step_len = delta / g_norm;
        for (int j = 0; j < n; j++) {
            g_data.p[j] = -step_len * g_data.g[j];
            g_data.x_new[j] = g_data.x[j] + g_data.p[j];
        }

        // Evaluate a
        double f_new = rosenbrock_fg(g_data.x_new, NULL, n);

        // Calculate actual vs. predicted reduction
        double actual_reduction = f_curr - f_new;
        // Predicted reduction: -g'p = -g'(-step_len*g) = step_len * ||g||^2 = delta * ||g||
        double pred_reduction = delta * g_norm;

        double rho = (pred_reduction > 1e-9) ? (actual_reduction / pred_reduction) : -1.0;

        // Update trust region radius
        if (rho < 0.25) {
            delta *= 0.25;
        } else if (rho > 0.75) {
            delta = fmin(2.0 * delta, max_delta);
        }

        // Accept or reject the step
        if (rho > eta) {
            memcpy(g_data.x, g_data.x_new, n * sizeof(double));
        }
    }

    // Store the final objective function value as the result.
    g_data.final_result = rosenbrock_fg(g_data.x, NULL, n);
}

// Frees all memory allocated in setup_benchmark.
void cleanup() {
    free(g_data.x);
    free(g_data.g);
    free(g_data.p);
    free(g_data.x_new);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print the final result to stdout
    printf("%f\n", g_data.final_result);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
