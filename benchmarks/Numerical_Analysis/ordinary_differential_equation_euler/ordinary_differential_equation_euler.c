#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <math.h>

// --- Mersenne Twister (verbatim) ---
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
// --- End Mersenne Twister ---

// Global structure to hold benchmark data
typedef struct {
    char function_str[32];
    int num_params;
    double initial_y;
    double start_x;
    double end_x;
    double step_size;
    uint32_t seed;
    long long num_steps;

    double* params;      // Parameters for the ODE function
    double result_accumulator; // Accumulator to prevent dead code elimination
} BenchmarkData;

static BenchmarkData g_data;

// The ODE function dy/dx = f(x, y)
// This benchmark implements a simple linear function: f(x, y) = p0 + p1*x + p2*y.
// This is the only function type supported, as specified by function_str='linear' and num_params=3.
static inline double ode_function(double x, double y) {
    return g_data.params[0] + g_data.params[1] * x + g_data.params[2] * y;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s function_str num_params initial_y start_x end_x step_size seed\n", argv[0]);
        exit(1);
    }

    // 1. Parse arguments
    strncpy(g_data.function_str, argv[1], sizeof(g_data.function_str)-1);
    g_data.function_str[sizeof(g_data.function_str)-1] = '\0';
    g_data.num_params = atoi(argv[2]);
    g_data.initial_y = atof(argv[3]);
    g_data.start_x = atof(argv[4]);
    g_data.end_x = atof(argv[5]);
    g_data.step_size = atof(argv[6]);
    g_data.seed = (uint32_t)strtoul(argv[7], NULL, 10);

    // 2. Validate arguments
    if (g_data.num_params <= 0 || strcmp(g_data.function_str, "linear") != 0 || g_data.num_params != 3) {
        fprintf(stderr, "Error: This benchmark only supports function_str='linear' with num_params=3.\n");
        exit(1);
    }
    if (g_data.step_size <= 0.0) {
        fprintf(stderr, "Error: step_size must be positive.\n");
        exit(1);
    }
    if (g_data.end_x <= g_data.start_x) {
        fprintf(stderr, "Error: end_x must be greater than start_x.\n");
        exit(1);
    }

    // 3. Calculate number of steps
    g_data.num_steps = (long long)((g_data.end_x - g_data.start_x) / g_data.step_size);

    // 4. Seed the random number generator
    mt_seed(g_data.seed);

    // 5. Allocate and initialize data
    g_data.params = (double*)malloc(g_data.num_params * sizeof(double));
    if (g_data.params == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for parameters.\n");
        exit(1);
    }

    // Generate random parameters between -1.0 and 1.0
    for (int i = 0; i < g_data.num_params; i++) {
        g_data.params[i] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }

    g_data.result_accumulator = 0.0;
}

void run_computation() {
    double y = g_data.initial_y;
    double x = g_data.start_x;
    const double h = g_data.step_size;
    const long long steps = g_data.num_steps;

    double accumulator = 0.0;

    for (long long i = 0; i < steps; i++) {
        // Euler's method: y_{n+1} = y_n + h * f(x_n, y_n)
        double derivative = ode_function(x, y);
        y += h * derivative;
        x += h;

        // Accumulate results to prevent compiler optimization.
        accumulator += y;
    }
    g_data.result_accumulator = accumulator;
}

void cleanup() {
    free(g_data.params);
    g_data.params = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Handle potential overflow/underflow or NaN in the result.
    if (isinf(g_data.result_accumulator) || isnan(g_data.result_accumulator)){
        printf("%.6f\n", 0.0);
    } else {
        printf("%.6f\n", g_data.result_accumulator);
    }

    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
