#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// For M_PI on some compilers
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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
// --- End of MT19937 --- 

// --- Benchmark Globals and Setup ---

typedef struct {
    int num_dimensions;
    int num_iterations;
    double *initial_point;
    double *current_point;
    double final_result;
} BenchmarkState;

static BenchmarkState state;

static const double RAST_A = 10.0;
static const double RAST_DOMAIN_MIN = -5.12;
static const double RAST_DOMAIN_MAX = 5.12;

// Helper to generate a random double in a given range
double random_double(double min, double max) {
    return min + (max - min) * (mt_rand() / (double)UINT32_MAX);
}

// Objective function to minimize (Rastrigin function)
static inline double objective_function(const double* x) {
    double sum = RAST_A * state.num_dimensions;
    for (int i = 0; i < state.num_dimensions; ++i) {
        sum += x[i] * x[i] - RAST_A * cos(2.0 * M_PI * x[i]);
    }
    return sum;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_dimensions> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    state.num_dimensions = atoi(argv[1]);
    state.num_iterations = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (state.num_dimensions <= 0 || state.num_iterations <= 0) {
        fprintf(stderr, "FATAL: Dimensions and iterations must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    state.initial_point = (double*)malloc(state.num_dimensions * sizeof(double));
    state.current_point = (double*)malloc(state.num_dimensions * sizeof(double));
    if (!state.initial_point || !state.current_point) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < state.num_dimensions; ++i) {
        state.initial_point[i] = random_double(RAST_DOMAIN_MIN, RAST_DOMAIN_MAX);
    }
}

void cleanup() {
    free(state.initial_point);
    free(state.current_point);
}

// --- Core Computation ---

void run_computation() {
    size_t vec_size = state.num_dimensions * sizeof(double);
    memcpy(state.current_point, state.initial_point, vec_size);

    double delta = 1.0; // Initial step size
    double current_value = objective_function(state.current_point);

    for (int i = 0; i < state.num_iterations; ++i) {
        int improved_in_iteration = 0;

        // Perform one pass of coordinate-wise search
        for (int j = 0; j < state.num_dimensions; ++j) {
            double original_coord = state.current_point[j];

            // Try positive step
            state.current_point[j] = original_coord + delta;
            double new_value = objective_function(state.current_point);

            if (new_value < current_value) {
                current_value = new_value;
                improved_in_iteration = 1;
                continue; // Keep the change and move to the next dimension
            }

            // Positive step was not better, try negative step
            state.current_point[j] = original_coord - delta;
            new_value = objective_function(state.current_point);

            if (new_value < current_value) {
                current_value = new_value;
                improved_in_iteration = 1;
                continue; // Keep the change
            }
            
            // Neither direction was an improvement, revert to original
            state.current_point[j] = original_coord;
        }

        if (!improved_in_iteration) {
            delta *= 0.5; // Reduce step size
            if (delta < 1e-9) { // Convergence criteria
                break;
            }
        }
    }

    state.final_result = current_value;
}


// --- Main and Timing ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", state.final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
