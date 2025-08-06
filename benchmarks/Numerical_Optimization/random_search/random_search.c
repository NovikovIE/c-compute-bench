#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <float.h>
#include <string.h>

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

// Benchmark-specific data
#define SEARCH_MIN -5.12
#define SEARCH_MAX 5.12

// Using a struct for clarity and organization
typedef struct {
    int num_dimensions;
    int num_samples;
    double* best_solution;
    double best_fitness;
    double result_checksum;
} BenchmarkData;

// Global pointer to the benchmark data
static BenchmarkData g_data;

// Objective function to minimize (Sphere function)
double sphere_function(const double* vector, int dimensions) {
    double sum = 0.0;
    for (int i = 0; i < dimensions; ++i) {
        sum += vector[i] * vector[i];
    }
    return sum;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_dimensions> <num_samples> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_dimensions = atoi(argv[1]);
    g_data.num_samples = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_dimensions <= 0 || g_data.num_samples <= 0) {
        fprintf(stderr, "FATAL: num_dimensions and num_samples must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.best_solution = (double*)malloc(g_data.num_dimensions * sizeof(double));
    if (g_data.best_solution == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for best_solution.\n");
        exit(1);
    }

    // Initialize with a known high value
    g_data.best_fitness = DBL_MAX;
    g_data.result_checksum = 0.0;
}

void run_computation() {
    double* candidate = (double*)malloc(g_data.num_dimensions * sizeof(double));
    if (candidate == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for candidate solution.\n");
        exit(1);
    }

    const double range = SEARCH_MAX - SEARCH_MIN;
    const double u32_max_inv = 1.0 / 4294967295.0; // 1.0 / (2^32 - 1)

    for (int i = 0; i < g_data.num_samples; ++i) {
        // Generate a random candidate solution
        for (int j = 0; j < g_data.num_dimensions; ++j) {
            double random_scaled = mt_rand() * u32_max_inv;
            candidate[j] = SEARCH_MIN + random_scaled * range;
        }

        // Evaluate the candidate
        double current_fitness = sphere_function(candidate, g_data.num_dimensions);

        // If it's a better solution, update the best known
        if (current_fitness < g_data.best_fitness) {
            g_data.best_fitness = current_fitness;
            memcpy(g_data.best_solution, candidate, g_data.num_dimensions * sizeof(double));
        }
    }

    free(candidate);

    // Calculate a checksum to prevent dead code elimination of the final result
    g_data.result_checksum = 0.0;
    for (int i = 0; i < g_data.num_dimensions; ++i) {
        g_data.result_checksum += g_data.best_solution[i];
    }
}

void cleanup() {
    if (g_data.best_solution != NULL) {
        free(g_data.best_solution);
        g_data.best_solution = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", g_data.result_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
