#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// Mersenne Twister (MT19937) generator - DO NOT MODIFY
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
// End of Mersenne Twister

// Global structure to hold benchmark data
typedef struct {
    int num_particles;
    int num_time_steps;
    double* pos_x;
    double* pos_y;
    double* pos_z;
    double final_result;
} BenchmarkData;

BenchmarkData g_data;

// Helper to generate a random double in [-0.5, 0.5]
static inline double rand_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) - 0.5;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_particles> <num_time_steps> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_particles = atoi(argv[1]);
    g_data.num_time_steps = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);
    
    if (g_data.num_particles <= 0 || g_data.num_time_steps <= 0) {
        fprintf(stderr, "FATAL: Parameters must be positive integers.\n");
        exit(1);
    }
    
    mt_seed(seed);

    // Allocate memory for particle positions, initialized to zero
    g_data.pos_x = (double*)calloc(g_data.num_particles, sizeof(double));
    g_data.pos_y = (double*)calloc(g_data.num_particles, sizeof(double));
    g_data.pos_z = (double*)calloc(g_data.num_particles, sizeof(double));

    if (!g_data.pos_x || !g_data.pos_y || !g_data.pos_z) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        free(g_data.pos_x);
        free(g_data.pos_y);
        free(g_data.pos_z);
        exit(1);
    }
    
    g_data.final_result = 0.0;
}

void run_computation() {
    // Simulate the random walk in 3D
    for (int t = 0; t < g_data.num_time_steps; ++t) {
        for (int i = 0; i < g_data.num_particles; ++i) {
            g_data.pos_x[i] += rand_double();
            g_data.pos_y[i] += rand_double();
            g_data.pos_z[i] += rand_double();
        }
    }

    // Calculate the final result to prevent dead code elimination.
    // The result is the sum of the squared distances from the origin.
    double sum_sq_dist = 0.0;
    for (int i = 0; i < g_data.num_particles; ++i) {
        double px = g_data.pos_x[i];
        double py = g_data.pos_y[i];
        double pz = g_data.pos_z[i];
        sum_sq_dist += px*px + py*py + pz*pz;
    }
    g_data.final_result = sum_sq_dist;
}

void cleanup() {
    free(g_data.pos_x);
    free(g_data.pos_y);
    free(g_data.pos_z);
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
