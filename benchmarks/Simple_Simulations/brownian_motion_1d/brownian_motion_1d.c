#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (Do Not Modify) ---
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
// --- End of Mersenne Twister ---

// Benchmark-specific global variables
int num_particles;
int num_steps;
double* particles_pos;
double final_result; // To store the accumulated result and prevent DCE

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_particles> <num_steps> <seed>\n", argv[0]);
        exit(1);
    }

    num_particles = atoi(argv[1]);
    num_steps = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_particles <= 0 || num_steps <= 0) {
        fprintf(stderr, "ERROR: num_particles and num_steps must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    particles_pos = (double*)malloc(num_particles * sizeof(double));
    if (particles_pos == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for particle positions.\n");
        exit(1);
    }

    // Initialize all particles at position 0.0
    for (int i = 0; i < num_particles; ++i) {
        particles_pos[i] = 0.0;
    }
}

void run_computation() {
    for (int s = 0; s < num_steps; ++s) {
        for (int p = 0; p < num_particles; ++p) {
            // Simulate a single step of a 1D random walk.
            // A step can be -1, 0, or +1.
            int move = (mt_rand() % 3) - 1;
            particles_pos[p] += (double)move;
        }
    }

    // To prevent dead-code elimination, calculate a final value.
    // We'll compute the sum of the final absolute positions of all particles.
    double total_displacement = 0.0;
    for (int p = 0; p < num_particles; ++p) {
        total_displacement += fabs(particles_pos[p]);
    }
    final_result = total_displacement;
}

void cleanup() {
    free(particles_pos);
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

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
