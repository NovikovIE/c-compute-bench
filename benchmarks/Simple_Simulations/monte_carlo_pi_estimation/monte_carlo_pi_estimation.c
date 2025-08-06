#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

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

// Global variables to share data between setup and computation
long long num_samples;
long long circle_count; // Stores the result of the computation

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_samples> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atoll(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_samples <= 0) {
        fprintf(stderr, "Error: num_samples must be a positive integer.\n");
        exit(1);
    }

    // Seed the random number generator
    mt_seed(seed);

    // Initialize result variable
    circle_count = 0;
}

void run_computation() {
    long long in_circle = 0;
    const double max_rand_val = (double)UINT32_MAX;

    for (long long i = 0; i < num_samples; ++i) {
        // Generate a random point (x, y) in the unit square [0, 1] x [0, 1]
        double x = (double)mt_rand() / max_rand_val;
        double y = (double)mt_rand() / max_rand_val;

        // Check if the point is inside the unit circle's first quadrant
        if (x * x + y * y <= 1.0) {
            in_circle++;
        }
    }

    // Store the final count in the global variable to be used in main
    circle_count = in_circle;
}

void cleanup() {
    // No heap-allocated data to free for this specific benchmark
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    // Calculate elapsed time in seconds
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Calculate the final Pi estimate
    // The ratio of areas is (pi * r^2 / 4) / r^2 = pi / 4
    // pi / 4 ≈ circle_count / num_samples
    double pi_estimate = 4.0 * (double)circle_count / (double)num_samples;

    // Print the computational result to stdout
    printf("%.10f\n", pi_estimate);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
