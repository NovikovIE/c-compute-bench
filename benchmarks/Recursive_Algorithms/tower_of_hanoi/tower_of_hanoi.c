#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (Provided) ---
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

// --- Global state for the benchmark ---
int g_num_disks;
unsigned long long g_total_moves;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_disks> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_disks = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_num_disks < 1 || g_num_disks > 63) {
        fprintf(stderr, "Number of disks must be between 1 and 63.\n");
        exit(1);
    }

    // Seed the random number generator (required by spec, though not used by this algorithm)
    mt_seed(seed);

    // Initialize benchmark state
    g_total_moves = 0;
    
    // Note: This benchmark does not require heap allocation as its primary
    // performance characteristic is call stack depth and recursive overhead.
}

// Recursive helper function for Tower of Hanoi
static void tower_of_hanoi_recursive(int n) {
    if (n <= 0) {
        return;
    }
    // Move n-1 disks from source to auxiliary (recursive call)
    tower_of_hanoi_recursive(n - 1);
    // Move the nth disk from source to destination (counted as one move)
    g_total_moves++;
    // Move n-1 disks from auxiliary to destination (recursive call)
    tower_of_hanoi_recursive(n - 1);
}

void run_computation() {
    tower_of_hanoi_recursive(g_num_disks);
}

void cleanup() {
    // No heap memory was allocated in setup_benchmark, so nothing to free.
    // This function is present to conform to the required structure.
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout.
    // The number of moves is 2^n - 1. We use unsigned long long to avoid overflow.
    printf("%llu\n", g_total_moves);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
