#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>


// START of Mersenne Twister (Do Not Modify)
// For a complete explanation of the Mersenne Twister, see:
// https://en.wikipedia.org/wiki/Mersenne_Twister

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
// END of Mersenne Twister


// --- Global Benchmark Data ---

// Input parameter: number of pairs of parentheses.
static int g_num_pairs;

// Output: Total count of valid sequences.
static unsigned long long g_count;

// --- Forward Declaration for Recursive Helper ---
void count_sequences_recursive(int open_count, int close_count);


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_pairs> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_pairs = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    // The number of valid sequences (Catalan numbers) grows very fast.
    // C(33) is the last one that fits in unsigned long long.
    if (g_num_pairs < 0 || g_num_pairs > 33) {
        fprintf(stderr, "Error: num_pairs must be between 0 and 33 inclusive.\n");
        exit(1);
    }

    // Seed the random number generator as required, though it is not used in this specific computation.
    mt_seed(seed);

    // Initialize the result counter.
    g_count = 0;
}

void run_computation() {
    // The number of well-formed parentheses strings is given by the Catalan numbers.
    // This benchmark computes them via recursive enumeration.
    if (g_num_pairs == 0) {
        g_count = 1; // There is one way for 0 pairs: the empty string.
    } else {
        count_sequences_recursive(0, 0);
    }
}

void cleanup() {
    // No heap memory was allocated for this benchmark.
}

// --- Recursive Helper Function ---

void count_sequences_recursive(int open_count, int close_count) {
    // Base case: A full, valid sequence has been formed.
    if (open_count == g_num_pairs && close_count == g_num_pairs) {
        g_count++;
        return;
    }

    // Recursive step 1: If we can add an opening parenthesis, do so.
    if (open_count < g_num_pairs) {
        count_sequences_recursive(open_count + 1, close_count);
    }

    // Recursive step 2: If we can add a closing parenthesis without invalidating the sequence,
    // (i.e., number of closed < number of open), do so.
    if (close_count < open_count) {
        count_sequences_recursive(open_count, close_count + 1);
    }
}


// --- Main Execution Logic ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (total count) to stdout.
    printf("%llu\n", g_count);

    // Print the execution time to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
