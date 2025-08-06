#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND GLOBALS ---
int n_term;
// Use unsigned long long to avoid overflow for Fibonacci numbers up to F(93)
unsigned long long* fib_table;
unsigned long long final_result;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <n_term> <seed>\n", argv[0]);
        exit(1);
    }

    n_term = atoi(argv[1]);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);

    // The Nth Fibonacci number for N > 93 exceeds the capacity of unsigned long long.
    if (n_term < 0 || n_term > 93) {
        fprintf(stderr, "Error: n_term must be between 0 and 93 for this benchmark.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the DP table (tabulation method)
    // We need space for n_term + 1 numbers (from 0 to n_term)
    fib_table = (unsigned long long*)malloc((n_term + 1) * sizeof(unsigned long long));
    if (fib_table == NULL) {
        fprintf(stderr, "Failed to allocate memory for the Fibonacci table.\n");
        exit(1);
    }

    final_result = 0;
}

void run_computation() {
    // This constant is tuned to make the computation take roughly 1 second
    // for the chosen n_term value.
    const int REPEATS = 15000000;

    for (int r = 0; r < REPEATS; ++r) {
        unsigned long long current_fib_n;
        
        if (n_term == 0) {
            current_fib_n = 0;
        } else if (n_term == 1) {
            current_fib_n = 1;
        } else {
            // Standard bottom-up DP (tabulation) for Fibonacci
            fib_table[0] = 0;
            fib_table[1] = 1;
            for (int i = 2; i <= n_term; ++i) {
                fib_table[i] = fib_table[i - 1] + fib_table[i - 2];
            }
            current_fib_n = fib_table[n_term];
        }

        // Accumulate result using XOR to prevent both overflow and dead code elimination
        // by the compiler. The result will be 0 if REPEATS is even, and current_fib_n if odd.
        final_result ^= current_fib_n;
    }
}

void cleanup() {
    free(fib_table);
    fib_table = NULL;
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    // The result is an accumulation to prevent dead code elimination
    printf("%llu\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
