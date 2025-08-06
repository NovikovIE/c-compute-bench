#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- END Mersenne Twister ---

// --- Benchmark Globals ---
typedef struct {
    int dividend_len;
    int divisor_len;
    int quotient_len;
    uint8_t *dividend;
    uint8_t *divisor;
    uint8_t *quotient;
    long long final_result; // Use long long for the accumulator
} BenchmarkData;

BenchmarkData g_data;

// --- Bignum Helper Functions ---

// Compares two big-endian numbers a and b of length n.
// Returns 1 if a < b, 0 otherwise.
static inline int is_less(const uint8_t* a, const uint8_t* b, int n) {
    for (int i = 0; i < n; i++) {
        if (a[i] < b[i]) return 1;
        if (a[i] > b[i]) return 0;
    }
    return 0; // a == b
}

// Subtracts big-endian number b from a (a -= b), both of length n.
// Assumes a >= b.
static inline void subtract(uint8_t* a, const uint8_t* b, int n) {
    int borrow = 0;
    for (int i = n - 1; i >= 0; i--) {
        int diff = (int)a[i] - b[i] - borrow;
        if (diff < 0) {
            diff += 256;
            borrow = 1;
        } else {
            borrow = 0;
        }
        a[i] = (uint8_t)diff;
    }
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <dividend_digits> <divisor_digits> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.dividend_len = atoi(argv[1]);
    g_data.divisor_len = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.dividend_len <= 0 || g_data.divisor_len <= 0) {
        fprintf(stderr, "FATAL: Digit counts must be positive.\n");
        exit(1);
    }
    if (g_data.divisor_len > g_data.dividend_len) {
        fprintf(stderr, "FATAL: Divisor cannot be longer than dividend.\n");
        exit(1);
    }
    
    mt_seed(seed);

    g_data.dividend = (uint8_t*)malloc(g_data.dividend_len * sizeof(uint8_t));
    g_data.divisor = (uint8_t*)malloc(g_data.divisor_len * sizeof(uint8_t));
    
    g_data.quotient_len = g_data.dividend_len - g_data.divisor_len + 1;
    g_data.quotient = (uint8_t*)malloc(g_data.quotient_len * sizeof(uint8_t));

    if (!g_data.dividend || !g_data.divisor || !g_data.quotient) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    
    // Generate dividend
    for (int i = 0; i < g_data.dividend_len; i++) {
        g_data.dividend[i] = mt_rand() % 256;
    }

    // Generate divisor
    for (int i = 0; i < g_data.divisor_len; i++) {
        g_data.divisor[i] = mt_rand() % 256;
    }
    
    // Ensure divisor's most significant digit is not zero
    // to avoid division by zero and simplify the algorithm.
    if (g_data.divisor_len > 0) {
        g_data.divisor[0] = (mt_rand() % 255) + 1;
    }

    // Ensure the most significant part of the dividend is larger than the divisor
    // to produce a non-zero quotient. This is not strictly necessary for correctness
    // but makes the benchmark more consistent.
    while (is_less(g_data.dividend, g_data.divisor, g_data.divisor_len)) {
        g_data.dividend[0] = (mt_rand() % 255) + 1;
    }
    
    g_data.final_result = 0;
}

void run_computation() {
    // This implements a basic but computationally intensive long division by repeated subtraction.
    // We compute Q by repeatedly subtracting the Divisor (shifted) from the Dividend.

    uint8_t* temp_dividend = (uint8_t*)malloc(g_data.dividend_len * sizeof(uint8_t));
    if (!temp_dividend) { fprintf(stderr, "FATAL: Malloc failed.\n"); exit(1); }
    memcpy(temp_dividend, g_data.dividend, g_data.dividend_len);

    memset(g_data.quotient, 0, g_data.quotient_len * sizeof(uint8_t));

    uint8_t* shifted_divisor = (uint8_t*)calloc(g_data.dividend_len, sizeof(uint8_t));
    if (!shifted_divisor) { fprintf(stderr, "FATAL: Malloc failed.\n"); free(temp_dividend); exit(1); }

    for (int i = 0; i <= g_data.dividend_len - g_data.divisor_len; i++) {
        // Create a representation of the divisor shifted to the current position.
        // E.g., for 9876/54, first we check 9876 vs 5400.
        memcpy(shifted_divisor + i, g_data.divisor, g_data.divisor_len);
        
        // Repeatedly subtract the shifted divisor from the temporary dividend.
        while (!is_less(temp_dividend, shifted_divisor, g_data.dividend_len)) {
            subtract(temp_dividend, shifted_divisor, g_data.dividend_len);
            g_data.quotient[i]++;
        }
        // Reset the just-used part of shifted_divisor to 0 for the next iteration.
        memset(shifted_divisor + i, 0, g_data.divisor_len);
    }

    free(shifted_divisor);
    free(temp_dividend);

    // Accumulate the result to prevent dead code elimination
    g_data.final_result = 0;
    for (int i = 0; i < g_data.quotient_len; i++) {
        g_data.final_result += g_data.quotient[i];
    }
}

void cleanup() {
    free(g_data.dividend);
    free(g_data.divisor);
    free(g_data.quotient);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
