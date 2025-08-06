#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) --- Do Not Modify ---
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

// --- Benchmark Specifics ---

// Global variables for benchmark data and parameters
int number_bits;
int precision_bits;
int num_words;

uint64_t *g_N_shifted;
uint64_t *g_x;
uint64_t *g_result;
uint64_t *g_temp1; // For division quotient
uint64_t *g_temp2; // For addition result
uint64_t *g_div_rem; // For division remainder
uint64_t *g_div_dvs; // For division shifted divisor

uint64_t final_checksum;

// --- Bignum Arithmetic Helper Functions ---

// Zero out a bignum
void bignum_zero(uint64_t *a, int n) {
    memset(a, 0, n * sizeof(uint64_t));
}

// Copy a bignum
void bignum_copy(uint64_t *dest, const uint64_t *src, int n) {
    memcpy(dest, src, n * sizeof(uint64_t));
}

// Compare two bignums: returns 1 if a > b, -1 if a < b, 0 if a == b
int bignum_cmp(const uint64_t *a, const uint64_t *b, int n) {
    for (int i = n - 1; i >= 0; i--) {
        if (a[i] > b[i]) return 1;
        if (a[i] < b[i]) return -1;
    }
    return 0;
}

// Get bit length of a bignum
int bignum_bit_length(const uint64_t *a, int n) {
    for (int i = n - 1; i >= 0; i--) {
        if (a[i] != 0) {
            return i * 64 + (64 - __builtin_clzll(a[i]));
        }
    }
    return 0;
}

// Set a specific bit k in a bignum
void bignum_set_bit(uint64_t *a, int k) {
    if (k < 0) return;
    int word_idx = k / 64;
    int bit_idx = k % 64;
    if (word_idx < num_words) {
        a[word_idx] |= (1ULL << bit_idx);
    }
}

// Add two bignums: r = a + b
void bignum_add(uint64_t *r, const uint64_t *a, const uint64_t *b, int n) {
    unsigned __int128 carry = 0;
    for (int i = 0; i < n; i++) {
        unsigned __int128 sum = (__int128)a[i] + b[i] + carry;
        r[i] = (uint64_t)sum;
        carry = sum >> 64;
    }
}

// Subtract two bignums: r = a - b (assumes a >= b)
void bignum_sub(uint64_t *r, const uint64_t *a, const uint64_t *b, int n) {
    unsigned __int128 borrow = 0;
    for (int i = 0; i < n; i++) {
        unsigned __int128 diff = (__int128)a[i] - b[i] - borrow;
        r[i] = (uint64_t)diff;
        borrow = (diff >> 127);
    }
}

// Shift right by 1 bit: a >>= 1
void bignum_rshift1(uint64_t *a, int n) {
    for (int i = 0; i < n - 1; i++) {
        a[i] = (a[i] >> 1) | (a[i + 1] << 63);
    }
    a[n - 1] >>= 1;
}

// Shift left by N bits
void bignum_lshift(uint64_t *r, const uint64_t *a, int n, int bits) {
    int word_shift = bits / 64;
    int bit_shift = bits % 64;
    
    bignum_zero(r, n);

    if (bit_shift == 0) {
        for (int i = 0; i < n - word_shift; i++) {
            r[i + word_shift] = a[i];
        }
    } else {
        for (int i = 0; i < n - word_shift; i++) {
            if (i + word_shift + 1 < n) {
                 r[i + word_shift + 1] |= (a[i] >> (64 - bit_shift));
            }
            if (i + word_shift < n) {
                r[i + word_shift] |= (a[i] << bit_shift);
            }
        }
    }
}

// Simple, slow restoring division: u / v. Quotient in q.
// Uses global buffers g_div_rem and g_div_dvs for performance.
void bignum_div(uint64_t *q, const uint64_t *u, const uint64_t *v, int n) {
    bignum_zero(q, n);
    bignum_copy(g_div_rem, u, n);
    
    int u_bits = bignum_bit_length(u, n);
    int v_bits = bignum_bit_length(v, n);

    if (v_bits == 0) return; // Division by zero
    int shift = u_bits - v_bits;
    if (shift < 0) return; // u < v, quotient is 0

    bignum_lshift(g_div_dvs, v, n, shift);

    for (int i = 0; i <= shift; i++) {
        if (bignum_cmp(g_div_rem, g_div_dvs, n) >= 0) {
            bignum_sub(g_div_rem, g_div_rem, g_div_dvs, n);
            bignum_set_bit(q, shift - i);
        }
        bignum_rshift1(g_div_dvs, n);
    }
}

// ----- Benchmark Functions -----

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <number_bits> <precision_bits> <seed>\n", argv[0]);
        exit(1);
    }

    number_bits = atoi(argv[1]);
    precision_bits = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    
    mt_seed(seed);

    int working_bits = number_bits + precision_bits;
    num_words = (working_bits + 63) / 64 + 1; // +1 for safety

    // Allocate memory
    g_N_shifted = (uint64_t*)malloc(num_words * sizeof(uint64_t));
    g_x = (uint64_t*)malloc(num_words * sizeof(uint64_t));
    g_result = (uint64_t*)malloc(num_words * sizeof(uint64_t));
    g_temp1 = (uint64_t*)malloc(num_words * sizeof(uint64_t));
    g_temp2 = (uint64_t*)malloc(num_words * sizeof(uint64_t));
    g_div_rem = (uint64_t*)malloc(num_words * sizeof(uint64_t));
    g_div_dvs = (uint64_t*)malloc(num_words * sizeof(uint64_t));

    uint64_t *N_temp = (uint64_t*)malloc(num_words * sizeof(uint64_t));

    // Check allocations
    if (!g_N_shifted || !g_x || !g_result || !g_temp1 || !g_temp2 || !g_div_rem || !g_div_dvs || !N_temp) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate a random 'number_bits' integer N
    bignum_zero(N_temp, num_words);
    int top_word = (number_bits - 1) / 64;
    for (int i = 0; i <= top_word; i++) {
        N_temp[i] = ((uint64_t)mt_rand() << 32) | mt_rand();
    }
    // Clear bits above number_bits to be precise
    int excess_bits = (top_word + 1) * 64 - number_bits;
    if (excess_bits > 0 && excess_bits < 64) {
        N_temp[top_word] >>= excess_bits;
    }
    // Ensure the top bit is set for full length
    bignum_set_bit(N_temp, number_bits - 1);

    // Shift N by precision_bits to make it a fixed-point number for calculation
    bignum_lshift(g_N_shifted, N_temp, num_words, precision_bits);
    
    free(N_temp);

    // Initial guess x_0 = 2^((number_bits + precision_bits) / 2)
    bignum_zero(g_x, num_words);
    bignum_set_bit(g_x, working_bits / 2);
}

void run_computation() {
    int working_bits = number_bits + precision_bits;
    // Number of iterations needed is related to log2 of the precision.
    // We add a few extra for convergence assurance.
    int max_iter = (int)(log2(working_bits) + 3);

    for (int i = 0; i < max_iter; i++) {
        // temp1 = N_shifted / x
        bignum_div(g_temp1, g_N_shifted, g_x, num_words);

        // temp2 = x + temp1
        bignum_add(g_temp2, g_x, g_temp1, num_words);

        // result = temp2 / 2
        bignum_copy(g_result, g_temp2, num_words);
        bignum_rshift1(g_result, num_words);

        // if x is unchanged, we have converged
        if (bignum_cmp(g_result, g_x, num_words) == 0) {
            break;
        }

        // x = result
        bignum_copy(g_x, g_result, num_words);
    }

    // Calculate a checksum to prevent dead code elimination
    final_checksum = 0;
    for (int i = 0; i < num_words; i++) {
        final_checksum ^= g_x[i];
    }
}

void cleanup() {
    free(g_N_shifted);
    free(g_x);
    free(g_result);
    free(g_temp1);
    free(g_temp2);
    free(g_div_rem);
    free(g_div_dvs);
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
    printf("%llu\n", (unsigned long long)final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
