#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator ---
// Do Not Modify
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


// --- BigNum Implementation ---

// Represents a large unsigned integer using a base-256 representation.
// digits[0] is the least significant byte.
typedef struct {
    unsigned char *digits;
    int size;      // Number of bytes currently used.
    int capacity;  // Allocated number of bytes.
} BigNum;

// Global variables for the benchmark
int g_num_digits;
int g_modulus_digits;
BigNum *g_number;
BigNum *g_modulus;
BigNum *g_result;
unsigned long long g_checksum;

// Function Prototypes for BigNum operations
BigNum* bignum_create(int capacity);
void bignum_free(BigNum* bn);
void bignum_normalize(BigNum *bn);
int bignum_compare(const BigNum *a, const BigNum *b);
void bignum_subtract_inplace(BigNum *a, const BigNum *b);
void bignum_mod(const BigNum *dividend, const BigNum *modulus, BigNum *result);


// Helper to create and initialize a BigNum
BigNum* bignum_create(int capacity) {
    if (capacity <= 0) capacity = 1;
    BigNum* bn = (BigNum*)malloc(sizeof(BigNum));
    if (!bn) exit(1);
    bn->digits = (unsigned char*)calloc(capacity, sizeof(unsigned char));
    if (!bn->digits) exit(1);
    bn->size = 1;
    bn->capacity = capacity;
    return bn;
}

// Helper to free a BigNum's memory
void bignum_free(BigNum* bn) {
    if (bn) {
        free(bn->digits);
        free(bn);
    }
}

// Trim leading zero bytes to maintain a canonical representation
void bignum_normalize(BigNum *bn) {
    while (bn->size > 1 && bn->digits[bn->size - 1] == 0) {
        bn->size--;
    }
}

// Compare two BigNums. Returns:
//  1 if a > b
// -1 if a < b
//  0 if a == b
int bignum_compare(const BigNum *a, const BigNum *b) {
    if (a->size > b->size) return 1;
    if (a->size < b->size) return -1;
    for (int i = a->size - 1; i >= 0; i--) {
        if (a->digits[i] > b->digits[i]) return 1;
        if (a->digits[i] < b->digits[i]) return -1;
    }
    return 0;
}

// Subtracts b from a (a = a - b). Assumes a >= b.
void bignum_subtract_inplace(BigNum *a, const BigNum *b) {
    int borrow = 0;
    for (int i = 0; i < b->size; i++) {
        int diff = a->digits[i] - b->digits[i] - borrow;
        if (diff < 0) {
            diff += 256;
            borrow = 1;
        } else {
            borrow = 0;
        }
        a->digits[i] = (unsigned char)diff;
    }

    for (int i = b->size; i < a->size && borrow > 0; i++) {
        int diff = a->digits[i] - borrow;
        if (diff < 0) {
            diff += 256;
            borrow = 1;
        } else {
            borrow = 0;
        }
        a->digits[i] = (unsigned char)diff;
    }
    bignum_normalize(a);
}

// Computes result = dividend % modulus using schoolbook long division.
void bignum_mod(const BigNum *dividend, const BigNum *modulus, BigNum *result) {
    // Create a temporary BigNum for the current remainder.
    // It needs to be one byte larger than the modulus to hold intermediate values.
    BigNum *rem = bignum_create(modulus->size + 1);
    
    // Iterate through dividend's digits from most significant to least
    for (int i = dividend->size - 1; i >= 0; i--) {
        // "Bring down" the next digit. This is equivalent to: rem = rem * 256 + digit
        // Shift rem left by one byte
        memmove(&rem->digits[1], &rem->digits[0], rem->size);
        if (rem->size < rem->capacity) {
            rem->size++;
        }
        rem->digits[0] = dividend->digits[i];
        bignum_normalize(rem);

        // Subtract modulus from remainder until remainder < modulus
        while (bignum_compare(rem, modulus) >= 0) {
            bignum_subtract_inplace(rem, modulus);
        }
    }

    // The final remainder is the result.
    memset(result->digits, 0, result->capacity);
    memcpy(result->digits, rem->digits, rem->size);
    result->size = rem->size;
    bignum_normalize(result);

    bignum_free(rem);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_digits> <modulus_digits> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_digits = atoi(argv[1]);
    g_modulus_digits = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if(g_num_digits <= 0 || g_modulus_digits <= 0 || g_num_digits < g_modulus_digits) {
        fprintf(stderr, "FATAL: Invalid arguments. num_digits and modulus_digits must be > 0 and num_digits >= modulus_digits.\n");
        exit(1);
    }

    mt_seed(seed);

    g_number = bignum_create(g_num_digits);
    g_modulus = bignum_create(g_modulus_digits);
    g_result = bignum_create(g_modulus_digits);

    // Generate random dividend
    g_number->size = g_num_digits;
    for (int i = 0; i < g_num_digits; i++) {
        g_number->digits[i] = mt_rand() & 0xFF;
    }

    // Generate random modulus
    g_modulus->size = g_modulus_digits;
    for (int i = 0; i < g_modulus_digits; i++) {
        g_modulus->digits[i] = mt_rand() & 0xFF;
    }
    // Ensure modulus's most significant byte is not zero to avoid trivial cases
    if (g_modulus_digits > 0) {
        g_modulus->digits[g_modulus_digits-1] = (mt_rand() % 255) + 1;
    }
    bignum_normalize(g_number);
    bignum_normalize(g_modulus);

    g_checksum = 0;
}

void run_computation() {
    bignum_mod(g_number, g_modulus, g_result);

    // Calculate a checksum of the result to prevent dead code elimination
    for (int i = 0; i < g_result->size; i++) {
        g_checksum = (g_checksum * 31 + g_result->digits[i]);
    }
}

void cleanup() {
    bignum_free(g_number);
    bignum_free(g_modulus);
    bignum_free(g_result);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    double time_taken;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%llu\n", g_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
