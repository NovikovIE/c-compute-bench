#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- Arbitrary Precision Arithmetic Structures ---
typedef struct {
    uint32_t *limbs;  // Array of 32-bit words
    int capacity;     // Allocated number of limbs
    int size;         // Used number of limbs
} BigInt;

typedef struct {
    BigInt *num;
    BigInt *den;
} Rational;

// --- Global Variables ---
// Parameters
static int NUMERATOR_BITS;
static int DENOMINATOR_BITS;
static int NUM_ITERATIONS;

// Data
static Rational *r_current;
static Rational *r_addend;

// Temporaries for computation
static BigInt *temp_ad, *temp_bc, *temp_sum, *temp_bd;

// Result accumulator
static unsigned long long final_result;

// --- Helper Functions ---
BigInt* bigint_create(int capacity) {
    BigInt* b = (BigInt*)malloc(sizeof(BigInt));
    if (!b) return NULL;
    b->limbs = (uint32_t*)malloc(capacity * sizeof(uint32_t));
    if (!b->limbs) {
        free(b);
        return NULL;
    }
    memset(b->limbs, 0, capacity * sizeof(uint32_t));
    b->capacity = capacity;
    b->size = 1;
    return b;
}

void bigint_destroy(BigInt *b) {
    if (b) {
        free(b->limbs);
        free(b);
    }
}

void bigint_normalize(BigInt* b) {
    int s = b->size;
    while (s > 1 && b->limbs[s - 1] == 0) {
        s--;
    }
    b->size = s;
}

void bigint_add(BigInt* res, const BigInt* a, const BigInt* b) {
    int max_size = a->size > b->size ? a->size : b->size;
    if (res->capacity < max_size + 1) {
        // This case should not be hit with proper pre-allocation
        return;
    }
    uint64_t carry = 0;
    for (int i = 0; i < max_size; ++i) {
        uint64_t val_a = (i < a->size) ? a->limbs[i] : 0;
        uint64_t val_b = (i < b->size) ? b->limbs[i] : 0;
        uint64_t sum = val_a + val_b + carry;
        res->limbs[i] = (uint32_t)sum;
        carry = sum >> 32;
    }
    if (carry > 0) {
        res->limbs[max_size] = (uint32_t)carry;
        res->size = max_size + 1;
    } else {
        res->size = max_size;
    }
    bigint_normalize(res);
}

void bigint_mul(BigInt* res, const BigInt* a, const BigInt* b) {
    if (res->capacity < a->size + b->size) {
        // This case should not be hit with proper pre-allocation
        return;
    }
    memset(res->limbs, 0, res->capacity * sizeof(uint32_t));
    for (int i = 0; i < b->size; ++i) {
        uint64_t carry = 0;
        for (int j = 0; j < a->size; ++j) {
            uint64_t prod = (uint64_t)a->limbs[j] * b->limbs[i] + res->limbs[i + j] + carry;
            res->limbs[i + j] = (uint32_t)prod;
            carry = prod >> 32;
        }
        res->limbs[i + a->size] += (uint32_t)carry;
    }
    res->size = a->size + b->size;
    bigint_normalize(res);
}

void bigint_copy(BigInt* dest, const BigInt* src) {
    if (dest->capacity < src->size) return; // Should not happen
    memcpy(dest->limbs, src->limbs, src->size * sizeof(uint32_t));
    dest->size = src->size;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s numerator_bits denominator_bits num_iterations seed\n", argv[0]);
        exit(1);
    }
    NUMERATOR_BITS = atoi(argv[1]);
    DENOMINATOR_BITS = atoi(argv[2]);
    NUM_ITERATIONS = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    int num_limbs_num = (NUMERATOR_BITS + 31) / 32;
    int num_limbs_den = (DENOMINATOR_BITS + 31) / 32;

    // Allocate Rational numbers
    r_current = (Rational*)malloc(sizeof(Rational));
    r_addend = (Rational*)malloc(sizeof(Rational));

    // Allocate BigInts. We add a generous buffer for growth during addition.
    int cap_num = 2 * (num_limbs_num + num_limbs_den);
    int cap_den = 2 * (num_limbs_den + num_limbs_den);
    r_current->num = bigint_create(cap_num);
    r_current->den = bigint_create(cap_den);
    r_addend->num = bigint_create(num_limbs_num);
    r_addend->den = bigint_create(num_limbs_den);

    // Initialize BigInts with random data
    for (int i = 0; i < num_limbs_num; ++i) r_current->num->limbs[i] = mt_rand();
    for (int i = 0; i < num_limbs_den; ++i) r_current->den->limbs[i] = mt_rand();
    for (int i = 0; i < num_limbs_num; ++i) r_addend->num->limbs[i] = mt_rand();
    for (int i = 0; i < num_limbs_den; ++i) r_addend->den->limbs[i] = mt_rand();
    
    r_current->num->size = num_limbs_num;
    r_current->den->size = num_limbs_den;
    r_addend->num->size = num_limbs_num;
    r_addend->den->size = num_limbs_den;
    
    // Ensure top limbs and denominators are not zero
    r_current->num->limbs[num_limbs_num -1] |= 1;
    r_current->den->limbs[num_limbs_den - 1] |= 1;
    r_addend->num->limbs[num_limbs_num - 1] |= 1;
    r_addend->den->limbs[num_limbs_den - 1] |= 1;

    bigint_normalize(r_current->num);
    bigint_normalize(r_current->den);
    bigint_normalize(r_addend->num);
    bigint_normalize(r_addend->den);

    // Allocate temporary BigInts for computation. Sizes are based on operation results.
    // ad, bc: size(a)+size(d). sum: size(ad)+1. bd: size(b)+size(d)
    temp_ad = bigint_create(cap_num);
    temp_bc = bigint_create(cap_num);
    temp_sum = bigint_create(cap_num + 1);
    temp_bd = bigint_create(cap_den);
}

void run_computation() {
    final_result = 0;
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        // Compute r_current = r_current + r_addend
        // Formula: a/b + c/d = (ad + bc) / bd

        // ad
        bigint_mul(temp_ad, r_current->num, r_addend->den);
        // bc
        bigint_mul(temp_bc, r_addend->num, r_current->den);
        // ad + bc
        bigint_add(temp_sum, temp_ad, temp_bc);
        // bd
        bigint_mul(temp_bd, r_current->den, r_addend->den);

        // Copy result back to r_current
        bigint_copy(r_current->num, temp_sum);
        bigint_copy(r_current->den, temp_bd);

        // Accumulate a value to prevent dead code elimination
        if (r_current->num->size > 0) {
            final_result += r_current->num->limbs[r_current->num->size - 1];
        }
    }
}

void cleanup() {
    bigint_destroy(r_current->num);
    bigint_destroy(r_current->den);
    free(r_current);

    bigint_destroy(r_addend->num);
    bigint_destroy(r_addend->den);
    free(r_addend);

    bigint_destroy(temp_ad);
    bigint_destroy(temp_bc);
    bigint_destroy(temp_sum);
    bigint_destroy(temp_bd);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print the final result to stdout
    printf("%llu\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
