/*
 * rsa_encrypt: A CPU benchmark simulating RSA encryption.
 *
 * This program benchmarks the performance of modular exponentiation with large integers,
 * which is the core mathematical operation in RSA encryption.
 *
 * The program generates a random large number 'n' (the modulus) and a fixed public
 * exponent 'e'. It then encrypts a block of random data by repeatedly computing
 * c = m^e mod n, where 'm' is a chunk of the input data.
 *
 * A minimal-dependency bignum library is implemented directly in this file to handle
 * arithmetic with numbers larger than standard integer types.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// --- Mersenne Twister (MT19937) PRNG --- VERBATIM, DO NOT MODIFY ---
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
                fprintf(stderr, "FATAL: Mersenne Twister not seeded.\n");
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
// --- End of MT19937 ---

// --- Minimal Bignum Library ---

#define BN_MAX_BITS 4096
#define BN_LIMB_BITS 32
#define BN_MAX_LIMBS (BN_MAX_BITS / BN_LIMB_BITS)

typedef struct {
    uint32_t limbs[BN_MAX_LIMBS];
    int len; // Number of limbs used
} Bignum;

void bignum_zero(Bignum* bn) {
    memset(bn->limbs, 0, sizeof(bn->limbs));
    bn->len = 0;
}

void bignum_trim(Bignum* bn) {
    while (bn->len > 0 && bn->limbs[bn->len - 1] == 0) {
        bn->len--;
    }
}

void bignum_copy(Bignum* dst, const Bignum* src) {
    dst->len = src->len;
    memcpy(dst->limbs, src->limbs, src->len * sizeof(uint32_t));
}

int bignum_is_zero(const Bignum* bn) {
    return bn->len == 0 || (bn->len == 1 && bn->limbs[0] == 0);
}

int bignum_compare(const Bignum* a, const Bignum* b) {
    if (a->len > b->len) return 1;
    if (a->len < b->len) return -1;
    for (int i = a->len - 1; i >= 0; i--) {
        if (a->limbs[i] > b->limbs[i]) return 1;
        if (a->limbs[i] < b->limbs[i]) return -1;
    }
    return 0;
}

void bignum_from_u64(Bignum* bn, uint64_t val) {
    bignum_zero(bn);
    if (val == 0) {
        bn->len = 0;
        return;
    }
    bn->limbs[0] = (uint32_t)(val & 0xFFFFFFFF);
    bn->limbs[1] = (uint32_t)(val >> 32);
    bn->len = (bn->limbs[1] > 0) ? 2 : 1;
}

void bignum_from_bytes(Bignum* bn, const unsigned char* bytes, size_t num_bytes) {
    bignum_zero(bn);
    bn->len = (num_bytes + 3) / 4;
    if(bn->len > BN_MAX_LIMBS) bn->len = BN_MAX_LIMBS;

    for (size_t i = 0; i < num_bytes; ++i) {
        size_t limb_idx = i / 4;
        size_t limb_shift = (i % 4) * 8;
        bn->limbs[limb_idx] |= ((uint32_t)bytes[i]) << limb_shift;
    }
    bignum_trim(bn);
}

// Basic subtraction: res = a - b. Assumes a >= b.
void bignum_sub(Bignum* res, const Bignum* a, const Bignum* b) {
    long long borrow = 0;
    int max_len = a->len > b->len ? a->len : b->len;
    for (int i = 0; i < max_len; i++) {
        long long a_limb = (i < a->len) ? a->limbs[i] : 0;
        long long b_limb = (i < b->len) ? b->limbs[i] : 0;
        long long diff = a_limb - b_limb - borrow;
        res->limbs[i] = (uint32_t)(diff & 0xFFFFFFFF);
        borrow = (diff < 0) ? 1 : 0;
    }
    res->len = max_len;
    bignum_trim(res);
}

// Basic addition: res = a + b.
void bignum_add(Bignum* res, const Bignum* a, const Bignum* b) {
    unsigned long long carry = 0;
    int max_len = a->len > b->len ? a->len : b->len;
    for (int i = 0; i < max_len; ++i) {
        unsigned long long sum = carry;
        if (i < a->len) sum += a->limbs[i];
        if (i < b->len) sum += b->limbs[i];
        res->limbs[i] = (uint32_t)(sum & 0xFFFFFFFF);
        carry = sum >> 32;
    }
    if (carry > 0 && max_len < BN_MAX_LIMBS) {
        res->limbs[max_len] = carry;
        res->len = max_len + 1;
    } else {
        res->len = max_len;
    }
    bignum_trim(res);
}

// Basic multiplication: res = a * b.
void bignum_mul(Bignum* res, const Bignum* a, const Bignum* b) {
    Bignum temp;
    bignum_zero(&temp);
    if (bignum_is_zero(a) || bignum_is_zero(b)) {
        bignum_zero(res);
        return;
    }

    for (int i = 0; i < b->len; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < a->len; j++) {
            if (i + j >= BN_MAX_LIMBS) continue;
            uint64_t prod = (uint64_t)b->limbs[i] * a->limbs[j] + temp.limbs[i + j] + carry;
            temp.limbs[i + j] = (uint32_t)(prod & 0xFFFFFFFF);
            carry = prod >> 32;
        }
        if (i + a->len < BN_MAX_LIMBS) {
            temp.limbs[i + a->len] += carry;
        }
    }
    temp.len = a->len + b->len;
    bignum_trim(&temp);
    bignum_copy(res, &temp);
}

// Shift left by one bit.
void bignum_lshift1(Bignum* bn) {
    uint32_t carry = 0;
    for (int i = 0; i < bn->len; i++) {
        uint32_t next_carry = (bn->limbs[i] >> 31);
        bn->limbs[i] = (bn->limbs[i] << 1) | carry;
        carry = next_carry;
    }
    if (carry && bn->len < BN_MAX_LIMBS) {
        bn->limbs[bn->len++] = carry;
    }
}

// Shift right by one bit.
void bignum_rshift1(Bignum* bn) {
    uint32_t carry = 0;
    for (int i = bn->len - 1; i >= 0; i--) {
        uint32_t next_carry = (bn->limbs[i] & 1);
        bn->limbs[i] = (bn->limbs[i] >> 1) | (carry << 31);
        carry = next_carry;
    }
    bignum_trim(bn);
}

// Modulo: res = a % mod. (Using slow restoring division)
void bignum_mod(Bignum* res, const Bignum* a, const Bignum* mod) {
    Bignum remainder, temp_mod;
    bignum_copy(&remainder, a);
    bignum_trim(&remainder);

    while (bignum_compare(&remainder, mod) >= 0) {
        bignum_copy(&temp_mod, mod);
        int shifts = 0;
        while(bignum_compare(&remainder, &temp_mod) >= 0) {
             bignum_lshift1(&temp_mod);
             shifts++;
        }
        if (shifts > 0) {
            bignum_rshift1(&temp_mod);
        }
        bignum_sub(&remainder, &remainder, &temp_mod);
    }
    bignum_copy(res, &remainder);
}

// Modular exponentiation: res = base^exp % mod.
void bignum_mod_exp(Bignum* res, const Bignum* base, const Bignum* exp, const Bignum* mod) {
    Bignum current_base, current_exp, temp_res, temp_mul;

    bignum_mod(&current_base, base, mod);
    bignum_copy(&current_exp, exp);
    bignum_from_u64(res, 1);

    while (!bignum_is_zero(&current_exp)) {
        if (current_exp.limbs[0] & 1) {
            bignum_mul(&temp_mul, res, &current_base);
            bignum_mod(res, &temp_mul, mod);
        }
        bignum_mul(&temp_mul, &current_base, &current_base);
        bignum_mod(&current_base, &temp_mul, mod);

        bignum_rshift1(&current_exp);
    }
}

// --- Global Benchmark Data ---

typedef struct {
    // RSA 'Key' components (public only)
    Bignum* n; // Modulus
    Bignum* e; // Public Exponent

    // Data to encrypt
    unsigned char* message;
    size_t data_size;
    size_t block_size;
    size_t num_blocks;

    // Output
    uint64_t result_checksum;
    
    // Parameters
    int key_bit_length;
} BenchmarkData;

static BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <key_bit_length> <data_size_bytes> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.key_bit_length = atoi(argv[1]);
    g_data.data_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    if (g_data.key_bit_length <= 0 || g_data.key_bit_length > BN_MAX_BITS) {
        fprintf(stderr, "Invalid key bit length. Must be > 0 and <= %d\n", BN_MAX_BITS);
        exit(1);
    }

    // Allocate memory for key components
    g_data.n = (Bignum*)malloc(sizeof(Bignum));
    g_data.e = (Bignum*)malloc(sizeof(Bignum));

    // Generate random modulus 'n'
    int n_limbs = (g_data.key_bit_length + BN_LIMB_BITS - 1) / BN_LIMB_BITS;
    bignum_zero(g_data.n);
    for (int i = 0; i < n_limbs; i++) {
        g_data.n->limbs[i] = mt_rand();
    }
    g_data.n->len = n_limbs;
    // Ensure top bit is set for correct length and bottom bit is set for oddness
    g_data.n->limbs[n_limbs - 1] |= (1UL << ((g_data.key_bit_length - 1) % BN_LIMB_BITS));
    g_data.n->limbs[0] |= 1;
    bignum_trim(g_data.n);

    // Set public exponent 'e' to the standard value 65537
    bignum_from_u64(g_data.e, 65537);

    // Generate random message data
    g_data.message = (unsigned char*)malloc(g_data.data_size);
    for (size_t i = 0; i < g_data.data_size; ++i) {
        g_data.message[i] = mt_rand() & 0xFF;
    }

    // Calculate block size for encryption (must be smaller than n)
    g_data.block_size = (g_data.key_bit_length / 8) - 1;
    if (g_data.block_size == 0) g_data.block_size = 1;
    g_data.num_blocks = g_data.data_size / g_data.block_size;

    // Initialize checksum
    g_data.result_checksum = 0;
}

void run_computation() {
    Bignum m, c; // Message and ciphertext bignums for one block

    for (size_t i = 0; i < g_data.num_blocks; ++i) {
        const unsigned char* block_start = g_data.message + (i * g_data.block_size);
        
        // Convert the raw byte block to a bignum
        bignum_from_bytes(&m, block_start, g_data.block_size);

        // The core cryptographic operation: c = m^e mod n
        bignum_mod_exp(&c, &m, g_data.e, g_data.n);
        
        // Accumulate a result to prevent dead code elimination
        g_data.result_checksum += c.limbs[0];
    }
}

void cleanup() {
    free(g_data.n);
    free(g_data.e);
    free(g_data.message);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print final result to stdout
    printf("%llu\n", (unsigned long long)g_data.result_checksum);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
