#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BIG NUMBER ARITHMETIC LIBRARY ---
#define BIGNUM_LIMB_BITS 64
typedef unsigned __int128 uint128_t;

typedef struct {
    uint64_t *limbs;  // Array of limbs (64-bit words)
    size_t capacity;  // Allocated number of limbs
    size_t count;     // Used number of limbs
} bignum_t;

// Forward declarations for bignum functions
static bignum_t* bignum_alloc(size_t capacity);
static void bignum_free(bignum_t *bn);
static void bignum_trim(bignum_t *bn);
static void bignum_copy(bignum_t *dst, const bignum_t *src);
static void bignum_set_ui(bignum_t *bn, uint64_t val);
static int bignum_is_zero(const bignum_t *bn);
static int bignum_cmp(const bignum_t *a, const bignum_t *b);
static int bignum_get_bit(const bignum_t *bn, size_t bit_idx);
static void bignum_mul(bignum_t *res, const bignum_t *a, const bignum_t *b);
static void bignum_mod(bignum_t *rem, const bignum_t *num, const bignum_t *den);
static void bignum_mod_pow(bignum_t *res, const bignum_t *base, const bignum_t *exp, const bignum_t *mod);
static void bignum_gen_rand(bignum_t *bn, size_t bit_length);

// --- GLOBAL DATA --- 
struct {
    int prime_bit_length;
    bignum_t *p, *g, *a_priv, *b_priv, *a_pub, *b_pub, *s_alice, *s_bob;
    uint64_t final_result;
} G_data;

// --- BENCHMARK FUNCTIONS --- 
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <prime_bit_length> <seed>\n", argv[0]);
        exit(1);
    }
    G_data.prime_bit_length = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);
    mt_seed(seed);

    if (G_data.prime_bit_length <= 0) {
        fprintf(stderr, "FATAL: prime_bit_length must be positive.\n");
        exit(1);
    }

    size_t limb_count = (G_data.prime_bit_length + BIGNUM_LIMB_BITS - 1) / BIGNUM_LIMB_BITS;

    G_data.p = bignum_alloc(limb_count);
    G_data.g = bignum_alloc(limb_count);
    G_data.a_priv = bignum_alloc(limb_count);
    G_data.b_priv = bignum_alloc(limb_count);
    G_data.a_pub = bignum_alloc(limb_count);
    G_data.b_pub = bignum_alloc(limb_count);
    G_data.s_alice = bignum_alloc(limb_count);
    G_data.s_bob = bignum_alloc(limb_count);

    // Generate a large prime p
    bignum_gen_rand(G_data.p, G_data.prime_bit_length);
    // Ensure it's odd and has the top bit set for correct bit length
    G_data.p->limbs[0] |= 1;
    G_data.p->limbs[G_data.p->count-1] |= (1ULL << ((G_data.prime_bit_length - 1) % BIGNUM_LIMB_BITS));

    // Use a common generator g=5
    bignum_set_ui(G_data.g, 5);

    // Generate private keys for Alice and Bob (half the bit length of p)
    bignum_gen_rand(G_data.a_priv, G_data.prime_bit_length / 2);
    bignum_gen_rand(G_data.b_priv, G_data.prime_bit_length / 2);

    G_data.final_result = 0;
}

void run_computation() {
    bignum_t *tmp_mul = bignum_alloc(G_data.p->capacity * 2);
    bignum_t *tmp_mod = bignum_alloc(G_data.p->capacity * 2);

    // A = g^a mod p
    bignum_mod_pow(G_data.a_pub, G_data.g, G_data.a_priv, G_data.p);

    // B = g^b mod p
    bignum_mod_pow(G_data.b_pub, G_data.g, G_data.b_priv, G_data.p);

    // s_alice = B^a mod p
    bignum_mod_pow(G_data.s_alice, G_data.b_pub, G_data.a_priv, G_data.p);

    // s_bob = A^b mod p
    bignum_mod_pow(G_data.s_bob, G_data.a_pub, G_data.b_priv, G_data.p);

    // Check that s_alice == s_bob (optional, but good practice)
    // if(bignum_cmp(G_data.s_alice, G_data.s_bob) != 0) {
    //     fprintf(stderr, "Error: Shared secrets do not match!\n");
    // }

    // Use a limb of the result to prevent dead code elimination
    G_data.final_result = G_data.s_alice->limbs[0];

    bignum_free(tmp_mul);
    bignum_free(tmp_mod);
}

void cleanup() {
    bignum_free(G_data.p);
    bignum_free(G_data.g);
    bignum_free(G_data.a_priv);
    bignum_free(G_data.b_priv);
    bignum_free(G_data.a_pub);
    bignum_free(G_data.b_pub);
    bignum_free(G_data.s_alice);
    bignum_free(G_data.s_bob);
}

// --- MAIN --- 
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%llu\n", (unsigned long long)G_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}

// --- BIGNUM IMPLEMENTATION --- 

static bignum_t* bignum_alloc(size_t capacity) {
    bignum_t *bn = (bignum_t*)malloc(sizeof(bignum_t));
    if (!bn) exit(1);
    bn->limbs = (uint64_t*)calloc(capacity, sizeof(uint64_t));
    if (!bn->limbs) { free(bn); exit(1); }
    bn->capacity = capacity;
    bn->count = 0;
    return bn;
}

static void bignum_free(bignum_t *bn) {
    if (bn) {
        free(bn->limbs);
        free(bn);
    }
}

static void bignum_trim(bignum_t *bn) {
    size_t i = bn->capacity;
    while (i > 0 && bn->limbs[i - 1] == 0) {
        i--;
    }
    bn->count = i == 0 ? 1 : i;
}

static void bignum_copy(bignum_t *dst, const bignum_t *src) {
    if (dst->capacity < src->count) {
        fprintf(stderr, "FATAL: bignum_copy destination too small.\n");
        exit(1);
    }
    dst->count = src->count;
    memcpy(dst->limbs, src->limbs, src->count * sizeof(uint64_t));
    if(dst->capacity > src->count) {
        memset(dst->limbs + src->count, 0, (dst->capacity - src->count) * sizeof(uint64_t));
    }
}

static void bignum_set_ui(bignum_t *bn, uint64_t val) {
    memset(bn->limbs, 0, bn->capacity * sizeof(uint64_t));
    bn->limbs[0] = val;
    bn->count = (val == 0) ? 1 : 1;
    bignum_trim(bn);
}

static int bignum_is_zero(const bignum_t *bn) {
    return bn->count == 1 && bn->limbs[0] == 0;
}

static int bignum_cmp(const bignum_t *a, const bignum_t *b) {
    if (a->count > b->count) return 1;
    if (a->count < b->count) return -1;
    for (size_t i = a->count; i > 0; --i) {
        if (a->limbs[i - 1] > b->limbs[i - 1]) return 1;
        if (a->limbs[i - 1] < b->limbs[i - 1]) return -1;
    }
    return 0;
}

static int bignum_get_bit(const bignum_t *bn, size_t bit_idx) {
    size_t limb_idx = bit_idx / BIGNUM_LIMB_BITS;
    size_t bit_in_limb = bit_idx % BIGNUM_LIMB_BITS;
    if (limb_idx >= bn->count) return 0;
    return (bn->limbs[limb_idx] >> bit_in_limb) & 1;
}

static uint64_t bignum_sub(bignum_t *res, const bignum_t *a, const bignum_t *b) {
    uint64_t borrow = 0;
    size_t max_count = a->count > b->count ? a->count : b->count;
    for (size_t i = 0; i < max_count; i++) {
        uint64_t limb_a = (i < a->count) ? a->limbs[i] : 0;
        uint64_t limb_b = (i < b->count) ? b->limbs[i] : 0;
        uint64_t diff = limb_a - limb_b - borrow;
        borrow = (limb_b > limb_a) || (limb_b == limb_a && borrow) || (diff > limb_a) ? 1 : 0;
        if (i < res->capacity) res->limbs[i] = diff;
    }
    res->count = max_count;
    bignum_trim(res);
    return borrow;
}

static void bignum_shl(bignum_t *bn, size_t bits) {
    size_t limb_shift = bits / BIGNUM_LIMB_BITS;
    size_t bit_shift = bits % BIGNUM_LIMB_BITS;

    if (limb_shift > 0) {
        for (size_t i = bn->count; i > 0; --i) {
            bn->limbs[i - 1 + limb_shift] = bn->limbs[i - 1];
        }
        memset(bn->limbs, 0, limb_shift * sizeof(uint64_t));
        bn->count += limb_shift;
    }

    if (bit_shift > 0) {
        uint64_t carry = 0;
        for (size_t i = 0; i < bn->count; i++) {
            uint64_t new_carry = bn->limbs[i] >> (BIGNUM_LIMB_BITS - bit_shift);
            bn->limbs[i] = (bn->limbs[i] << bit_shift) | carry;
            carry = new_carry;
        }
        if (carry > 0) {
            bn->limbs[bn->count++] = carry;
        }
    }
    bignum_trim(bn);
}

static void bignum_mul(bignum_t *res, const bignum_t *a, const bignum_t *b) {
    if (res->capacity < a->count + b->count) {
        // This simplified implementation doesn't support reallocation.
        // The caller must provide a large enough result bignum.
        fprintf(stderr, "FATAL: bignum_mul result capacity too small.\n");
        exit(1);
    }
    memset(res->limbs, 0, res->capacity * sizeof(uint64_t));
    
    for (size_t i = 0; i < b->count; ++i) {
        uint128_t carry = 0;
        for (size_t j = 0; j < a->count; ++j) {
            uint128_t product = (uint128_t)b->limbs[i] * a->limbs[j] + res->limbs[i + j] + carry;
            res->limbs[i + j] = (uint64_t)product;
            carry = product >> BIGNUM_LIMB_BITS;
        }
        if (i + a->count < res->capacity) {
            res->limbs[i + a->count] = carry;
        }
    }
    res->count = a->count + b->count;
    bignum_trim(res);
}

// Basic restoring division for modulus
static void bignum_mod(bignum_t *rem, const bignum_t *num, const bignum_t *den) {
    if (bignum_cmp(num, den) < 0) {
        bignum_copy(rem, num);
        return;
    }

    bignum_t *current_rem = bignum_alloc(num->capacity);
    bignum_copy(current_rem, num);

    size_t den_bits = 0;
    if (!bignum_is_zero(den)) {
        den_bits = (den->count - 1) * BIGNUM_LIMB_BITS;
        uint64_t ms_limb = den->limbs[den->count - 1];
        while (ms_limb > 0) {
            den_bits++;
            ms_limb >>= 1;
        }
    }

    bignum_t *d = bignum_alloc(current_rem->capacity + 1);
    bignum_copy(d, den);

    size_t rem_bits = 0;
    if (!bignum_is_zero(current_rem)) {
        rem_bits = (current_rem->count - 1) * BIGNUM_LIMB_BITS;
        uint64_t ms_limb = current_rem->limbs[current_rem->count - 1];
        while (ms_limb > 0) {
            rem_bits++;
            ms_limb >>= 1;
        }
    }

    if (rem_bits < den_bits) {
         bignum_copy(rem, current_rem);
         bignum_free(d);
         bignum_free(current_rem);
         return;
    }

    int shift = rem_bits - den_bits;
    bignum_shl(d, shift);

    while (shift >= 0) {
        if (bignum_cmp(current_rem, d) >= 0) {
            bignum_sub(current_rem, current_rem, d);
        }
        // Right shift d by 1
        uint64_t carry = 0;
        for (size_t i = d->count; i > 0; --i) {
            uint64_t new_carry = (d->limbs[i - 1] & 1) << (BIGNUM_LIMB_BITS - 1);
            d->limbs[i - 1] = (d->limbs[i - 1] >> 1) | carry;
            carry = new_carry;
        }
        bignum_trim(d);
        shift--;
    }
    
    bignum_copy(rem, current_rem);
    bignum_free(d);
    bignum_free(current_rem);
}

// Computes res = (base^exp) % mod
static void bignum_mod_pow(bignum_t *res, const bignum_t *base, const bignum_t *exp, const bignum_t *mod) {
    bignum_t *b = bignum_alloc(mod->capacity);
    bignum_mod(b, base, mod);

    bignum_set_ui(res, 1);

    bignum_t *temp_mul = bignum_alloc(mod->capacity * 2);

    size_t num_bits = 0;
    if (!bignum_is_zero(exp)) {
        num_bits = (exp->count - 1) * BIGNUM_LIMB_BITS;
        uint64_t ms_limb = exp->limbs[exp->count - 1];
        while (ms_limb > 0) {
            num_bits++;
            ms_limb >>= 1;
        }
    }
    
    for (size_t i = num_bits; i > 0; --i) {
        // res = (res * res) % mod
        bignum_mul(temp_mul, res, res);
        bignum_mod(res, temp_mul, mod);

        if (bignum_get_bit(exp, i - 1)) {
            // res = (res * b) % mod
            bignum_mul(temp_mul, res, b);
            bignum_mod(res, temp_mul, mod);
        }
    }

    bignum_free(b);
    bignum_free(temp_mul);
}

static void bignum_gen_rand(bignum_t *bn, size_t bit_length) {
    size_t limb_count = (bit_length + BIGNUM_LIMB_BITS - 1) / BIGNUM_LIMB_BITS;
    if (bn->capacity < limb_count) {
        fprintf(stderr, "FATAL: bignum_gen_rand capacity too small.\n");
        exit(1);
    }
    memset(bn->limbs, 0, bn->capacity * sizeof(uint64_t));
    
    for(size_t i = 0; i < limb_count; ++i) {
        bn->limbs[i] = ((uint64_t)mt_rand() << 32) | mt_rand();
    }
 
    size_t remaining_bits = bit_length % BIGNUM_LIMB_BITS;
    if (remaining_bits > 0) {
        bn->limbs[limb_count - 1] &= (1ULL << remaining_bits) - 1;
    }
    bn->count = limb_count;
    bignum_trim(bn);
}
