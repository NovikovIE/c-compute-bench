#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// Mersenne Twister (MT19937) Generator - DO NOT MODIFY
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

// --- BENCHMARK SPECIFIC --- //

#define LIMB_TYPE uint32_t
#define LIMB_BITS 32
#define MAX_LIMBS 256 // Corresponds to 8192 bits, max for this benchmark

// --- Data Structures ---
typedef struct {
    LIMB_TYPE limbs[MAX_LIMBS];
    int size; // Number of active limbs
} BigInt;

typedef struct {
    BigInt x;
    BigInt y;
    int is_infinity;
} ECPoint;

typedef struct {
    BigInt p; // Prime modulus
    BigInt a; // Curve parameter a
    ECPoint G; // Base point
    BigInt n; // Order of G
} ECCurve;

// --- Global State ---
typedef struct {
    int curve_bit_length;
    int num_limbs;
    ECCurve curve;
    BigInt private_key;
    ECPoint public_key;
} BenchmarkData;

static BenchmarkData* g_data;
static volatile uint32_t g_result_accumulator;

// --- BigInt Utility Functions ---
void bigint_init(BigInt* bi, int num_limbs) {
    memset(bi->limbs, 0, MAX_LIMBS * sizeof(LIMB_TYPE));
    bi->size = num_limbs;
}

void bigint_copy(BigInt* dest, const BigInt* src) {
    memcpy(dest->limbs, src->limbs, src->size * sizeof(LIMB_TYPE));
    dest->size = src->size;
}

int bigint_is_zero(const BigInt* bi) {
    for (int i = 0; i < bi->size; i++) {
        if (bi->limbs[i] != 0) return 0;
    }
    return 1;
}

void bigint_set_word(BigInt* bi, uint64_t word) {
    bigint_init(bi, bi->size);
    bi->limbs[0] = (LIMB_TYPE)(word & 0xFFFFFFFF);
    if (LIMB_BITS < 64 && bi->size > 1) {
        bi->limbs[1] = (LIMB_TYPE)(word >> 32);
    }
}

int bigint_get_bit(const BigInt* bi, int k) {
    if (k >= bi->size * LIMB_BITS) return 0;
    int limb_idx = k / LIMB_BITS;
    int bit_idx = k % LIMB_BITS;
    return (bi->limbs[limb_idx] >> bit_idx) & 1;
}

void bigint_rand(BigInt* bi) {
    for (int i = 0; i < bi->size; i++) {
        bi->limbs[i] = mt_rand();
    }
    // Ensure it's of the correct bit length
    bi->limbs[bi->size - 1] |= (1UL << (LIMB_BITS - 1));
}

// --- BigInt Modular Arithmetic (Simplified for benchmark) ---
// These are NOT cryptographically secure, just compute-intensive proxies.

void bigint_mod_add(BigInt* res, const BigInt* a, const BigInt* b, const BigInt* m) {
    uint64_t carry = 0;
    for (int i = 0; i < m->size; i++) {
        uint64_t sum = (uint64_t)a->limbs[i] + b->limbs[i] + carry;
        res->limbs[i] = (LIMB_TYPE)sum;
        carry = sum >> LIMB_BITS;
    }
    // Simple modular reduction: if result > m, subtract m.
    if (carry) {
        uint64_t borrow = 0;
        for (int i = 0; i < m->size; i++) {
            uint64_t diff = (uint64_t)res->limbs[i] - m->limbs[i] - borrow;
            res->limbs[i] = (LIMB_TYPE)diff;
            borrow = (diff >> 63) & 1;
        }
    }
}

void bigint_mod_sub(BigInt* res, const BigInt* a, const BigInt* b, const BigInt* m) {
    int64_t borrow = 0;
    for (int i = 0; i < m->size; i++) {
        int64_t diff = (int64_t)a->limbs[i] - b->limbs[i] - borrow;
        res->limbs[i] = (LIMB_TYPE)diff;
        borrow = (diff < 0);
    }
    if (borrow) {
        uint64_t carry = 0;
        for (int i = 0; i < m->size; i++) {
            uint64_t sum = (uint64_t)res->limbs[i] + m->limbs[i] + carry;
            res->limbs[i] = (LIMB_TYPE)sum;
            carry = sum >> LIMB_BITS;
        }
    }
}

// In-place right shift by 1 bit
void bigint_rshift1(BigInt *a) {
    for (int i = 0; i < a->size - 1; i++) {
        a->limbs[i] = (a->limbs[i] >> 1) | (a->limbs[i+1] << (LIMB_BITS-1));
    }
    a->limbs[a->size - 1] >>= 1;
}

void bigint_mul_mod(BigInt* res, const BigInt* a, const BigInt* b, const BigInt* m) {
    BigInt temp_res;
    bigint_init(&temp_res, m->size);
    
    BigInt temp_a;
    bigint_copy(&temp_a, a);
    
    for (int i = 0; i < m->size * LIMB_BITS; i++) {
        if (bigint_get_bit(b, i)) {
            bigint_mod_add(&temp_res, &temp_res, &temp_a, m);
        }
        uint64_t carry = (temp_a.limbs[m->size-1] >> (LIMB_BITS-1)) & 1;
        for (int j = m->size - 1; j > 0; j--) {
            temp_a.limbs[j] = (temp_a.limbs[j] << 1) | (temp_a.limbs[j-1] >> (LIMB_BITS-1));
        }
        temp_a.limbs[0] <<= 1;
        if (carry) {
            bigint_mod_sub(&temp_a, &temp_a, m, m);
        }
    }
    bigint_copy(res, &temp_res);
}

void bigint_pow_mod(BigInt* res, const BigInt* base, const BigInt* exp, const BigInt* m) {
    BigInt temp_res, temp_base, temp_exp;
    bigint_init(&temp_res, m->size);
    bigint_init(&temp_base, m->size);
    bigint_init(&temp_exp, m->size);

    bigint_set_word(&temp_res, 1);
    bigint_copy(&temp_base, base);
    bigint_copy(&temp_exp, exp);

    while (!bigint_is_zero(&temp_exp)) {
        if (temp_exp.limbs[0] & 1) {
            bigint_mul_mod(&temp_res, &temp_res, &temp_base, m);
        }
        bigint_mul_mod(&temp_base, &temp_base, &temp_base, m);
        bigint_rshift1(&temp_exp);
    }
    bigint_copy(res, &temp_res);
}

// Using Fermat's Little Theorem: a^(p-2) mod p
void bigint_inv_mod(BigInt* res, const BigInt* a, const BigInt* p) {
    BigInt exp;
    bigint_copy(&exp, p);
    bigint_mod_sub(&exp, &exp, &(BigInt){.limbs={2,0}, .size=p->size}, p);
    bigint_pow_mod(res, a, &exp, p);
}

// --- Elliptic Curve Arithmetic ---
void ec_point_double(ECPoint* res, const ECPoint* p, const ECCurve* curve) {
    if (p->is_infinity) {
        res->is_infinity = 1;
        return;
    }
    BigInt lambda, temp1, temp2, temp3;
    int num_limbs = curve->p.size;
    bigint_init(&lambda, num_limbs);
    bigint_init(&temp1, num_limbs);
    bigint_init(&temp2, num_limbs);
    bigint_init(&temp3, num_limbs);

    // lambda = (3*x^2 + a) / (2*y)
    // Denominator: 2*y
    bigint_mod_add(&temp1, &p->y, &p->y, &curve->p);
    bigint_inv_mod(&temp1, &temp1, &curve->p);

    // Numerator: 3*x^2 + a
    bigint_mul_mod(&temp2, &p->x, &p->x, &curve->p);
    bigint_mod_add(&temp3, &temp2, &temp2, &curve->p);
    bigint_mod_add(&temp2, &temp3, &temp2, &curve->p);
    bigint_mod_add(&temp2, &temp2, &curve->a, &curve->p);

    bigint_mul_mod(&lambda, &temp2, &temp1, &curve->p);

    // rx = lambda^2 - 2*x
    bigint_mul_mod(&temp1, &lambda, &lambda, &curve->p);
    bigint_mod_add(&temp2, &p->x, &p->x, &curve->p);
    bigint_mod_sub(&res->x, &temp1, &temp2, &curve->p);

    // ry = lambda(px - rx) - py
    bigint_mod_sub(&temp1, &p->x, &res->x, &curve->p);
    bigint_mul_mod(&temp2, &lambda, &temp1, &curve->p);
    bigint_mod_sub(&res->y, &temp2, &p->y, &curve->p);
    res->is_infinity = 0;
}

void ec_point_add(ECPoint* res, const ECPoint* p1, const ECPoint* p2, const ECCurve* curve) {
    if (p1->is_infinity) { *res = *p2; return; }
    if (p2->is_infinity) { *res = *p1; return; }

    if (memcmp(p1->x.limbs, p2->x.limbs, p1->x.size * sizeof(LIMB_TYPE)) == 0 &&
        memcmp(p1->y.limbs, p2->y.limbs, p1->y.size * sizeof(LIMB_TYPE)) == 0) {
        ec_point_double(res, p1, curve);
        return;
    }

    BigInt lambda, temp1, temp2;
    int num_limbs = curve->p.size;
    bigint_init(&lambda, num_limbs);
    bigint_init(&temp1, num_limbs);
    bigint_init(&temp2, num_limbs);

    // lambda = (y2 - y1) / (x2 - x1)
    bigint_mod_sub(&temp1, &p2->x, &p1->x, &curve->p);
    bigint_inv_mod(&temp1, &temp1, &curve->p);

    bigint_mod_sub(&temp2, &p2->y, &p1->y, &curve->p);
    bigint_mul_mod(&lambda, &temp2, &temp1, &curve->p);

    // rx = lambda^2 - x1 - x2
    bigint_mul_mod(&temp1, &lambda, &lambda, &curve->p);
    bigint_mod_sub(&temp1, &temp1, &p1->x, &curve->p);
    bigint_mod_sub(&res->x, &temp1, &p2->x, &curve->p);

    // ry = lambda(x1 - rx) - y1
    bigint_mod_sub(&temp1, &p1->x, &res->x, &curve->p);
    bigint_mul_mod(&temp2, &lambda, &temp1, &curve->p);
    bigint_mod_sub(&res->y, &temp2, &p1->y, &curve->p);
    res->is_infinity = 0;
}

// --- BENCHMARK CORE ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <curve_bit_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("malloc");
        exit(1);
    }

    g_data->curve_bit_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);

    if (g_data->curve_bit_length <= 0 || g_data->curve_bit_length > MAX_LIMBS * LIMB_BITS) {
        fprintf(stderr, "Invalid curve_bit_length. Max: %d\n", MAX_LIMBS * LIMB_BITS);
        exit(1);
    }

    g_data->num_limbs = (g_data->curve_bit_length + LIMB_BITS - 1) / LIMB_BITS;

    mt_seed(seed);

    bigint_init(&g_data->curve.p, g_data->num_limbs);
    bigint_init(&g_data->curve.a, g_data->num_limbs);
    bigint_init(&g_data->curve.G.x, g_data->num_limbs);
    bigint_init(&g_data->curve.G.y, g_data->num_limbs);
    bigint_init(&g_data->private_key, g_data->num_limbs);
    bigint_init(&g_data->public_key.x, g_data->num_limbs);
    bigint_init(&g_data->public_key.y, g_data->num_limbs);
    g_data->curve.G.is_infinity = 0;
    g_data->public_key.is_infinity = 1;

    // Generate random curve parameters and keys
    bigint_rand(&g_data->curve.p);
    g_data->curve.p.limbs[0] |= 1; // Make it odd
    bigint_rand(&g_data->curve.a);
    bigint_rand(&g_data->curve.G.x);
    bigint_rand(&g_data->curve.G.y);
    bigint_rand(&g_data->private_key);
}

void run_computation() {
    ECPoint current_point = g_data->curve.G;
    ECPoint result_point;
    result_point.is_infinity = 1;
    bigint_init(&result_point.x, g_data->num_limbs);
    bigint_init(&result_point.y, g_data->num_limbs);

    // Double-and-add algorithm for scalar multiplication
    for (int i = 0; i < g_data->curve_bit_length; ++i) {
        if (bigint_get_bit(&g_data->private_key, i)) {
            ec_point_add(&result_point, &result_point, &current_point, &g_data->curve);
        }
        ec_point_double(&current_point, &current_point, &g_data->curve);
    }

    g_data->public_key = result_point;

    // Calculate a checksum to prevent dead code elimination
    g_result_accumulator = 0;
    for (int i = 0; i < g_data->num_limbs; ++i) {
        g_result_accumulator ^= g_data->public_key.x.limbs[i];
        g_result_accumulator ^= g_data->public_key.y.limbs[i];
    }
}

void cleanup() {
    free(g_data);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%u\n", g_result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}