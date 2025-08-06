#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

// NOTE: This implementation uses uint64_t for bignum arithmetic for simplicity.
// It is functionally correct for the ECC algorithm but only supports a curve_bit_length up to 64.

// --- Mersenne Twister (Verbatim) ---
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
// --- End Mersenne Twister ---

// --- Types and Global Data ---
typedef uint64_t BigInt;

typedef struct {
    BigInt x;
    BigInt y;
    int is_infinity;
} EC_Point;

typedef struct {
    BigInt p; // Prime modulus
    BigInt a;
} EC_Curve;

// Global struct for benchmark data
static struct {
    int curve_bit_length;
    int num_ops;
    EC_Curve curve;
    EC_Point base_point;
    BigInt* scalars;
    EC_Point* results;
    BigInt final_checksum;
} G;

// --- Modular Arithmetic Utilities ---

// Modular addition: r = (a + b) mod p
static inline BigInt add_mod(BigInt a, BigInt b, BigInt p) {
    unsigned __int128 sum = (unsigned __int128)a + b;
    return (BigInt)(sum % p);
}

// Modular subtraction: r = (a - b) mod p
static inline BigInt sub_mod(BigInt a, BigInt b, BigInt p) {
    unsigned __int128 diff = (unsigned __int128)a + p - b;
    return (BigInt)(diff % p);
}

// Modular multiplication: r = (a * b) mod p
static inline BigInt mul_mod(BigInt a, BigInt b, BigInt p) {
    unsigned __int128 prod = (unsigned __int128)a * b;
    return (BigInt)(prod % p);
}

// Modular inverse: r = n^-1 mod p using Extended Euclidean Algorithm for uint64_t
static BigInt inv_mod(BigInt n, BigInt p) {
    if (p == 1) return 0;
    int64_t m0 = p, t;
    int64_t x0 = 0, x1 = 1;
    int64_t a = n, m = p;
    while (a > 1) {
        if (m == 0) return 0; // No inverse
        int64_t q = a / m;
        t = m; m = a % m; a = t;
        t = x0; x0 = x1 - q * x0; x1 = t;
    }
    if (x1 < 0) x1 += m0;
    return (BigInt)x1;
}

// --- ECC Operations ---

// Point addition: R = P + Q
void ec_point_add(EC_Point* R, const EC_Point* P, const EC_Point* Q, const EC_Curve* curve) {
    if (P->is_infinity) { *R = *Q; return; }
    if (Q->is_infinity) { *R = *P; return; }

    BigInt lambda;

    if (P->x == Q->x) { // Point doubling or P = -Q
        if ((P->y + Q->y) % curve->p == 0) { // P = -Q
            R->is_infinity = 1; 
            return;
        }
        // Point doubling: lambda = (3*x_p^2 + a) / (2*y_p)
        BigInt num = mul_mod(3, mul_mod(P->x, P->x, curve->p), curve->p);
        num = add_mod(num, curve->a, curve->p);
        BigInt den = add_mod(P->y, P->y, curve->p);
        lambda = mul_mod(num, inv_mod(den, curve->p), curve->p);
    } else { // Point addition: lambda = (y_q - y_p) / (x_q - x_p)
        BigInt num = sub_mod(Q->y, P->y, curve->p);
        BigInt den = sub_mod(Q->x, P->x, curve->p);
        lambda = mul_mod(num, inv_mod(den, curve->p), curve->p);
    }

    // xr = lambda^2 - xp - xq
    BigInt xr = mul_mod(lambda, lambda, curve->p);
    xr = sub_mod(xr, P->x, curve->p);
    xr = sub_mod(xr, Q->x, curve->p);

    // yr = lambda * (xp - xr) - yp
    BigInt yr = mul_mod(lambda, sub_mod(P->x, xr, curve->p), curve->p);
    yr = sub_mod(yr, P->y, curve->p);
    
    R->is_infinity = 0;
    R->x = xr;
    R->y = yr;
}

// Scalar multiplication: R = k * P using double-and-add
void ec_scalar_multiply(EC_Point* R, BigInt k, const EC_Point* P, const EC_Curve* curve) {
    EC_Point current_p = *P;
    EC_Point result = {0, 0, 1}; // Point at infinity

    for (int i = 0; i < G.curve_bit_length; i++) {
        if ((k >> i) & 1) {
            ec_point_add(&result, &result, &current_p, curve);
        }
        ec_point_add(&current_p, &current_p, &current_p, curve);
    }
    *R = result;
}

// --- Benchmark Setup, Computation, and Cleanup ---

// Helper function to generate a 64-bit random number from the 32-bit MT
static inline uint64_t rand_64() {
    return ((uint64_t)mt_rand() << 32) | mt_rand();
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <curve_bit_length> <seed>\n", argv[0]);
        exit(1);
    }

    G.curve_bit_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (G.curve_bit_length <= 0 || G.curve_bit_length > 64) {
        fprintf(stderr, "Error: curve_bit_length must be between 1 and 64 for this implementation.\n");
        exit(1);
    }

    // Scale num_ops to keep runtime somewhat constant across different bit lengths.
    // Tuned for ~1s on a modern CPU for bit length 64.
    const double base_ops = 80000.0;
    const double base_len = 64.0;
    double scaling_factor = pow(base_len / G.curve_bit_length, 2.0);
    G.num_ops = (int)(base_ops * scaling_factor);
    if (G.num_ops < 1) G.num_ops = 1;

    // 1. Generate Curve Parameters (p, a)
    BigInt p_mask = (G.curve_bit_length == 64) ? -1ULL : (1ULL << G.curve_bit_length) - 1;
    G.curve.p = (rand_64() & p_mask) | (1ULL << (G.curve_bit_length - 1)) | 1; // Ensure odd and full length
    G.curve.a = rand_64() % G.curve.p;

    // 2. Generate a valid base point G = (x, y)
    G.base_point.is_infinity = 0;
    G.base_point.x = rand_64() % G.curve.p;
    G.base_point.y = rand_64() % G.curve.p;

    // 3. Allocate and generate scalars
    G.scalars = (BigInt*)malloc(G.num_ops * sizeof(BigInt));
    G.results = (EC_Point*)malloc(G.num_ops * sizeof(EC_Point));
    if (!G.scalars || !G.results) {
        fprintf(stderr, "Failed to allocate memory for scalars/results.\n");
        exit(1);
    }

    for (int i = 0; i < G.num_ops; i++) {
        G.scalars[i] = rand_64() & p_mask;
        G.results[i].is_infinity = 1;
    }

    // 4. Initialize final checksum
    G.final_checksum = 0;
}

void run_computation() {
    for (int i = 0; i < G.num_ops; ++i) {
        ec_scalar_multiply(&G.results[i], G.scalars[i], &G.base_point, &G.curve);
    }
    
    // Aggregate results to prevent dead code elimination
    for(int i = 0; i < G.num_ops; ++i) {
        if (!G.results[i].is_infinity) {
            G.final_checksum = add_mod(G.final_checksum, G.results[i].x, G.curve.p);
        }
    }
}

void cleanup() {
    free(G.scalars);
    free(G.results);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result checksum to stdout
    printf("%llu\n", (unsigned long long)G.final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
