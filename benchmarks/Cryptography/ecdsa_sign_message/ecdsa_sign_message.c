#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- BENCHMARK DATA AND GLOBALS ---

int curve_bit_length_g;
int num_limbs_g;

// Big number representations
uint32_t* base_g;     // Represents the message hash 'z'
uint32_t* exp_g;      // Represents the private key 'd'
uint32_t* res_g;      // Stores the final result of the modular exponentiation

// Workspace arrays for computation
uint32_t* temp1_g;
uint32_t* temp2_g;

volatile uint32_t final_result_g;

// --- HELPER FUNCTIONS ---

static void bignum_copy(uint32_t* dst, const uint32_t* src, int n) {
    memcpy(dst, src, n * sizeof(uint32_t));
}

static void bignum_set(uint32_t* dst, uint32_t val, int n) {
    memset(dst, 0, n * sizeof(uint32_t));
    dst[0] = val;
}

static int bignum_is_bit_set(const uint32_t* num, int bit_idx) {
    int limb_idx = bit_idx / 32;
    int bit_in_limb = bit_idx % 32;
    if (limb_idx >= num_limbs_g) return 0;
    return (num[limb_idx] >> bit_in_limb) & 1;
}

// This is the computational kernel. It's a CPU-bound workload designed to have
// a data access pattern and instruction mix similar to big number multiplication.
// It is NOT a cryptographically correct modular multiplication.
static void compute_intensive_kernel(uint32_t* res, const uint32_t* a, const uint32_t* b, int n) {
    uint64_t c1 = 0xface;
    uint64_t c2 = 0xbeef;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            c1 += (uint64_t)a[j] * b[i] + res[(i + j) % n];
            c2 += (uint64_t)a[i] * b[j] + res[(i * j * 3) % n];
            res[(i * j) % n] = (uint32_t)(c1 >> 32) ^ (uint32_t)c2;
            res[(i + j) % n] = (uint32_t)c1 + (uint32_t)(c2 >> 32);
        }
    }
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <curve_bit_length> <hash_size_bytes> <seed>\n", argv[0]);
        exit(1);
    }

    curve_bit_length_g = atoi(argv[1]);
    int hash_size_bytes = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (curve_bit_length_g <= 0 || curve_bit_length_g % 32 != 0) {
        fprintf(stderr, "curve_bit_length must be a positive multiple of 32.\n");
        exit(1);
    }

    mt_seed(seed);

    num_limbs_g = curve_bit_length_g / 32;
    int hash_limbs = hash_size_bytes / 4;
    if (hash_limbs <= 0) hash_limbs = 1;
    if (hash_limbs > num_limbs_g) hash_limbs = num_limbs_g;

    base_g = (uint32_t*)malloc(num_limbs_g * sizeof(uint32_t));
    exp_g = (uint32_t*)malloc(num_limbs_g * sizeof(uint32_t));
    res_g = (uint32_t*)malloc(num_limbs_g * sizeof(uint32_t));
    temp1_g = (uint32_t*)malloc(num_limbs_g * sizeof(uint32_t));
    temp2_g = (uint32_t*)malloc(num_limbs_g * sizeof(uint32_t));
    
    if (!base_g || !exp_g || !res_g || !temp1_g || !temp2_g) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize data with random values
    for (int i = 0; i < num_limbs_g; ++i) {
        base_g[i] = (i < hash_limbs) ? mt_rand() : 0;
        exp_g[i] = mt_rand();
    }

    final_result_g = 0;
}

void run_computation() {
    // Simulate modular exponentiation: res = (base ^ exp) mod N
    // The workload is representative of the core operations in an ECDSA signature,
    // specifically modular inversion and modular multiplications.
    uint32_t* power = temp1_g;
    uint32_t* result_acc = res_g;
    uint32_t* temp_b = temp2_g;

    bignum_copy(power, base_g, num_limbs_g);
    bignum_set(result_acc, 1, num_limbs_g);

    int max_bits = curve_bit_length_g;
    for (int i = 0; i < max_bits; ++i) {
        if (bignum_is_bit_set(exp_g, i)) {
            // result_acc = (result_acc * power)
            bignum_copy(temp_b, result_acc, num_limbs_g);
            compute_intensive_kernel(result_acc, temp_b, power, num_limbs_g);
        }
        // power = (power * power)
        bignum_copy(temp_b, power, num_limbs_g);
        compute_intensive_kernel(power, temp_b, temp_b, num_limbs_g);
    }

    // Accumulate a result to prevent dead code elimination
    final_result_g = result_acc[0];
}

void cleanup() {
    free(base_g);
    free(exp_g);
    free(res_g);
    free(temp1_g);
    free(temp2_g);
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
    printf("%u\n", final_result_g);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
