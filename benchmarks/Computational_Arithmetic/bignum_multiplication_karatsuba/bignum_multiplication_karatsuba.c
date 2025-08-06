#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

// Mersenne Twister (MT19937) PRNG (verbatim as requested)
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

// --- Benchmark-specific code ---

// Karatsuba multiplication cutoff. Below this size, classic multiplication is used.
#define KARATSUBA_CUTOFF 32

// Struct to hold all benchmark data
typedef struct {
    int num_digits_a;
    int num_digits_b;
    int *num_a;
    int *num_b;
    int *product;
    int final_result; // Checksum of the product
} BenchmarkData;

static BenchmarkData g_data;

// Forward declarations
static int bignum_trim(const int *num, int n);
static void bignum_mult_classic(const int *a, int n_a, const int *b, int n_b, int *res, int *n_res);
static void karatsuba_recursive(const int *a, const int *b, int *res, int n);

// Trims leading zeros from a bignum and returns its new effective length
static int bignum_trim(const int *num, int n) {
    int i = n - 1;
    while (i > 0 && num[i] == 0) {
        i--;
    }
    return i + 1;
}

// Adds two big numbers: res = a + b
static void bignum_add(const int *a, int n_a, const int *b, int n_b, int *res, int *n_res) {
    int max_n = n_a > n_b ? n_a : n_b;
    int carry = 0;
    for (int i = 0; i < max_n; i++) {
        int digit_a = (i < n_a) ? a[i] : 0;
        int digit_b = (i < n_b) ? b[i] : 0;
        int sum = digit_a + digit_b + carry;
        res[i] = sum % 10;
        carry = sum / 10;
    }
    if (carry > 0) {
        res[max_n] = carry;
        *n_res = max_n + 1;
    } else {
        *n_res = max_n;
    }
}

// Subtracts two big numbers: res = a - b (assumes a >= b)
static void bignum_sub(const int *a, int n_a, const int *b, int n_b, int *res, int *n_res) {
    int borrow = 0;
    for (int i = 0; i < n_a; i++) {
        int digit_b = (i < n_b) ? b[i] : 0;
        int diff = a[i] - digit_b - borrow;
        if (diff < 0) {
            diff += 10;
            borrow = 1;
        } else {
            borrow = 0;
        }
        res[i] = diff;
    }
    *n_res = bignum_trim(res, n_a);
}

// Adds 'source' to 'target' in-place at a given shift, handling carries
static void bignum_add_to_shifted(int *target, int target_size, const int *source, int source_len, int shift) {
    int carry = 0;
    int i;
    for (i = 0; i < source_len; i++) {
        if (i + shift >= target_size) break; // Bounds check
        int sum = target[i + shift] + source[i] + carry;
        target[i + shift] = sum % 10;
        carry = sum / 10;
    }
    while (carry > 0 && i + shift < target_size) {
        int sum = target[i + shift] + carry;
        target[i + shift] = sum % 10;
        carry = sum / 10;
        i++;
    }
}

// Classic O(n^2) schoolbook multiplication
static void bignum_mult_classic(const int *a, int n_a, const int *b, int n_b, int *res, int *n_res) {
    int len = n_a + n_b;
    memset(res, 0, sizeof(int) * len);
    for (int i = 0; i < n_a; i++) {
        for (int j = 0; j < n_b; j++) {
            res[i + j] += a[i] * b[j];
        }
    }
    int carry = 0;
    for (int i = 0; i < len; i++) {
        int temp = res[i] + carry;
        res[i] = temp % 10;
        carry = temp / 10;
    }
    *n_res = bignum_trim(res, len);
}

// Karatsuba multiplication algorithm (recursive implementation)
static void karatsuba_recursive(const int *a, const int *b, int *res, int n) {
    memset(res, 0, sizeof(int) * (2 * n + 2));
    if (n <= KARATSUBA_CUTOFF) {
        int dummy_len;
        bignum_mult_classic(a, n, b, n, res, &dummy_len);
        return;
    }

    int m = n / 2;
    int m_high = n - m;

    const int *a_low = a;
    const int *a_high = a + m;
    const int *b_low = b;
    const int *b_high = b + m;

    // z0 = a_low * b_low
    // z2 = a_high * b_high
    int *z0 = (int *)calloc(2 * m + 2, sizeof(int));
    int *z2 = (int *)calloc(2 * m_high + 2, sizeof(int));
    karatsuba_recursive(a_low, b_low, z0, m);
    karatsuba_recursive(a_high, b_high, z2, m_high);

    // p1 = (a_low + a_high) * (b_low + b_high)
    int *sum_a = (int *)calloc(m_high + 1, sizeof(int));
    int *sum_b = (int *)calloc(m_high + 1, sizeof(int));
    int n_sum_a, n_sum_b;
    bignum_add(a_low, m, a_high, m_high, sum_a, &n_sum_a);
    bignum_add(b_low, m, b_high, m_high, sum_b, &n_sum_b);
    
    int n_sum = n_sum_a > n_sum_b ? n_sum_a : n_sum_b;
    int *p1 = (int *)calloc(2 * n_sum + 2, sizeof(int));
    karatsuba_recursive(sum_a, sum_b, p1, n_sum);
    free(sum_a);
    free(sum_b);

    // z1 = p1 - z2 - z0
    int n_z0 = bignum_trim(z0, 2 * m + 2);
    int n_z2 = bignum_trim(z2, 2 * m_high + 2);
    int n_p1 = bignum_trim(p1, 2 * n_sum + 2);

    int sub_len1, sub_len2;
    bignum_sub(p1, n_p1, z2, n_z2, p1, &sub_len1);
    bignum_sub(p1, sub_len1, z0, n_z0, p1, &sub_len2);
    int n_z1 = sub_len2;

    // Result = z2*B^2m + z1*B^m + z0
    bignum_add_to_shifted(res, 2 * n + 2, z0, n_z0, 0);
    bignum_add_to_shifted(res, 2 * n + 2, p1, n_z1, m);
    bignum_add_to_shifted(res, 2 * n + 2, z2, n_z2, 2 * m);

    free(z0);
    free(z2);
    free(p1);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_digits_a> <num_digits_b> <seed>\n", argv[0]); exit(1);
    }
    g_data.num_digits_a = atoi(argv[1]);
    g_data.num_digits_b = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_digits_a <= 0 || g_data.num_digits_b <= 0) {
        fprintf(stderr, "FATAL: Number of digits must be positive.\n"); exit(1);
    }
    mt_seed(seed);

    g_data.num_a = (int*)malloc(g_data.num_digits_a * sizeof(int));
    g_data.num_b = (int*)malloc(g_data.num_digits_b * sizeof(int));
    int max_prod_digits = g_data.num_digits_a + g_data.num_digits_b + 2;
    g_data.product = (int*)malloc(max_prod_digits * sizeof(int));

    if (!g_data.num_a || !g_data.num_b || !g_data.product) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n"); exit(1);
    }

    for (int i = 0; i < g_data.num_digits_a; i++) g_data.num_a[i] = mt_rand() % 10;
    if (g_data.num_a[g_data.num_digits_a - 1] == 0) g_data.num_a[g_data.num_digits_a - 1] = 1;

    for (int i = 0; i < g_data.num_digits_b; i++) g_data.num_b[i] = mt_rand() % 10;
    if (g_data.num_b[g_data.num_digits_b - 1] == 0) g_data.num_b[g_data.num_digits_b - 1] = 1;

    g_data.final_result = 0;
}

void run_computation() {
    int n_a = g_data.num_digits_a;
    int n_b = g_data.num_digits_b;
    int n = n_a > n_b ? n_a : n_b;

    // Pad smaller number to match length of larger number
    int *a_pad = (int*)calloc(n, sizeof(int));
    int *b_pad = (int*)calloc(n, sizeof(int));
    memcpy(a_pad, g_data.num_a, n_a * sizeof(int));
    memcpy(b_pad, g_data.num_b, n_b * sizeof(int));

    // Result buffer needs to be large enough for the product
    int res_size = 2 * n + 2;
    int *res_pad = (int*)calloc(res_size, sizeof(int));
    
    karatsuba_recursive(a_pad, b_pad, res_pad, n);
    
    int final_len = bignum_trim(res_pad, res_size);
    memcpy(g_data.product, res_pad, final_len * sizeof(int));
    
    // Calculate a checksum of the result digits to prevent dead-code elimination
    for (int i = 0; i < final_len; i++) {
        g_data.final_result = (g_data.final_result + g_data.product[i]) % 1000;
    }
    
    free(a_pad);
    free(b_pad);
    free(res_pad);
}

void cleanup() {
    free(g_data.num_a);
    free(g_data.num_b);
    free(g_data.product);
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
    printf("%d\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
