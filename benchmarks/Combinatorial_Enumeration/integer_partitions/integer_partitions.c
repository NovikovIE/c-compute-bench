#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator - Do Not Modify ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
typedef struct {
    int target_integer;
    unsigned long long* p_table; // Table for memoization of partition counts
    unsigned long long final_result;
} benchmark_data;

static benchmark_data g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <target_integer> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.target_integer = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_data.target_integer < 0) {
        fprintf(stderr, "Error: target_integer must be non-negative.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.p_table = (unsigned long long*)malloc((g_data.target_integer + 1) * sizeof(unsigned long long));
    if (!g_data.p_table) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    g_data.final_result = 0;
}

void run_computation() {
    int n = g_data.target_integer;
    unsigned long long* p = g_data.p_table;

    p[0] = 1; // There's one way to partition 0 (the empty sum)

    // Using Euler's pentagonal number theorem for recurrence relation:
    // p(n) = p(n-1) + p(n-2) - p(n-5) - p(n-7) + ...
    // The terms subtracted are p(n-k) where k are generalized pentagonal numbers.
    for (int i = 1; i <= n; i++) {
        p[i] = 0;
        for (int k = 1; ; k++) {
            long long pent1 = (long long)k * (3 * k - 1) / 2;
            long long pent2 = (long long)k * (3 * k + 1) / 2;
            int term1 = i - (int)pent1;
            int term2 = i - (int)pent2;

            if (term1 < 0 && term2 < 0) {
                break;
            }
            
            // The sign is positive for k=1, -1, 3, -3... and negative for k=2, -2, 4, -4...
            // This corresponds to k being odd or even in our loop.
            if ((k - 1) % 2 == 0) { // k is odd, sign is +
                if (term1 >= 0) p[i] += p[term1];
                if (term2 >= 0) p[i] += p[term2];
            } else { // k is even, sign is -
                if (term1 >= 0) p[i] -= p[term1];
                if (term2 >= 0) p[i] -= p[term2];
            }
        }
    }
    g_data.final_result = p[n];
}

void cleanup() {
    free(g_data.p_table);
    g_data.p_table = NULL;
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
    printf("%llu\n", g_data.final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
