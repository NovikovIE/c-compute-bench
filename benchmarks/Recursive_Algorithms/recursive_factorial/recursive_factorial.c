#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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
// --- End of MT19937 ---

// --- Benchmark Globals ---
typedef struct {
    int n;
    long long iterations;
    unsigned long long result;
} BenchmarkData;

BenchmarkData *g_data;

// --- Core Recursive Function ---
static unsigned long long recursive_factorial(int n) {
    if (n <= 1) {
        return 1;
    }
    return (unsigned long long)n * recursive_factorial(n - 1);
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <n_value> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData *)malloc(sizeof(BenchmarkData));
    if (g_data == NULL) {
        fprintf(stderr, "Failed to allocate memory for benchmark data\n");
        exit(1);
    }

    g_data->n = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed); // Seed the generator, although it's not used in this specific benchmark

    // Parameter validation
    if (g_data->n < 0 || g_data->n > 20) {
        fprintf(stderr, "Error: n_value must be between 0 and 20 to avoid unsigned long long overflow.\n");
        free(g_data);
        exit(1);
    }

    // Fixed number of iterations to make the benchmark run for a meaningful duration.
    // Total recursive calls = n * iterations. For n=19, this is 19 * 25,000,000 = 475M calls.
    // This should take approximately 1 second on a modern CPU core.
    g_data->iterations = 25000000;
    g_data->result = 0;
}

void run_computation() {
    unsigned long long total_sum = 0;
    int n = g_data->n;
    long long iterations = g_data->iterations;

    for (long long i = 0; i < iterations; i++) {
        total_sum += recursive_factorial(n);
    }

    g_data->result = total_sum;
}

void cleanup() {
    free(g_data);
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print result to stdout to prevent dead code elimination & provide a check value
    printf("%llu\n", g_data->result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
