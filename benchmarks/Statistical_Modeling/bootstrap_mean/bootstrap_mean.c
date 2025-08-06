#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Benchmark Globals ---
typedef struct {
    int num_data_points;
    int num_bootstrap_samples;
    uint32_t seed;
    double *data_points;
    double *bootstrap_means;
    double final_result_sum;
} BenchmarkData;

BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_data_points> <num_bootstrap_samples> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_data_points = atoi(argv[1]);
    g_data.num_bootstrap_samples = atoi(argv[2]);
    g_data.seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_data_points <= 0 || g_data.num_bootstrap_samples <= 0) {
        fprintf(stderr, "FATAL: Parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(g_data.seed);

    g_data.data_points = (double*)malloc(g_data.num_data_points * sizeof(double));
    if (g_data.data_points == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for data_points.\n");
        exit(1);
    }

    g_data.bootstrap_means = (double*)malloc(g_data.num_bootstrap_samples * sizeof(double));
    if (g_data.bootstrap_means == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for bootstrap_means.\n");
        free(g_data.data_points);
        exit(1);
    }

    // Generate initial data from a pseudo-normal distribution (sum of uniform)
    for (int i = 0; i < g_data.num_data_points; ++i) {
        double sum = 0.0;
        for(int j = 0; j < 12; ++j) { // Central Limit Theorem approximation
            sum += (double)mt_rand() / UINT32_MAX;
        }
        g_data.data_points[i] = sum - 6.0;
    }

    g_data.final_result_sum = 0.0;
}

void run_computation() {
    for (int i = 0; i < g_data.num_bootstrap_samples; ++i) {
        double current_sample_sum = 0.0;
        // Create a bootstrap sample by resampling with replacement
        for (int j = 0; j < g_data.num_data_points; ++j) {
            int random_index = mt_rand() % g_data.num_data_points;
            current_sample_sum += g_data.data_points[random_index];
        }
        double mean = current_sample_sum / g_data.num_data_points;
        g_data.bootstrap_means[i] = mean;
    }

    // Accumulate results into a single value to prevent dead code elimination
    double total_sum = 0.0;
    for (int i = 0; i < g_data.num_bootstrap_samples; ++i) {
        total_sum += g_data.bootstrap_means[i];
    }
    g_data.final_result_sum = total_sum;
}

void cleanup() {
    free(g_data.data_points);
    free(g_data.bootstrap_means);
    g_data.data_points = NULL;
    g_data.bootstrap_means = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", g_data.final_result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
