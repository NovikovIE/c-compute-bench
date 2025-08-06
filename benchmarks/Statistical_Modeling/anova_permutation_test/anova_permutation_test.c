#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) --- DO NOT MODIFY ---
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
    int num_groups;
    int num_samples_per_group;
    int num_permutations;
    int total_samples;

    double *original_data; // N = num_groups * num_samples_per_group
    double *permuted_data; // A buffer for shuffling
    double *group_means;   // A buffer for calculating means

    int significant_permutations_count; // The final result
} BenchmarkData;

BenchmarkData g_data;

// --- Helper Functions ---

double calculate_f_statistic(const double* data) {
    int k = g_data.num_groups;
    int n_per_group = g_data.num_samples_per_group;
    int N = g_data.total_samples;

    if (k <= 1 || N <= k) return 0.0;

    // 1. Calculate grand mean
    double grand_mean = 0.0;
    for (int i = 0; i < N; ++i) {
        grand_mean += data[i];
    }
    grand_mean /= N;

    // 2. Calculate group means and Sum of Squares Within (SSW)
    double ssw = 0.0;
    for (int i = 0; i < k; ++i) {
        double group_sum = 0.0;
        for (int j = 0; j < n_per_group; ++j) {
            group_sum += data[i * n_per_group + j];
        }
        g_data.group_means[i] = group_sum / n_per_group;
        
        for (int j = 0; j < n_per_group; ++j) {
            double diff = data[i * n_per_group + j] - g_data.group_means[i];
            ssw += diff * diff;
        }
    }

    // 3. Calculate Sum of Squares Between (SSB)
    double ssb = 0.0;
    for (int i = 0; i < k; ++i) {
        double diff = g_data.group_means[i] - grand_mean;
        ssb += n_per_group * (diff * diff);
    }
    
    // 4. Calculate F-statistic
    double msw = ssw / (N - k);
    double msb = ssb / (k - 1);

    if (msw < 1e-12) {
        return (msb > 1e-12) ? 1e12 : 0.0;
    }

    return msb / msw;
}

void shuffle(double *array, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = mt_rand() % (i + 1);
        double temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_groups num_samples_per_group num_permutations seed\n", argv[0]);
        exit(1);
    }

    g_data.num_groups = atoi(argv[1]);
    g_data.num_samples_per_group = atoi(argv[2]);
    g_data.num_permutations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    g_data.total_samples = g_data.num_groups * g_data.num_samples_per_group;
    g_data.significant_permutations_count = 0;

    g_data.original_data = (double*)malloc(g_data.total_samples * sizeof(double));
    g_data.permuted_data = (double*)malloc(g_data.total_samples * sizeof(double));
    g_data.group_means = (double*)malloc(g_data.num_groups * sizeof(double));

    if (!g_data.original_data || !g_data.permuted_data || !g_data.group_means) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate data where different groups have different means to make the F-statistic meaningful
    for (int i = 0; i < g_data.num_groups; ++i) {
        for (int j = 0; j < g_data.num_samples_per_group; ++j) {
            double val = (double)mt_rand() / (double)UINT32_MAX + (double)i;
            g_data.original_data[i * g_data.num_samples_per_group + j] = val;
        }
    }
}

void run_computation() {
    double observed_f = calculate_f_statistic(g_data.original_data);

    for (int i = 0; i < g_data.total_samples; ++i) {
        g_data.permuted_data[i] = g_data.original_data[i];
    }

    for (int p = 0; p < g_data.num_permutations; ++p) {
        shuffle(g_data.permuted_data, g_data.total_samples);
        double permuted_f = calculate_f_statistic(g_data.permuted_data);
        if (permuted_f >= observed_f) {
            g_data.significant_permutations_count++;
        }
    }
}

void cleanup() {
    free(g_data.original_data);
    free(g_data.permuted_data);
    free(g_data.group_means);
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
    printf("%d\n", g_data.significant_permutations_count);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
