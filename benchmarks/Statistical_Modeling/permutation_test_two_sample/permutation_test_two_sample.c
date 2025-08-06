#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator ---
// Do not modify
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
int g_sample1_size;
int g_sample2_size;
int g_num_permutations;
int g_total_size;

double *g_sample1;
double *g_sample2;
double *g_pooled_data;
double *g_shuffled_buffer;

int g_final_result;
// --- End of Globals ---

// Function to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <sample1_size> <sample2_size> <num_permutations> <seed>\n", argv[0]);
        exit(1);
    }

    g_sample1_size = atoi(argv[1]);
    g_sample2_size = atoi(argv[2]);
    g_num_permutations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    g_total_size = g_sample1_size + g_sample2_size;

    g_sample1 = (double*)malloc(g_sample1_size * sizeof(double));
    g_sample2 = (double*)malloc(g_sample2_size * sizeof(double));
    g_pooled_data = (double*)malloc(g_total_size * sizeof(double));
    g_shuffled_buffer = (double*)malloc(g_total_size * sizeof(double));

    if (!g_sample1 || !g_sample2 || !g_pooled_data || !g_shuffled_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate two samples from slightly different distributions
    // to ensure a non-zero observed difference.
    for (int i = 0; i < g_sample1_size; ++i) {
        g_sample1[i] = rand_double();
    }
    for (int i = 0; i < g_sample2_size; ++i) {
        g_sample2[i] = rand_double() + 0.05; // Small effect size
    }

    // Pool the data
    memcpy(g_pooled_data, g_sample1, g_sample1_size * sizeof(double));
    memcpy(g_pooled_data + g_sample1_size, g_sample2, g_sample2_size * sizeof(double));
}

void run_computation() {
    // 1. Calculate the observed difference in means
    double sum1_obs = 0.0;
    for (int i = 0; i < g_sample1_size; ++i) {
        sum1_obs += g_sample1[i];
    }
    double mean1_obs = sum1_obs / g_sample1_size;

    double sum2_obs = 0.0;
    for (int i = 0; i < g_sample2_size; ++i) {
        sum2_obs += g_sample2[i];
    }
    double mean2_obs = sum2_obs / g_sample2_size;
    
    double observed_diff = fabs(mean1_obs - mean2_obs);

    int more_extreme_count = 0;

    // 2. Perform permutations
    for (int p = 0; p < g_num_permutations; ++p) {
        // Create a fresh copy of pooled data to shuffle
        memcpy(g_shuffled_buffer, g_pooled_data, g_total_size * sizeof(double));

        // Shuffle the buffer using Fisher-Yates algorithm
        for (int i = g_total_size - 1; i > 0; --i) {
            int j = mt_rand() % (i + 1);
            double temp = g_shuffled_buffer[i];
            g_shuffled_buffer[i] = g_shuffled_buffer[j];
            g_shuffled_buffer[j] = temp;
        }

        // Calculate means of the permuted groups
        double sum1_perm = 0.0;
        for (int i = 0; i < g_sample1_size; ++i) {
            sum1_perm += g_shuffled_buffer[i];
        }
        double mean1_perm = sum1_perm / g_sample1_size;

        double sum2_perm = 0.0;
        // The second group starts after the first one ends
        for (int i = g_sample1_size; i < g_total_size; ++i) {
            sum2_perm += g_shuffled_buffer[i];
        }
        double mean2_perm = sum2_perm / g_sample2_size;
        
        double permuted_diff = fabs(mean1_perm - mean2_perm);

        // 3. Count how many times the permuted diff is >= observed diff
        if (permuted_diff >= observed_diff) {
            more_extreme_count++;
        }
    }
    g_final_result = more_extreme_count;
}

void cleanup() {
    free(g_sample1);
    free(g_sample2);
    free(g_pooled_data);
    free(g_shuffled_buffer);
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
    printf("%d\n", g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
