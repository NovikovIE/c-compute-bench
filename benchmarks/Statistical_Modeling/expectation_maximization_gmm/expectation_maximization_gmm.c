#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

// Benchmark parameters
static int num_data_points;
static int num_features;
static int num_clusters;
static int max_iterations;
static uint32_t seed;

// Data structures for GMM
static double** data;             // [num_data_points][num_features]
static double** means;            // [num_clusters][num_features]
static double** covariances;      // [num_clusters][num_features] (diagonal)
static double* weights;           // [num_clusters]
static double** responsibilities; // [num_data_points][num_clusters]

// Final result
static double final_log_likelihood = 0.0;

// Helper to generate a random double in [0, 1]
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Helper to allocate a 2D double array
double** allocate_2d_double(int rows, int cols) {
    double** array = (double**)malloc(rows * sizeof(double*));
    if (array == NULL) return NULL;

    double* data_block = (double*)malloc((size_t)rows * cols * sizeof(double));
    if (data_block == NULL) {
        free(array);
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        array[i] = &data_block[i * cols];
    }
    return array;
}

// Helper to free a 2D double array allocated with allocate_2d_double
void free_2d_double(double** array) {
    if (array != NULL) {
        free(array[0]); // Free the contiguous block of data
        free(array);    // Free the array of pointers
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_data_points num_features num_clusters max_iterations seed\n", argv[0]);
        exit(1);
    }

    num_data_points = atoi(argv[1]);
    num_features = atoi(argv[2]);
    num_clusters = atoi(argv[3]);
    max_iterations = atoi(argv[4]);
    seed = (uint32_t)strtoul(argv[5], NULL, 10);

    mt_seed(seed);

    // Allocate memory
    data = allocate_2d_double(num_data_points, num_features);
    means = allocate_2d_double(num_clusters, num_features);
    covariances = allocate_2d_double(num_clusters, num_features);
    responsibilities = allocate_2d_double(num_data_points, num_clusters);
    weights = (double*)malloc(num_clusters * sizeof(double));

    if (!data || !means || !covariances || !responsibilities || !weights) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize data with random points
    for (int i = 0; i < num_data_points; i++) {
        for (int j = 0; j < num_features; j++) {
            data[i][j] = rand_double();
        }
    }

    // Initialize GMM parameters
    // 1. Weights: uniform
    for (int j = 0; j < num_clusters; j++) {
        weights[j] = 1.0 / num_clusters;
    }

    // 2. Means: random subset of data points
    for (int j = 0; j < num_clusters; j++) {
        int data_idx = mt_rand() % num_data_points;
        for (int f = 0; f < num_features; f++) {
            means[j][f] = data[data_idx][f];
        }
    }

    // 3. Covariances: identity (1.0 along diagonal)
    for (int j = 0; j < num_clusters; j++) {
        for (int f = 0; f < num_features; f++) {
            covariances[j][f] = 1.0;
        }
    }
}

void run_computation() {
    const double PI = 3.14159265358979323846;
    const double LOG_2PI = log(2.0 * PI);
    const double EPSILON = 1e-9; // For numerical stability

    for (int iter = 0; iter < max_iterations; ++iter) {
        double current_log_likelihood = 0.0;

        // E-Step: Calculate responsibilities
        for (int i = 0; i < num_data_points; ++i) {
            double log_probs[num_clusters]; // VLA for per-point cluster log probabilities
            double max_log_prob = -INFINITY;

            for (int j = 0; j < num_clusters; ++j) {
                double log_pdf = 0.0;
                double log_det_term = 0.0;

                for (int f = 0; f < num_features; ++f) {
                    double cov_val = covariances[j][f] + EPSILON;
                    double diff = data[i][f] - means[j][f];
                    log_pdf -= (diff * diff) / (2.0 * cov_val);
                    log_det_term += log(cov_val);
                }
                log_pdf -= 0.5 * (num_features * LOG_2PI + log_det_term);
                
                log_probs[j] = log(weights[j] + EPSILON) + log_pdf;
                if (log_probs[j] > max_log_prob) {
                    max_log_prob = log_probs[j];
                }
            }

            double sum_exp = 0.0;
            for (int j = 0; j < num_clusters; ++j) {
                sum_exp += exp(log_probs[j] - max_log_prob);
            }

            double log_sum = max_log_prob + log(sum_exp);
            current_log_likelihood += log_sum;

            for (int j = 0; j < num_clusters; ++j) {
                responsibilities[i][j] = exp(log_probs[j] - log_sum);
            }
        }

        // M-Step: Update parameters
        for (int j = 0; j < num_clusters; ++j) {
            double soft_count = 0.0;
            for (int i = 0; i < num_data_points; ++i) {
                soft_count += responsibilities[i][j];
            }

            weights[j] = soft_count / num_data_points;
            
            double effective_soft_count = soft_count + EPSILON;

            for (int f = 0; f < num_features; ++f) {
                double mean_sum = 0.0;
                for (int i = 0; i < num_data_points; ++i) {
                    mean_sum += responsibilities[i][j] * data[i][f];
                }
                means[j][f] = mean_sum / effective_soft_count;
            }

            for (int f = 0; f < num_features; ++f) {
                double cov_sum = 0.0;
                for (int i = 0; i < num_data_points; ++i) {
                    double diff = data[i][f] - means[j][f];
                    cov_sum += responsibilities[i][j] * diff * diff;
                }
                covariances[j][f] = cov_sum / effective_soft_count;
            }
        }

        final_log_likelihood = current_log_likelihood;
    }
}

void cleanup() {
    free_2d_double(data);
    free_2d_double(means);
    free_2d_double(covariances);
    free_2d_double(responsibilities);
    free(weights);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_log_likelihood);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
