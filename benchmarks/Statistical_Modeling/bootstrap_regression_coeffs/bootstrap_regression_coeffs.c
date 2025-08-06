#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- BEGIN MERSENNE TWISTER (MT19937) ---
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
// --- END MERSENNE TWISTER ---

typedef struct {
    int num_data_points;
    int num_predictors;
    int num_bootstrap_samples;
    double *X_data; // Design matrix, size: num_data_points x num_predictors
    double *y_data; // Response vector, size: num_data_points
    double result_accumulator;
} BenchmarkData;

BenchmarkData g_data;

void solve_linear_system(double *A, double *b, double *x, int n);

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_data_points> <num_predictors> <num_bootstrap_samples> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_data_points = atoi(argv[1]);
    g_data.num_predictors = atoi(argv[2]);
    g_data.num_bootstrap_samples = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    g_data.X_data = (double *)malloc((size_t)g_data.num_data_points * g_data.num_predictors * sizeof(double));
    g_data.y_data = (double *)malloc((size_t)g_data.num_data_points * sizeof(double));
    double *beta_true = (double *)malloc((size_t)g_data.num_predictors * sizeof(double));

    if (!g_data.X_data || !g_data.y_data || !beta_true) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_predictors; ++i) {
        beta_true[i] = ((double)mt_rand() / UINT32_MAX) * 2.0 - 1.0; // in [-1, 1]
    }

    for (int i = 0; i < g_data.num_data_points; ++i) {
        double y_val = 0.0;
        for (int j = 0; j < g_data.num_predictors; ++j) {
            double u1 = (double)mt_rand() / (UINT32_MAX + 1.0);
            double u2 = (double)mt_rand() / (UINT32_MAX + 1.0);
            double rand_std_normal = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2); // Box-Muller transform
            g_data.X_data[i * g_data.num_predictors + j] = rand_std_normal;
            y_val += g_data.X_data[i * g_data.num_predictors + j] * beta_true[j];
        }
        double noise = ((double)mt_rand() / UINT32_MAX - 0.5) * 0.2;
        g_data.y_data[i] = y_val + noise;
    }

    free(beta_true);
    g_data.result_accumulator = 0.0;
}

void solve_linear_system(double *A, double *b, double *x, int n) {
    double *aug_matrix = (double *)malloc((size_t)n * (n + 1) * sizeof(double));
    if (!aug_matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed in solver.\n");
        exit(1);
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            aug_matrix[i * (n + 1) + j] = A[i * n + j];
        }
        aug_matrix[i * (n + 1) + n] = b[i];
    }

    for (int k = 0; k < n; ++k) {
        int max_row = k;
        for (int i = k + 1; i < n; ++i) {
            if (fabs(aug_matrix[i * (n + 1) + k]) > fabs(aug_matrix[max_row * (n + 1) + k])) {
                max_row = i;
            }
        }
        if (max_row != k) {
            for (int j = k; j < n + 1; ++j) {
                double temp = aug_matrix[k * (n + 1) + j];
                aug_matrix[k * (n + 1) + j] = aug_matrix[max_row * (n + 1) + j];
                aug_matrix[max_row * (n + 1) + j] = temp;
            }
        }

        double pivot = aug_matrix[k * (n + 1) + k];
        if (fabs(pivot) < 1e-12) {
            free(aug_matrix);
            for(int i = 0; i < n; ++i) x[i] = 0.0;
            return;
        }

        for (int i = k + 1; i < n; ++i) {
            double factor = aug_matrix[i * (n + 1) + k] / pivot;
            for (int j = k; j < n + 1; ++j) {
                aug_matrix[i * (n + 1) + j] -= factor * aug_matrix[k * (n + 1) + j];
            }
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += aug_matrix[i * (n + 1) + j] * x[j];
        }
        x[i] = (aug_matrix[i * (n + 1) + n] - sum) / aug_matrix[i * (n + 1) + i];
    }

    free(aug_matrix);
}

void run_computation() {
    int n = g_data.num_data_points;
    int k = g_data.num_predictors;

    double *XtX = (double *)malloc((size_t)k * k * sizeof(double));
    double *Xty = (double *)malloc((size_t)k * sizeof(double));
    double *beta_hat = (double *)malloc((size_t)k * sizeof(double));
    int *indices = (int *)malloc((size_t)n * sizeof(int));

    if (!XtX || !Xty || !beta_hat || !indices) {
        fprintf(stderr, "FATAL: Memory allocation failed in computation.\n");
        exit(1);
    }

    for (int b = 0; b < g_data.num_bootstrap_samples; ++b) {
        for (int i = 0; i < n; ++i) {
            indices[i] = mt_rand() % n;
        }

        for (int i = 0; i < k; ++i) {
            for (int j = i; j < k; ++j) {
                double sum = 0.0;
                for (int row = 0; row < n; ++row) {
                    int sample_idx = indices[row];
                    sum += g_data.X_data[sample_idx * k + i] * g_data.X_data[sample_idx * k + j];
                }
                XtX[i * k + j] = sum;
                XtX[j * k + i] = sum;
            }
            double y_sum = 0.0;
            for (int row = 0; row < n; ++row) {
                int sample_idx = indices[row];
                y_sum += g_data.X_data[sample_idx * k + i] * g_data.y_data[sample_idx];
            }
            Xty[i] = y_sum;
        }

        solve_linear_system(XtX, Xty, beta_hat, k);

        for (int i = 0; i < k; ++i) {
            g_data.result_accumulator += beta_hat[i];
        }
    }

    free(XtX);
    free(Xty);
    free(beta_hat);
    free(indices);
}

void cleanup() {
    free(g_data.X_data);
    free(g_data.y_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.result_accumulator);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
