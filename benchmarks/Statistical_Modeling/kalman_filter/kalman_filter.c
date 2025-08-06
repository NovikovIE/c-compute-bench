#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
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

// --- Global Benchmark Variables ---
int num_time_steps;
int n; // num_state_variables
int m; // num_observation_variables

double *x_hat; // State estimate (n x 1)
double *P;     // Error covariance (n x n)
double *F;     // State transition matrix (n x n)
double *Q;     // Process noise covariance (n x n)
double *H;     // Observation matrix (m x n)
double *R;     // Measurement noise covariance (m x m)
double *z_all; // All observation vectors (num_time_steps * m)

double result_accumulator;

// --- Helper Functions ---
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void matrix_multiply(double *A, double *B, double *C, int rA, int cA, int cB) {
    for (int i = 0; i < rA; ++i) {
        for (int j = 0; j < cB; ++j) {
            double sum = 0.0;
            for (int k = 0; k < cA; ++k) {
                sum += A[i * cA + k] * B[k * cB + j];
            }
            C[i * cB + j] = sum;
        }
    }
}

void matrix_add(double *A, double *B, double *C, int r, int c) {
    for (int i = 0; i < r * c; ++i) C[i] = A[i] + B[i];
}

void matrix_subtract(double *A, double *B, double *C, int r, int c) {
    for (int i = 0; i < r * c; ++i) C[i] = A[i] - B[i];
}

void matrix_transpose(double *A, double *AT, int r, int c) {
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            AT[j * r + i] = A[i * c + j];
        }
    }
}

int matrix_inverse(double *A, double *A_inv, int size) {
    int augmented_cols = 2 * size;
    double *aug = (double *)malloc(size * augmented_cols * sizeof(double));
    if (!aug) return 1;

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            aug[i * augmented_cols + j] = A[i * size + j];
            aug[i * augmented_cols + size + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < size; ++i) {
        int pivot_row = i;
        for (int k = i + 1; k < size; ++k) {
            if (fabs(aug[k * augmented_cols + i]) > fabs(aug[pivot_row * augmented_cols + i])) {
                pivot_row = k;
            }
        }
        if (pivot_row != i) {
            for (int j = 0; j < augmented_cols; ++j) {
                double temp = aug[i * augmented_cols + j];
                aug[i * augmented_cols + j] = aug[pivot_row * augmented_cols + j];
                aug[pivot_row * augmented_cols + j] = temp;
            }
        }
        if (fabs(aug[i * augmented_cols + i]) < 1e-12) {
            free(aug);
            return 1;
        }
        double pivot_val = aug[i * augmented_cols + i];
        for (int j = i; j < augmented_cols; ++j) {
            aug[i * augmented_cols + j] /= pivot_val;
        }
        for (int k = 0; k < size; ++k) {
            if (k != i) {
                double factor = aug[k * augmented_cols + i];
                for (int j = i; j < augmented_cols; ++j) {
                    aug[k * augmented_cols + j] -= factor * aug[i * augmented_cols + j];
                }
            }
        }
    }

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            A_inv[i * size + j] = aug[i * augmented_cols + size + j];
        }
    }
    free(aug);
    return 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_time_steps> <num_state> <num_obs> <seed>\n", argv[0]);
        exit(1);
    }

    num_time_steps = atoi(argv[1]);
    n = atoi(argv[2]);
    m = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    x_hat = (double*)malloc(n * sizeof(double));
    P = (double*)malloc(n * n * sizeof(double));
    F = (double*)malloc(n * n * sizeof(double));
    Q = (double*)malloc(n * n * sizeof(double));
    H = (double*)malloc(m * n * sizeof(double));
    R = (double*)malloc(m * m * sizeof(double));
    z_all = (double*)malloc(num_time_steps * m * sizeof(double));

    for (int i = 0; i < n; i++) x_hat[i] = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            F[i * n + j] = (i == j) ? 1.0 : rand_double() * 0.1;
            P[i * n + j] = (i == j) ? 1.0 : 0.0;
            Q[i * n + j] = (i == j) ? (rand_double() * 0.1 + 0.01) : 0.0;
        }
    }
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            R[i * m + j] = (i == j) ? (rand_double() * 0.5 + 0.1) : 0.0;
        }
    }
    for (int i = 0; i < m * n; i++) H[i] = rand_double();
    for (int i = 0; i < num_time_steps * m; i++) z_all[i] = rand_double() * 10.0;

    result_accumulator = 0.0;
}

void run_computation() {
    double *x_hat_pred = (double*)malloc(n * sizeof(double));
    double *P_pred = (double*)malloc(n * n * sizeof(double));
    double *F_T = (double*)malloc(n * n * sizeof(double));
    double *H_T = (double*)malloc(n * m * sizeof(double));
    double *y = (double*)malloc(m * sizeof(double));
    double *S = (double*)malloc(m * m * sizeof(double));
    double *S_inv = (double*)malloc(m * m * sizeof(double));
    double *K = (double*)malloc(n * m * sizeof(double));
    double *I_n = (double*)malloc(n * n * sizeof(double));
    double *temp_n_n_1 = (double*)malloc(n * n * sizeof(double));
    double *temp_n_n_2 = (double*)malloc(n * n * sizeof(double));
    double *temp_n_m = (double*)malloc(n * m * sizeof(double));
    double *temp_m_n = (double*)malloc(m * n * sizeof(double));
    double *temp_m_m = (double*)malloc(m * m * sizeof(double));
    double *temp_n_1 = (double*)malloc(n * sizeof(double));
    double *temp_m_1 = (double*)malloc(m * sizeof(double));

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            I_n[i*n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    matrix_transpose(F, F_T, n, n);
    matrix_transpose(H, H_T, m, n);

    for (int t = 0; t < num_time_steps; ++t) {
        // Prediction
        matrix_multiply(F, x_hat, x_hat_pred, n, n, 1);
        matrix_multiply(F, P, temp_n_n_1, n, n, n);
        matrix_multiply(temp_n_n_1, F_T, P_pred, n, n, n);
        matrix_add(P_pred, Q, P_pred, n, n);

        // Update
        double *z_k = &z_all[t * m];
        matrix_multiply(H, x_hat_pred, temp_m_1, m, n, 1);
        matrix_subtract(z_k, temp_m_1, y, m, 1);
        matrix_multiply(H, P_pred, temp_m_n, m, n, n);
        matrix_multiply(temp_m_n, H_T, temp_m_m, m, n, m);
        matrix_add(temp_m_m, R, S, m, m);

        if (matrix_inverse(S, S_inv, m) != 0) {
            for(int i = 0; i < n * n; i++) P[i] = P_pred[i];
            for(int i = 0; i < n; i++) x_hat[i] = x_hat_pred[i];
            continue;
        }

        matrix_multiply(P_pred, H_T, temp_n_m, n, n, m);
        matrix_multiply(temp_n_m, S_inv, K, n, m, m);
        matrix_multiply(K, y, temp_n_1, n, m, 1);
        matrix_add(x_hat_pred, temp_n_1, x_hat, n, 1);
        matrix_multiply(K, H, temp_n_n_1, n, m, n);
        matrix_subtract(I_n, temp_n_n_1, temp_n_n_2, n, n);
        matrix_multiply(temp_n_n_2, P_pred, P, n, n, n);

        result_accumulator += x_hat[0];
    }

    free(x_hat_pred); free(P_pred); free(F_T); free(H_T); free(y); free(S);
    free(S_inv); free(K); free(I_n); free(temp_n_n_1); free(temp_n_n_2);
    free(temp_n_m); free(temp_m_n); free(temp_m_m); free(temp_n_1); free(temp_m_1);
}

void cleanup() {
    free(x_hat);
    free(P);
    free(F);
    free(Q);
    free(H);
    free(R);
    free(z_all);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", result_accumulator);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
