#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// Mersenne Twister (verbatim)
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
// End of Mersenne Twister

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

// Global data structure
typedef struct {
    int time_series_length;
    int p_order;
    int q_order;
    double* time_series;    // The observed data X_t
    double* phi_params;     // Autoregressive coefficients (p params)
    double* theta_params;   // Moving average coefficients (q params)
    double* residuals;      // Computed error terms e_t
    double final_result;    // Sum of squared residuals
} BenchmarkData;

static BenchmarkData* g_data;

// Function to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <time_series_length> <p_order> <q_order> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("Failed to allocate memory for g_data");
        exit(1);
    }

    g_data->time_series_length = atoi(argv[1]);
    g_data->p_order = atoi(argv[2]);
    g_data->q_order = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    mt_seed(seed);

    // Allocate memory for arrays
    g_data->time_series = (double*)malloc(g_data->time_series_length * sizeof(double));
    g_data->phi_params = (double*)malloc(g_data->p_order * sizeof(double));
    g_data->theta_params = (double*)malloc(g_data->q_order * sizeof(double));
    g_data->residuals = (double*)calloc(g_data->time_series_length, sizeof(double)); // Use calloc to initialize to zero

    if (!g_data->time_series || !g_data->phi_params || !g_data->theta_params || !g_data->residuals) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate synthetic time series data
    for (int i = 0; i < g_data->time_series_length; ++i) {
        g_data->time_series[i] = rand_double() * 2.0 - 1.0; // Values between -1.0 and 1.0
    }

    // Generate random ARMA parameters
    for (int i = 0; i < g_data->p_order; ++i) {
        g_data->phi_params[i] = rand_double() * 0.5 - 0.25; // Small random values
    }
    for (int i = 0; i < g_data->q_order; ++i) {
        g_data->theta_params[i] = rand_double() * 0.5 - 0.25; // Small random values
    }

    g_data->final_result = 0.0;
}

void run_computation() {
    int n = g_data->time_series_length;
    int p = g_data->p_order;
    int q = g_data->q_order;

    double* x = g_data->time_series;
    double* phi = g_data->phi_params;
    double* theta = g_data->theta_params;
    double* residuals = g_data->residuals;

    int start_index = MAX(p, q);

    // The core computation a single pass of calculating residuals in ARMA model fitting.
    // e_t = X_t - (SUM_{i=1 to p} phi_i * X_{t-i}) - (SUM_{j=1 to q} theta_j * e_{t-j})
    
    for (int t = start_index; t < n; ++t) {
        double ar_sum = 0.0;
        for (int i = 0; i < p; ++i) {
            ar_sum += phi[i] * x[t - i - 1];
        }

        double ma_sum = 0.0;
        for (int j = 0; j < q; ++j) {
            ma_sum += theta[j] * residuals[t - j - 1];
        }

        residuals[t] = x[t] - ar_sum - ma_sum;
    }

    // To prevent dead code elimination, calculate the Sum of Squared Residuals (SSR).
    double ssr = 0.0;
    for (int t = start_index; t < n; ++t) {
        ssr += residuals[t] * residuals[t];
    }
    g_data->final_result = ssr;
}

void cleanup() {
    if (g_data) {
        free(g_data->time_series);
        free(g_data->phi_params);
        free(g_data->theta_params);
        free(g_data->residuals);
        free(g_data);
    }
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data->final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
