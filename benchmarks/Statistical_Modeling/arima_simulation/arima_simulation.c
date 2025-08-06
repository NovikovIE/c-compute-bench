#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- MERSENNE TWISTER (Verbatim) ---
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

// --- BENCHMARK DATA AND FUNCTIONS ---

typedef struct {
    int num_time_steps;
    int p_order;
    int d_order;
    int q_order;
    double* ar_coeffs;
    double* ma_coeffs;
    double* noise;
    double* series;
    double final_sum;
} BenchmarkData;

BenchmarkData g_data;

// Function to generate a random double between -0.5 and 0.5
double random_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) - 0.5;
}

// Function to generate a random double between 0.0 and 1.0 (for noise)
double random_noise() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_time_steps> <p_order> <d_order> <q_order> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_time_steps = atoi(argv[1]);
    g_data.p_order = atoi(argv[2]);
    g_data.d_order = atoi(argv[3]);
    g_data.q_order = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    if(g_data.num_time_steps <= 0 || g_data.p_order < 0 || g_data.d_order < 0 || g_data.q_order < 0) {
        fprintf(stderr, "FATAL: Invalid parameters. All orders must be non-negative and steps must be positive.\n");
        exit(1);
    }
    
    int max_lag = g_data.p_order > g_data.q_order ? g_data.p_order : g_data.q_order;
    if (g_data.num_time_steps <= max_lag) {
        fprintf(stderr, "FATAL: num_time_steps must be greater than max(p_order, q_order).\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.ar_coeffs = (double*)malloc(g_data.p_order * sizeof(double));
    g_data.ma_coeffs = (double*)malloc(g_data.q_order * sizeof(double));
    g_data.noise = (double*)malloc(g_data.num_time_steps * sizeof(double));
    g_data.series = (double*)malloc(g_data.num_time_steps * sizeof(double));

    if (!g_data.ar_coeffs || !g_data.ma_coeffs || !g_data.noise || !g_data.series) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    
    for (int i = 0; i < g_data.p_order; ++i) {
        g_data.ar_coeffs[i] = random_double();
    }
    
    for (int i = 0; i < g_data.q_order; ++i) {
        g_data.ma_coeffs[i] = random_double();
    }
    
    for (int i = 0; i < g_data.num_time_steps; ++i) {
        g_data.noise[i] = random_noise();
    }

    g_data.final_sum = 0.0;
}

void run_computation() {
    int p = g_data.p_order;
    int q = g_data.q_order;
    int d = g_data.d_order;
    int n = g_data.num_time_steps;
    
    double* series = g_data.series;
    double* noise = g_data.noise;
    double* ar_coeffs = g_data.ar_coeffs;
    double* ma_coeffs = g_data.ma_coeffs;

    int max_lag = p > q ? p : q;

    // Generate ARMA(p, q) series
    // Initialize first max_lag values
    for (int t = 0; t < max_lag; ++t) {
        series[t] = noise[t];
    }
    
    // Generate the rest of the series
    for (int t = max_lag; t < n; ++t) {
        double ar_sum = 0.0;
        for (int i = 0; i < p; ++i) {
            ar_sum += ar_coeffs[i] * series[t - 1 - i];
        }

        double ma_sum = 0.0;
        for (int i = 0; i < q; ++i) {
            ma_sum += ma_coeffs[i] * noise[t - 1 - i];
        }

        series[t] = ar_sum + ma_sum + noise[t];
    }

    // Integrate the series d times (in-place cumulative sum)
    for (int i = 0; i < d; ++i) {
        for (int t = 1; t < n; ++t) {
            series[t] += series[t - 1];
        }
    }

    // Calculate a checksum to prevent dead code elimination
    double sum = 0.0;
    for (int t = 0; t < n; ++t) {
        sum += series[t];
    }
    g_data.final_sum = sum;
}

void cleanup() {
    free(g_data.ar_coeffs);
    free(g_data.ma_coeffs);
    free(g_data.noise);
    free(g_data.series);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_sum);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
