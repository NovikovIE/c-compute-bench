#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End Mersenne Twister ---

// Benchmark data structure
typedef struct {
    int num_points;
    int num_interp_points;
    double *x;
    double *y;
    double *interp_x;
    double *interp_y;
    double *y2; // Stores the second derivatives
    double final_result; // Accumulated result
} BenchmarkData;

static BenchmarkData *g_data;

// Generates a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <num_interp_points> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("Failed to allocate memory for BenchmarkData");
        exit(1);
    }

    g_data->num_points = atoi(argv[1]);
    g_data->num_interp_points = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (g_data->num_points < 3) {
        fprintf(stderr, "Error: num_points must be at least 3 for cubic spline.\n");
        exit(1);
    }

    // Allocate memory for data arrays
    g_data->x = (double*)malloc(g_data->num_points * sizeof(double));
    g_data->y = (double*)malloc(g_data->num_points * sizeof(double));
    g_data->y2 = (double*)malloc(g_data->num_points * sizeof(double));
    g_data->interp_x = (double*)malloc(g_data->num_interp_points * sizeof(double));
    g_data->interp_y = (double*)malloc(g_data->num_interp_points * sizeof(double));

    if (!g_data->x || !g_data->y || !g_data->y2 || !g_data->interp_x || !g_data->interp_y) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        free(g_data);
        exit(1);
    }

    // Generate sorted x points and corresponding random y points
    for (int i = 0; i < g_data->num_points; i++) {
        g_data->x[i] = (double)i; // Simple, sorted x values for easy interval finding
        g_data->y[i] = rand_double() * 200.0 - 100.0; // y values in [-100, 100]
    }

    // Generate random points at which to interpolate
    double x_max = g_data->x[g_data->num_points - 1];
    for (int i = 0; i < g_data->num_interp_points; i++) {
        g_data->interp_x[i] = rand_double() * x_max;
    }

    g_data->final_result = 0.0;
}

// Helper to find the interval k such that x[k] <= val < x[k+1]
static int find_interval(double val, const int n, const double* x_arr) {
    int klo = 0;
    int khi = n - 1;
    while (khi - klo > 1) {
        int k = (khi + klo) >> 1;
        if (x_arr[k] > val) {
            khi = k;
        } else {
            klo = k;
        }
    }
    return klo;
}

void run_computation() {
    int n = g_data->num_points;
    double *x = g_data->x;
    double *y = g_data->y;
    double *y2 = g_data->y2;
    
    // Step 1: Compute the second derivatives (y2) using a tridiagonal solver
    // This implements a natural cubic spline (y2[0] = y2[n-1] = 0)
    double *u = (double*)malloc((n - 1) * sizeof(double));
    if (!u) {
        perror("Failed to allocate temporary memory in run_computation");
        exit(1);
    }
    
    y2[0] = 0.0;
    u[0] = 0.0;

    for (int i = 1; i < n - 1; i++) {
        double sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
        double p = sig * y2[i-1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        u[i] = (y[i+1] - y[i]) / (x[i+1] - x[i]) - (y[i] - y[i-1]) / (x[i] - x[i-1]);
        u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p;
    }

    y2[n-1] = 0.0;
    for (int k = n - 2; k >= 0; k--) {
        y2[k] = y2[k] * y2[k+1] + u[k];
    }

    free(u);

    // Step 2: Perform interpolation for each point in interp_x
    double total_sum = 0.0;
    for (int i = 0; i < g_data->num_interp_points; i++) {
        double px = g_data->interp_x[i];

        int k = find_interval(px, n, x);
        
        double h = x[k+1] - x[k];
        double a = (x[k+1] - px) / h;
        double b = (px - x[k]) / h;
        g_data->interp_y[i] = a * y[k] + b * y[k+1] + 
                              ((a*a*a - a) * y2[k] + (b*b*b - b) * y2[k+1]) * (h * h) / 6.0;
        total_sum += g_data->interp_y[i];
    }
    
    g_data->final_result = total_sum;
}

void cleanup() {
    free(g_data->x);
    free(g_data->y);
    free(g_data->y2);
    free(g_data->interp_x);
    free(g_data->interp_y);
    free(g_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print final result to stdout
    printf("%f\n", g_data->final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
