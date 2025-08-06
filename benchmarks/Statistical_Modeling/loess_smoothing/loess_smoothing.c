#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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

// Global structure to hold benchmark data
struct BenchmarkData {
    int num_data_points;
    double span;
    double *x;          // Independent variable
    double *y;          // Dependent variable (with noise)
    double *y_smooth;   // Smoothed output from LOESS
};

struct BenchmarkData g_data;
double g_result = 0.0;

// Setup: parse arguments, allocate memory, generate data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_data_points> <span> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_data_points = atoi(argv[1]);
    g_data.span = atof(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_data_points <= 0 || g_data.span <= 0.0 || g_data.span > 1.0) {
        fprintf(stderr, "Invalid arguments.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.x = (double *)malloc(g_data.num_data_points * sizeof(double));
    g_data.y = (double *)malloc(g_data.num_data_points * sizeof(double));
    g_data.y_smooth = (double *)malloc(g_data.num_data_points * sizeof(double));

    if (!g_data.x || !g_data.y || !g_data.y_smooth) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate a simple sine wave with random noise
    for (int i = 0; i < g_data.num_data_points; i++) {
        g_data.x[i] = (double)i;
        double noise = (mt_rand() / (double)UINT32_MAX - 0.5) * 0.5;
        g_data.y[i] = sin(g_data.x[i] / (g_data.num_data_points / 30.0)) + noise;
    }
}

// Computation: perform LOESS smoothing
void run_computation() {
    int n = g_data.num_data_points;
    // Calculate window size based on span
    int k = (int)(g_data.span * n);
    if (k < 2) k = 2; // Need at least 2 points for a linear fit

    for (int i = 0; i < n; i++) {
        // Define the neighborhood window for the current point x[i]
        // This is a simplified approach for sorted x data
        int window_start = (i > k / 2) ? (i - k / 2) : 0;
        int window_end = window_start + k;
        if (window_end > n) {
            window_end = n;
            window_start = n - k;
        }

        // Find the maximum distance in the window for tricube weight calculation
        double x_i = g_data.x[i];
        double dist1 = x_i - g_data.x[window_start];
        double dist2 = g_data.x[window_end - 1] - x_i;
        double max_dist = (dist1 > dist2) ? dist1 : dist2;

        if (max_dist < 1e-9) {
            g_data.y_smooth[i] = g_data.y[i]; // Fallback for single-point window
            continue;
        }

        // Perform weighted linear regression for the neighborhood
        double Sw = 0, Swx = 0, Swy = 0, Swxx = 0, Swxy = 0;
        for (int j = window_start; j < window_end; j++) {
            double x_j = g_data.x[j];
            double y_j = g_data.y[j];
            double dist = fabs(x_i - x_j);
            double u = dist / max_dist;
            double w_val = 1.0 - u * u * u;
            double w = w_val * w_val * w_val; // Tricube weight

            Sw   += w;
            Swx  += w * x_j;
            Swy  += w * y_j;
            Swxx += w * x_j * x_j;
            Swxy += w * x_j * y_j;
        }

        // Calculate coefficients a (intercept) and b (slope)
        double a, b;
        double denom = Sw * Swxx - Swx * Swx;
        if (fabs(denom) > 1e-12) {
            b = (Sw * Swxy - Swx * Swy) / denom;
            a = (Swy - b * Swx) / Sw;
        } else {
            // Fallback for a singular matrix (e.g., all x in window are the same)
            b = 0;
            a = g_data.y[i]; 
        }
        
        g_data.y_smooth[i] = a + b * x_i;
    }

    // Accumulate a final result to prevent dead code elimination
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += g_data.y_smooth[i];
    }
    g_result = sum;
}

// Cleanup: free allocated memory
void cleanup() {
    free(g_data.x);
    free(g_data.y);
    free(g_data.y_smooth);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%f\n", g_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}