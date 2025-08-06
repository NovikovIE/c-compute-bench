#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
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

// Helper to generate a random double from 0.0 to 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Global structure to hold benchmark data
typedef struct {
    // Parameters
    int num_umbrella_windows;
    int data_points_per_window;
    int num_bins;
    int num_iterations;

    // Data arrays for setup
    double* window_centers;         // Center of the harmonic potential for each window
    double* window_k;
    int** histograms;               // N_ij: histogram counts for window i in bin j

    // WHAM arrays for computation
    double* f_offsets;              // Free energy offset for window i
    double** V_bias;                // V_ij: bias potential of window i at center of bin j

    // Final result
    double final_pmf_checksum;
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_umbrella_windows> <data_points_per_window> <num_bins> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_umbrella_windows = atoi(argv[1]);
    g_data.data_points_per_window = atoi(argv[2]);
    g_data.num_bins = atoi(argv[3]);
    g_data.num_iterations = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);
    mt_seed(seed);

    // --- Memory Allocation ---
    g_data.window_centers = (double*)malloc(g_data.num_umbrella_windows * sizeof(double));
    g_data.window_k = (double*)malloc(g_data.num_umbrella_windows * sizeof(double));
    g_data.f_offsets = (double*)malloc(g_data.num_umbrella_windows * sizeof(double));
    
    g_data.histograms = (int**)malloc(g_data.num_umbrella_windows * sizeof(int*));
    for (int i = 0; i < g_data.num_umbrella_windows; i++) {
        g_data.histograms[i] = (int*)calloc(g_data.num_bins, sizeof(int));
    }

    g_data.V_bias = (double**)malloc(g_data.num_umbrella_windows * sizeof(double*));
    for (int i = 0; i < g_data.num_umbrella_windows; i++) {
        g_data.V_bias[i] = (double*)malloc(g_data.num_bins * sizeof(double));
    }

    // --- Data Generation ---
    double coord_min = 0.0;
    double coord_max = 100.0;
    double bin_width = (coord_max - coord_min) / g_data.num_bins;
    double fixed_force_constant = 0.2;

    // Initialize window parameters and f_offsets
    for (int i = 0; i < g_data.num_umbrella_windows; i++) {
        g_data.window_centers[i] = coord_min + (coord_max - coord_min) * (double)i / (g_data.num_umbrella_windows - 1);
        g_data.window_k[i] = fixed_force_constant;
        g_data.f_offsets[i] = 0.0;
    }

    // Pre-calculate biasing potentials V_ij for each window i and bin j
    for (int i = 0; i < g_data.num_umbrella_windows; i++) {
        for (int j = 0; j < g_data.num_bins; j++) {
            double bin_center = coord_min + (j + 0.5) * bin_width;
            double delta_x = bin_center - g_data.window_centers[i];
            g_data.V_bias[i][j] = 0.5 * g_data.window_k[i] * delta_x * delta_x;
        }
    }
    
    // Generate reaction coordinate data and build histograms
    // Using Box-Muller transform to generate normally distributed data points around window centers
    for (int i = 0; i < g_data.num_umbrella_windows; i++) {
        double mu = g_data.window_centers[i];
        double sigma = sqrt(1.0 / g_data.window_k[i]); // std dev from harmonic potential

        for (int p = 0; p < g_data.data_points_per_window; p += 2) {
            double u1 = rand_double();
            double u2 = rand_double();
            double mag = sigma * sqrt(-2.0 * log(u1));
            double z1 = mag * cos(2.0 * M_PI * u2) + mu;
            double z2 = mag * sin(2.0 * M_PI * u2) + mu;

            int bin_idx1 = (int)((z1 - coord_min) / bin_width);
            if (bin_idx1 >= 0 && bin_idx1 < g_data.num_bins) {
                g_data.histograms[i][bin_idx1]++;
            }
            if (p + 1 < g_data.data_points_per_window) {
                int bin_idx2 = (int)((z2 - coord_min) / bin_width);
                if (bin_idx2 >= 0 && bin_idx2 < g_data.num_bins) {
                    g_data.histograms[i][bin_idx2]++;
                }
            }
        }
    }
    g_data.final_pmf_checksum = 0.0;
}

void run_computation() {
    double *P_unnormalized = (double*)malloc(g_data.num_bins * sizeof(double));
    double *new_f_offsets = (double*)malloc(g_data.num_umbrella_windows * sizeof(double));
    const double tolerance = 1e-7;

    for (int iter = 0; iter < g_data.num_iterations; iter++) {
        // Step 1: Calculate unnormalized probabilities P_j for each bin
        for (int j = 0; j < g_data.num_bins; j++) {
            double numerator = 0.0;
            double denominator = 0.0;
            for (int i = 0; i < g_data.num_umbrella_windows; i++) {
                numerator += g_data.histograms[i][j];
            }

            if (numerator < tolerance) {
                P_unnormalized[j] = 0.0;
                continue;
            }

            for (int k = 0; k < g_data.num_umbrella_windows; k++) {
                denominator += (double)g_data.data_points_per_window * exp(g_data.f_offsets[k] - g_data.V_bias[k][j]);
            }
            P_unnormalized[j] = numerator / denominator;
        }

        // Step 2: Update free energy offsets f_i for each window
        for (int i = 0; i < g_data.num_umbrella_windows; i++) {
            double sum_arg = 0.0;
            for (int j = 0; j < g_data.num_bins; j++) {
                sum_arg += P_unnormalized[j] * exp(-g_data.V_bias[i][j]);
            }

            if (sum_arg < tolerance) {
                new_f_offsets[i] = 0.0; // Avoid log(0)
            } else {
                new_f_offsets[i] = -log(sum_arg);
            }
        }
        memcpy(g_data.f_offsets, new_f_offsets, g_data.num_umbrella_windows * sizeof(double));
    }

    // Final Calculation: Potential of Mean Force (PMF)
    double* F = (double*)malloc(g_data.num_bins * sizeof(double));
    double min_F = 1.0e100; 
    for (int j = 0; j < g_data.num_bins; j++) {
        if (P_unnormalized[j] > tolerance) {
            F[j] = -log(P_unnormalized[j]);
            if (F[j] < min_F) {
                min_F = F[j];
            }
        } else {
            F[j] = 1.0e100; // Effectively infinity for empty bins
        }
    }

    // Normalize PMF to have minimum of 0 and compute a checksum
    double checksum = 0.0;
    for (int j = 0; j < g_data.num_bins; j++) {
        if (F[j] < 1.0e99) {
            double normalized_F = F[j] - min_F;
            checksum += normalized_F * (j + 1); // Weighted sum to ensure value depends on order
        }
    }
    g_data.final_pmf_checksum = checksum;

    free(P_unnormalized);
    free(new_f_offsets);
    free(F);
}

void cleanup() {
    free(g_data.window_centers);
    free(g_data.window_k);
    free(g_data.f_offsets);

    for (int i = 0; i < g_data.num_umbrella_windows; i++) {
        free(g_data.histograms[i]);
        free(g_data.V_bias[i]);
    }
    free(g_data.histograms);
    free(g_data.V_bias);
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
    printf("%f\n", g_data.final_pmf_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
