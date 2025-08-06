#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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


// --- Benchmark Globals ---
int N_DATA;
int N_GRID;
double BANDWIDTH;

double* data_points;       // Input data points
double* grid_points;       // Points to evaluate the density at
double* density_estimates; // Output density estimates
double final_checksum = 0.0;

// Gaussian kernel constant: 1 / sqrt(2*PI)
const double KERNEL_CONST = 0.3989422804014327;


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_data_points num_grid_points bandwidth seed\n", argv[0]);
        exit(1);
    }

    N_DATA = atoi(argv[1]);
    N_GRID = atoi(argv[2]);
    BANDWIDTH = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if(N_DATA <= 0 || N_GRID <= 0 || BANDWIDTH <= 0.0) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }
    
    mt_seed(seed);

    data_points = (double*)malloc(N_DATA * sizeof(double));
    grid_points = (double*)malloc(N_GRID * sizeof(double));
    density_estimates = (double*)malloc(N_GRID * sizeof(double));

    if (!data_points || !grid_points || !density_estimates) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate random data points from a bimodal distribution
    // Mix of two normal distributions N(0.25, 0.1^2) and N(0.75, 0.1^2)
    for (int i = 0; i < N_DATA; i++) {
        double u1, u2;
        // Box-Muller transform requires u1 in (0, 1] for log()
        do {
            u1 = (double)mt_rand() / UINT32_MAX;
        } while (u1 == 0.0);
        u2 = (double)mt_rand() / UINT32_MAX;
        
        double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);

        // Mix the distributions to create a bimodal sample
        if (i % 2 == 0) {
            data_points[i] = 0.25 + z0 * 0.1; // First mode
        } else {
            data_points[i] = 0.75 + z0 * 0.1; // Second mode
        }
    }
    
    // Create a linear grid of points to evaluate the KDE, covering the data range
    double grid_min = -0.5;
    double grid_max = 1.5;
    double step = (grid_max - grid_min) / (N_GRID - 1);
    for (int i = 0; i < N_GRID; i++) {
        grid_points[i] = grid_min + i * step;
        density_estimates[i] = 0.0;
    }
}

void run_computation() {
    double total_sum = 0.0;
    const double inv_n_h = 1.0 / (N_DATA * BANDWIDTH);
    
    // For each grid point, calculate the kernel density estimate
    for (int i = 0; i < N_GRID; i++) {
        double grid_point = grid_points[i];
        double current_sum = 0.0;
        
        // Sum the kernel function values over all data points
        for (int j = 0; j < N_DATA; j++) {
            double u = (grid_point - data_points[j]) / BANDWIDTH;
            // Gaussian kernel contribution (unscaled)
            current_sum += exp(-0.5 * u * u);
        }
        
        // Final estimate for the grid point
        double estimate = KERNEL_CONST * inv_n_h * current_sum;
        density_estimates[i] = estimate;
        total_sum += estimate; // Accumulate for checksum
    }
    final_checksum = total_sum;
}

void cleanup() {
    free(data_points);
    free(grid_points);
    free(density_estimates);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_checksum);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
