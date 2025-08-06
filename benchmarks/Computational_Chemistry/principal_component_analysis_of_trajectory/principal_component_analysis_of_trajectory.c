#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) BEGIN ---
// Do not modify this section.
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
                fprintf(stderr, "FATAL: Mersenne Twister not seeded.\n");
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
// --- Mersenne Twister (MT19937) END ---

// --- Benchmark Globals ---
int num_frames;
int num_atoms;
int num_dof; // Degrees of Freedom (3 * num_atoms)

// Flattened 2D array: [frame * num_dof + dof_idx]
double *trajectory;

// 1D array: [dof_idx]
double *avg_coords;

// Flattened 2D array: [row_dof * num_dof + col_dof]
double *covariance_matrix;

double final_result = 0.0;

// --- Utility Function ---
double rand_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_frames num_atoms seed\n", argv[0]);
        exit(1);
    }

    num_frames = atoi(argv[1]);
    num_atoms = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_frames <= 0 || num_atoms <= 0) {
        fprintf(stderr, "ERROR: num_frames and num_atoms must be positive integers.\n");
        exit(1);
    }

    num_dof = num_atoms * 3; // Each atom has x, y, z coordinates

    mt_seed(seed);

    // Allocate memory
    trajectory = (double*)malloc((size_t)num_frames * num_dof * sizeof(double));
    avg_coords = (double*)malloc((size_t)num_dof * sizeof(double));
    covariance_matrix = (double*)malloc((size_t)num_dof * num_dof * sizeof(double));

    if (!trajectory || !avg_coords || !covariance_matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize trajectory with random coordinates (e.g., in Angstroms)
    for (int i = 0; i < num_frames; ++i) {
        for (int j = 0; j < num_dof; ++j) {
            trajectory[i * num_dof + j] = rand_double(-10.0, 10.0);
        }
    }
}

void run_computation() {
    // 1. Calculate the average coordinates over all frames
    for (int j = 0; j < num_dof; ++j) {
        double sum = 0.0;
        for (int i = 0; i < num_frames; ++i) {
            sum += trajectory[i * num_dof + j];
        }
        avg_coords[j] = sum / (double)num_frames;
    }

    // 2. Build the covariance matrix.
    // The element C(i,j) is the covariance between degree of freedom i and j.
    // C(i,j) = <(pos_i - <pos_i>) * (pos_j - <pos_j>)>_frames
    for (int i = 0; i < num_dof; ++i) {
        // Exploit symmetry: C(i,j) == C(j,i), so only compute the upper-triangle
        for (int j = i; j < num_dof; ++j) {
            double cov_sum = 0.0;
            for (int k = 0; k < num_frames; ++k) {
                double dev_i = trajectory[k * num_dof + i] - avg_coords[i];
                double dev_j = trajectory[k * num_dof + j] - avg_coords[j];
                cov_sum += dev_i * dev_j;
            }
            double cov_val = cov_sum / (double)num_frames;
            covariance_matrix[i * num_dof + j] = cov_val;
            if (i != j) {
                covariance_matrix[j * num_dof + i] = cov_val;
            }
        }
    }

    // 3. To prevent dead-code elimination, calculate a final result.
    // The trace of the covariance matrix is a physically meaningful quantity
    // representing the total mean square fluctuation of the system.
    double trace = 0.0;
    for (int i = 0; i < num_dof; ++i) {
        trace += covariance_matrix[i * num_dof + i];
    }
    final_result = trace;
}

void cleanup() {
    free(trajectory);
    free(avg_coords);
    free(covariance_matrix);
}

// --- Main and Timing ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%.6f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}