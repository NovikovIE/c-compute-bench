#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
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

// Benchmark parameters
int g_num_atoms;
int g_max_iterations;
double g_alpha;
double g_epsilon;
double g_sigma;

// Global data structures
double* g_positions; // Size: num_atoms * 3 (x, y, z)
double* g_forces;    // Size: num_atoms * 3 (fx, fy, fz)

// Final result
double g_final_energy;

// Function to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_atoms max_iterations alpha epsilon sigma seed\n", argv[0]);
        exit(1);
    }

    g_num_atoms = atoi(argv[1]);
    g_max_iterations = atoi(argv[2]);
    g_alpha = atof(argv[3]);
    g_epsilon = atof(argv[4]);
    g_sigma = atof(argv[5]);
    uint32_t seed = atoi(argv[6]);

    if (g_num_atoms <= 1 || g_max_iterations <= 0) {
        fprintf(stderr, "Error: num_atoms must be > 1 and max_iterations must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    g_positions = (double*)malloc(g_num_atoms * 3 * sizeof(double));
    g_forces = (double*)malloc(g_num_atoms * 3 * sizeof(double));

    if (!g_positions || !g_forces) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize atom positions randomly in a box scaled by atom count
    double box_size_L = cbrt(g_num_atoms) * g_sigma * 1.5;
    for (int i = 0; i < g_num_atoms; ++i) {
        g_positions[i * 3 + 0] = rand_double() * box_size_L; // x
        g_positions[i * 3 + 1] = rand_double() * box_size_L; // y
        g_positions[i * 3 + 2] = rand_double() * box_size_L; // z
    }
}

void run_computation() {
    double sigma_sq = g_sigma * g_sigma;
    double const_24_epsilon = 24.0 * g_epsilon;

    for (int iter = 0; iter < g_max_iterations; ++iter) {
        // Zero out forces from previous iteration
        for (int i = 0; i < g_num_atoms * 3; ++i) {
            g_forces[i] = 0.0;
        }

        // Calculate pairwise forces (Lennard-Jones)
        for (int i = 0; i < g_num_atoms; ++i) {
            for (int j = i + 1; j < g_num_atoms; ++j) {
                double dx = g_positions[i * 3 + 0] - g_positions[j * 3 + 0];
                double dy = g_positions[i * 3 + 1] - g_positions[j * 3 + 1];
                double dz = g_positions[i * 3 + 2] - g_positions[j * 3 + 2];

                double r_sq = dx * dx + dy * dy + dz * dz;

                // Avoid division by zero with a small cutoff
                if (r_sq < 0.01) r_sq = 0.01;

                double s_over_r_sq = sigma_sq / r_sq;
                double s_over_r_6 = s_over_r_sq * s_over_r_sq * s_over_r_sq;
                double s_over_r_12 = s_over_r_6 * s_over_r_6;

                // F/r = (24 * epsilon / r^2) * [2 * (sigma/r)^12 - (sigma/r)^6]
                double force_over_r = const_24_epsilon / r_sq * (2.0 * s_over_r_12 - s_over_r_6);

                double fx = force_over_r * dx;
                double fy = force_over_r * dy;
                double fz = force_over_r * dz;

                // Add forces according to Newton's 3rd law
                g_forces[i * 3 + 0] += fx;
                g_forces[i * 3 + 1] += fy;
                g_forces[i * 3 + 2] += fz;

                g_forces[j * 3 + 0] -= fx;
                g_forces[j * 3 + 1] -= fy;
                g_forces[j * 3 + 2] -= fz;
            }
        }

        // Update positions using steepest descent method
        for (int i = 0; i < g_num_atoms; ++i) {
            g_positions[i * 3 + 0] += g_alpha * g_forces[i * 3 + 0];
            g_positions[i * 3 + 1] += g_alpha * g_forces[i * 3 + 1];
            g_positions[i * 3 + 2] += g_alpha * g_forces[i * 3 + 2];
        }
    }

    // Calculate final total potential energy as the result to prevent dead code elimination
    g_final_energy = 0.0;
    double const_4_epsilon = 4.0 * g_epsilon;
    for (int i = 0; i < g_num_atoms; ++i) {
        for (int j = i + 1; j < g_num_atoms; ++j) {
            double dx = g_positions[i * 3 + 0] - g_positions[j * 3 + 0];
            double dy = g_positions[i * 3 + 1] - g_positions[j * 3 + 1];
            double dz = g_positions[i * 3 + 2] - g_positions[j * 3 + 2];

            double r_sq = dx * dx + dy * dy + dz * dz;

            if (r_sq < 0.01) r_sq = 0.01;
            
            double s_over_r_sq = sigma_sq / r_sq;
            double s_over_r_6 = s_over_r_sq * s_over_r_sq * s_over_r_sq;
            double s_over_r_12 = s_over_r_6 * s_over_r_6;

            g_final_energy += const_4_epsilon * (s_over_r_12 - s_over_r_6);
        }
    }
}

void cleanup() {
    free(g_positions);
    free(g_forces);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final energy to stdout
    printf("%f\n", g_final_energy);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
