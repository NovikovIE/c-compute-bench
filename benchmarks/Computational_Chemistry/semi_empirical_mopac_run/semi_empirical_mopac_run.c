#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// Benchmark parameters and data structures
static int NUM_ATOMS;
static double *coords;        // 3D coordinates for each atom (size: NUM_ATOMS * 3)
static double *h_core;          // Core Hamiltonian matrix (size: NUM_ATOMS * NUM_ATOMS)
static double *density_matrix;  // Density matrix (size: NUM_ATOMS * NUM_ATOMS)
static double *fock_matrix;     // Fock matrix (size: NUM_ATOMS * NUM_ATOMS)

static double total_energy;    // Final result

#define SCF_ITERATIONS 10

// Helper to generate a random double in [min, max]
double rand_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

// Helper to calculate squared distance between two atoms
double dist_sq(int i, int j) {
    double dx = coords[i * 3 + 0] - coords[j * 3 + 0];
    double dy = coords[i * 3 + 1] - coords[j * 3 + 1];
    double dz = coords[i * 3 + 2] - coords[j * 3 + 2];
    return dx * dx + dy * dy + dz * dz;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_atoms> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_ATOMS = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (NUM_ATOMS <= 0) {
        fprintf(stderr, "FATAL: num_atoms must be positive.\n");
        exit(1);
    }

    size_t n = (size_t)NUM_ATOMS;
    coords = (double*)malloc(n * 3 * sizeof(double));
    h_core = (double*)malloc(n * n * sizeof(double));
    density_matrix = (double*)malloc(n * n * sizeof(double));
    fock_matrix = (double*)malloc(n * n * sizeof(double));

    if (!coords || !h_core || !density_matrix || !fock_matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize atomic coordinates in a box of side 20 Angstroms
    for (size_t i = 0; i < n * 3; ++i) {
        coords[i] = rand_double(-10.0, 10.0);
    }

    // Initialize matrices with random data
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // Core Hamiltonian elements are typically negative and large
            h_core[i * n + j] = rand_double(-50.0, -1.0);
            // Density matrix elements are typically smaller
            density_matrix[i * n + j] = rand_double(0.0, 2.0 / n);
        }
    }

    total_energy = 0.0;
}

void run_computation() {
    size_t n = (size_t)NUM_ATOMS;

    // Simulate a Self-Consistent Field (SCF) procedure for a fixed number of iterations.
    // This is the core computational loop in many quantum chemistry methods.
    for (int iter = 0; iter < SCF_ITERATIONS; ++iter) {
        // 1. Build the Fock matrix from the density matrix.
        // This is the most computationally expensive step, scaling as O(N^3) or O(N^4)
        // in real methods. We simulate a O(N^3) kernel.
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) { // Exploit F_ij = F_ji symmetry
                double two_electron_term = 0.0;
                // This O(N) inner loop makes the total complexity O(N^3).
                // It simulates the evaluation of two-electron repulsion integrals.
                for (size_t k = 0; k < n; ++k) {
                     double r_ik_sq = dist_sq(i, k);
                     double r_jk_sq = dist_sq(j, k);
                     // A fictitious term that depends on i,j,k to ensure O(N^3) work
                     if (r_ik_sq < 1e-9 || r_jk_sq < 1e-9) continue;
                     double interaction = (density_matrix[k * n + k]) / (sqrt(r_ik_sq) + sqrt(r_jk_sq));
                     two_electron_term += interaction;
                }
                // Fock = Core Hamiltonian + Two-Electron Term
                double val = h_core[i * n + j] + two_electron_term;
                fock_matrix[i * n + j] = val;
                fock_matrix[j * n + i] = val; // Set symmetric element
            }
        }

        // 2. Calculate the electronic energy for this iteration.
        // Energy = 0.5 * Sum(P_ij * (H_core_ij + F_ij)) (O(N^2) work)
        double current_energy = 0.0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                current_energy += density_matrix[i * n + j] * (h_core[i * n + j] + fock_matrix[i * n + j]);
            }
        }
        total_energy += 0.5 * current_energy;

        // 3. Update the density matrix for the next iteration.
        // In a real calculation, this involves diagonalizing the Fock matrix (O(N^3)).
        // Here, we simulate the update with a simple mixing scheme (damping) to create O(N^2) work.
        for (size_t i = 0; i < n * n; ++i) {
            density_matrix[i] = 0.8 * density_matrix[i] + 0.2 * (fock_matrix[i] / (n * 100.0));
        }
    }
}

void cleanup() {
    free(coords);
    free(h_core);
    free(density_matrix);
    free(fock_matrix);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated energy to stdout to prevent dead-code elimination
    printf("%.4f\n", total_energy);

    // Print the timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
