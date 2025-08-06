#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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

// --- Benchmark Globals ---
int N; // num_basis_functions
int MAX_ITER; // max_scf_iterations

// Matrices stored in 1D row-major layout
double* H_core; // Core Hamiltonian (N*N)
double* P;      // Density matrix (N*N)
double* F;      // Fock matrix (N*N)
double* ERI;    // Two-electron repulsion integrals (N*N*N*N)

double final_energy; // Result of the computation

// Helper to generate a random double in [0, 1)
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_basis_functions> <max_scf_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    MAX_ITER = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (N <= 0 || MAX_ITER <= 0) {
        fprintf(stderr, "Parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    size_t n2 = (size_t)N * N;
    size_t n4 = n2 * n2;

    H_core = (double*)malloc(n2 * sizeof(double));
    P      = (double*)malloc(n2 * sizeof(double));
    F      = (double*)malloc(n2 * sizeof(double));
    ERI    = (double*)malloc(n4 * sizeof(double));

    if (!H_core || !P || !F || !ERI) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        free(H_core); free(P); free(F); free(ERI);
        exit(1);
    }

    // Initialize matrices with random data
    // This simulates pre-calculating integrals and an initial guess for the density matrix
    for (size_t i = 0; i < n2; ++i) {
        H_core[i] = random_double() * 2.0 - 1.0; // [-1, 1)
        P[i] = random_double(); // [0, 1)
    }

    for (size_t i = 0; i < n4; ++i) {
        ERI[i] = random_double() * 0.1; // Integrals are often small
    }

    final_energy = 0.0;
}

void run_computation() {
    double current_energy = 0.0;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        // 1. Build the Fock matrix: F = H_core + G
        // G_ij = sum_{k,l} P_lk * (2*(ij|kl) - (ik|jl))
        // This is the dominant O(N^4) step of the calculation.
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double g_ij = 0.0;
                for (int k = 0; k < N; ++k) {
                    for (int l = 0; l < N; ++l) {
                        size_t p_idx_kl = (size_t)k * N + l;
                        size_t eri_idx_ijkl = ((((size_t)i * N + j) * N + k) * N + l);
                        size_t eri_idx_ikjl = ((((size_t)i * N + k) * N + j) * N + l);
                        g_ij += P[p_idx_kl] * (2.0 * ERI[eri_idx_ijkl] - ERI[eri_idx_ikjl]);
                    }
                }
                size_t f_idx_ij = (size_t)i * N + j;
                F[f_idx_ij] = H_core[f_idx_ij] + g_ij;
            }
        }

        // 2. Compute the new energy
        // E = 0.5 * sum_{i,j} P_ij * (H_core_ij + F_ij)
        current_energy = 0.0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                size_t idx = (size_t)i * N + j;
                current_energy += 0.5 * P[idx] * (H_core[idx] + F[idx]);
            }
        }

        // 3. Update the density matrix for the next iteration.
        // A real SCF would use the eigenvectors of F, but we simulate this
        // by mixing the old P with a function of F (damping).
        // This prevents the calculation from being static.
        for (size_t i = 0; i < (size_t)N * N; ++i) {
             P[i] = 0.5 * P[i] + 0.5 * (F[i] / (double)N); // Simple mixing
        }
    }

    final_energy = current_energy;
}

void cleanup() {
    free(H_core);
    free(P);
    free(F);
    free(ERI);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", final_energy);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
