#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
int NUM_ATOMS;
int NUM_BASIS_FUNCTIONS;
double* density_matrix;      // Size: N_bf x N_bf
double* overlap_integrals;   // Size: N_bf x N_bf
double* shielding_tensors;   // Size: N_at x 3 x 3
double final_result;

// --- Utility Functions ---
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_atoms> <num_basis_functions> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_ATOMS = atoi(argv[1]);
    NUM_BASIS_FUNCTIONS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    if (NUM_ATOMS <= 0 || NUM_BASIS_FUNCTIONS <= 0) {
         fprintf(stderr, "FATAL: Invalid parameters. Must be positive integers.\n");
         exit(1);
    }

    long N = NUM_BASIS_FUNCTIONS;
    density_matrix = (double*)malloc(N * N * sizeof(double));
    overlap_integrals = (double*)malloc(N * N * sizeof(double));
    shielding_tensors = (double*)malloc(NUM_ATOMS * 3 * 3 * sizeof(double));

    if (!density_matrix || !overlap_integrals || !shielding_tensors) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (long i = 0; i < N * N; ++i) {
        density_matrix[i] = random_double() * 2.0 - 1.0; // Random values in [-1, 1]
        overlap_integrals[i] = random_double();
    }

    for (int i = 0; i < NUM_ATOMS * 9; ++i) {
        shielding_tensors[i] = 0.0;
    }
}

void run_computation() {
    double checksum = 0.0;
    long N = NUM_BASIS_FUNCTIONS;

    // This loop simulates the expensive part of calculating NMR shielding tensors,
    // which involves a contraction of a density matrix with four-index two-electron integrals.
    // The complexity scales as O(num_atoms * num_basis_functions^4).
    for (int atom = 0; atom < NUM_ATOMS; ++atom) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double tensor_element = 0.0;
                // A mock perturbation factor, constant for this tensor element
                double perturbation = (double)(atom + i + j + 1) * 1e-5;

                for (long mu = 0; mu < N; ++mu) {
                    for (long nu = 0; nu < N; ++nu) {
                        double D_munu = density_matrix[mu * N + nu];
                        double term1 = D_munu * perturbation;
                        for (long lambda = 0; lambda < N; ++lambda) {
                             for (long sigma = 0; sigma < N; ++sigma) {
                                // Simulate calculation of (mu nu|lambda sigma) integral contribution
                                double V = (overlap_integrals[mu * N + lambda] * overlap_integrals[nu * N + sigma])
                                       - 0.5 * (overlap_integrals[mu * N + sigma] * overlap_integrals[nu * N + lambda]);
                                tensor_element += term1 * V;
                            }
                        }
                    }
                }
                shielding_tensors[atom * 9 + i * 3 + j] = tensor_element;
            }
        }
    }

    // Accumulate the final result to prevent dead-code elimination by the compiler
    for (int i = 0; i < NUM_ATOMS * 9; ++i) {
        checksum += shielding_tensors[i];
    }
    final_result = checksum;
}

void cleanup() {
    free(density_matrix);
    free(overlap_integrals);
    free(shielding_tensors);
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
