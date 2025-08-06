#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- START MERSENNE TWISTER (MT19937) ---
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
// --- END MERSENNE TWISTER (MT19937) ---

// --- BENCHMARK DATA AND PARAMETERS ---
// Parameters
int num_atoms;
int num_basis_functions_per_atom;
int integration_grid_quality;
int N; // Total number of basis functions: num_atoms * num_basis_functions_per_atom
int num_grid_points; // Total number of integration grid points: num_atoms * integration_grid_quality
const int SCF_ITERATIONS = 5; // Fixed number of SCF iterations for the benchmark

// Data arrays - matrices are stored in row-major order
double* density_matrix;      // P (N x N)
double* fock_matrix;         // F (N x N)
double* h_core_matrix;       // H_core (N x N), one-electron Hamiltonian
double* temp_matrix;         // A temporary matrix for intermediate calculations (N x N)
double* grid_basis_values; // phi(r_g, i) -> values of N basis functions on each grid point (num_grid_points x N)
double* grid_weights;        // w_g -> weights for each grid point (num_grid_points)
double* grid_density;        // rho(r_g) -> electron density on each grid point (num_grid_points)

double final_result = 0.0;

// --- BENCHMARK FUNCTIONS ---

// Helper to generate a random double between 0.0 and 1.0
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_atoms> <num_basis_functions_per_atom> <integration_grid_quality> <seed>\n", argv[0]);
        exit(1);
    }

    num_atoms = atoi(argv[1]);
    num_basis_functions_per_atom = atoi(argv[2]);
    integration_grid_quality = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    N = num_atoms * num_basis_functions_per_atom;
    num_grid_points = num_atoms * integration_grid_quality;

    if (N <= 0 || num_grid_points <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters resulting in zero-sized arrays.\n");
        exit(1);
    }

    // Allocate memory
    density_matrix = (double*)malloc(N * N * sizeof(double));
    fock_matrix = (double*)malloc(N * N * sizeof(double));
    h_core_matrix = (double*)malloc(N * N * sizeof(double));
    temp_matrix = (double*)malloc(N * N * sizeof(double));
    grid_basis_values = (double*)malloc(num_grid_points * N * sizeof(double));
    grid_weights = (double*)malloc(num_grid_points * sizeof(double));
    grid_density = (double*)malloc(num_grid_points * sizeof(double));

    if (!density_matrix || !fock_matrix || !h_core_matrix || !temp_matrix || !grid_basis_values || !grid_weights || !grid_density) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize matrices and vectors with random data
    // Initial guess for density matrix - must be symmetric and trace = num_electrons/2
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
             density_matrix[i * N + j] = (i==j) ? 0.5 : random_double() * 0.01;
        }
    }

    // H_core matrix - symmetric
    for (int i = 0; i < N; ++i) {
        for (int j = i; j < N; ++j) {
            double val = -random_double() * 10.0; // Core hamiltonian elements are typically negative
            h_core_matrix[i * N + j] = val;
            h_core_matrix[j * N + i] = val;
        }
    }

    // Basis function values on the grid
    for (int i = 0; i < num_grid_points * N; ++i) {
        grid_basis_values[i] = random_double() * 2.0 - 1.0; // [-1, 1]
    }

    // Grid weights
    for (int i = 0; i < num_grid_points; ++i) {
        grid_weights[i] = random_double() * 0.01;
    }
}

void run_computation() {
    // Main SCF loop - a fixed number of iterations for benchmarking
    for (int iter = 0; iter < SCF_ITERATIONS; ++iter) {

        // 1. Calculate electron density on the integration grid from the density matrix.
        //    rho(r_g) = sum_{i,j} P_ij * phi_i(r_g) * phi_j(r_g)
        //    This is an O(num_grid_points * N^2) operation.
        for (int g = 0; g < num_grid_points; ++g) {
            double rho_g = 0.0;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    rho_g += density_matrix[i * N + j] * grid_basis_values[g * N + i] * grid_basis_values[g * N + j];
                }
            }
            grid_density[g] = rho_g > 1e-12 ? rho_g : 1e-12; // Avoid zero or negative density
        }

        // 2. Build the Fock matrix: F = H_core + G(P), where G represents the two-electron term.
        //    We model G in two parts:
        //    a) The Exchange-Correlation (XC) part, built via numerical integration.
        //    b) The Coulomb (J) part, which we simulate with a placeholder O(N^3) operation.
        
        // 2a. Calculate the XC contribution and add to H_core to start building the Fock matrix.
        //     V_xc_ij = integral( phi_i(r) * v_xc[rho(r)] * phi_j(r) dr )
        //               ~= sum_g w_g * phi_i(r_g) * v_xc[rho(r_g)] * phi_j(r_g)
        //     This is another O(num_grid_points * N^2) operation.
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double v_xc_ij = 0.0;
                for (int g = 0; g < num_grid_points; ++g) {
                    // Simple Slater exchange functional approximation: v_xc ~ -rho^(1/3)
                    double v_xc_potential = -pow(grid_density[g], 1.0 / 3.0);
                    v_xc_ij += grid_weights[g] * grid_basis_values[g * N + i] * v_xc_potential * grid_basis_values[g * N + j];
                }
                fock_matrix[i * N + j] = h_core_matrix[i * N + j] + v_xc_ij;
            }
        }

        // 2b. Add a placeholder for the computationally expensive Coulomb (J) matrix build.
        //     This is typically O(N^4) but reduced to O(N^3) with modern techniques.
        //     We simulate the O(N^3) cost with a matrix-matrix multiplication.
        for(int i = 0; i < N; ++i) {
            for(int j = 0; j < N; ++j) {
                double g_ij = 0.0;
                for(int k = 0; k < N; ++k) {
                     g_ij += density_matrix[i * N + k] * h_core_matrix[k * N + j] * 0.1; // Using H_core as a stand-in
                }
                fock_matrix[i * N + j] += g_ij; 
            }
        }

        // 3. Update the density matrix. This normally involves diagonalizing the Fock matrix.
        //    Diagonalization is O(N^3). We simulate this cost by performing a matrix multiplication
        //    and forming a new density matrix.
        //    P_new = C * C^T, where C are eigenvectors. We simulate with P_new = F * F^T.
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double p_new_ij = 0.0;
                for (int k = 0; k < N; ++k) {
                    p_new_ij += fock_matrix[i * N + k] * fock_matrix[j * N + k]; // F * F_transpose
                }
                temp_matrix[i * N + j] = p_new_ij;
            }
        }
        
        // Normalize the new density matrix to conserve the number of electrons.
        double trace = 0.0;
        for (int i = 0; i < N; ++i) {
            trace += temp_matrix[i * N + i];
        }
        double norm_factor = (N / 2.0) / (trace + 1e-9); // Target trace is N/2 (num_electron_pairs)

        for(int i = 0; i < N * N; ++i) {
            density_matrix[i] = temp_matrix[i] * norm_factor;
        }
    }

    // Calculate the final result to prevent dead-code elimination.
    // The trace of the density matrix gives the number of electron pairs.
    double total_density_trace = 0.0;
    for (int i = 0; i < N; ++i) {
        total_density_trace += density_matrix[i * N + i];
    }
    final_result = total_density_trace;
}

void cleanup() {
    free(density_matrix);
    free(fock_matrix);
    free(h_core_matrix);
    free(temp_matrix);
    free(grid_basis_values);
    free(grid_weights);
    free(grid_density);
}

int main(int argc, char* argv[]) {
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
