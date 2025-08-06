#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// Benchmark parameters
int num_atoms;
int num_basis_functions;

// Derived parameters
int hessian_dim; // 3 * num_atoms

// Data structures
double *atomic_masses;      // Size: num_atoms
double *hessian;            // Size: hessian_dim * hessian_dim
double *eigenvectors;       // Size: hessian_dim * hessian_dim. Stores eigenvector matrix.

// Final result accumulator
double result_sum;

// Generates a random double between -1.0 and 1.0.
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_atoms> <num_basis_functions> <seed>\n", argv[0]);
        exit(1);
    }

    num_atoms = atoi(argv[1]);
    num_basis_functions = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    hessian_dim = 3 * num_atoms;
    mt_seed(seed);

    // Allocate memory
    atomic_masses = (double *)malloc(num_atoms * sizeof(double));
    hessian = (double *)malloc(hessian_dim * hessian_dim * sizeof(double));
    eigenvectors = (double *)malloc(hessian_dim * hessian_dim * sizeof(double));

    if (!atomic_masses || !hessian || !eigenvectors) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize atomic masses (e.g., between 1.0 for H and 200 for heavier elements)
    for (int i = 0; i < num_atoms; ++i) {
        atomic_masses[i] = 1.0 + (mt_rand() / (double)UINT32_MAX) * 199.0;
    }

    // Initialize the Hessian matrix with random symmetric values.
    // The inner loop simulates the expensive part of the calculation, scaling with basis functions.
    long long fake_work_load = (long long)num_basis_functions * num_basis_functions;
    for (int i = 0; i < hessian_dim; ++i) {
        for (int j = i; j < hessian_dim; ++j) {
            double temp_sum = 0.0;
            // Simulate O(B^2) work for calculating each Hessian element.
            for(int k=0; k < num_basis_functions; ++k) {
                 temp_sum += rand_double() * rand_double(); 
            }
            double val = temp_sum / num_basis_functions;
            hessian[i * hessian_dim + j] = val; 
            hessian[j * hessian_dim + i] = val; // Ensure symmetry
        }
    }
    
    // Initialize eigenvectors matrix to identity
    for(int i = 0; i < hessian_dim; ++i) {
        for(int j = 0; j < hessian_dim; ++j) {
            eigenvectors[i * hessian_dim + j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

void run_computation() {
    // --- Step 1: Mass-weight the Hessian matrix ---
    // H_mw[i,j] = H[i,j] / sqrt(mass(i) * mass(j))
    for (int i = 0; i < hessian_dim; ++i) {
        for (int j = 0; j < hessian_dim; ++j) {
            hessian[i * hessian_dim + j] /= sqrt(atomic_masses[i / 3] * atomic_masses[j / 3]);
        }
    }

    // --- Step 2: Diagonalize the matrix using Jacobi method ---
    // For a benchmark, a fixed number of sweeps is sufficient.
    const int num_sweeps = 8;
    for (int sweep = 0; sweep < num_sweeps; ++sweep) {
        for (int p = 0; p < hessian_dim - 1; ++p) {
            for (int q = p + 1; q < hessian_dim; ++q) {
                double app = hessian[p * hessian_dim + p];
                double aqq = hessian[q * hessian_dim + q];
                double apq = hessian[p * hessian_dim + q];

                if (fabs(apq) < 1e-12) continue;

                double tau = (aqq - app) / (2.0 * apq);
                double t = (tau >= 0) ? (1.0 / (tau + sqrt(1.0 + tau * tau))) : (-1.0 / (-tau + sqrt(1.0 + tau * tau)));
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;

                // Update Hessian matrix A' = R^T * A * R
                for (int i = 0; i < hessian_dim; ++i) {
                    double aip = hessian[i * hessian_dim + p];
                    double aiq = hessian[i * hessian_dim + q];
                    hessian[i * hessian_dim + p] = c * aip - s * aiq;
                    hessian[p * hessian_dim + i] = c * aip - s * aiq;
                    hessian[i * hessian_dim + q] = s * aip + c * aiq;
                    hessian[q * hessian_dim + i] = s * aip + c * aiq;
                }

                hessian[p * hessian_dim + p] = app - t * apq;
                hessian[q * hessian_dim + q] = aqq + t * apq;
                hessian[p * hessian_dim + q] = 0.0;
                hessian[q * hessian_dim + p] = 0.0;
            }
        }
    }

    // --- Step 3: Calculate frequencies from eigenvalues ---
    // Eigenvalues are now on the diagonal of the modified Hessian.
    // Frequency is proportional to sqrt(eigenvalue).
    // Sum them up to produce a single result.
    result_sum = 0.0;
    for (int i = 0; i < hessian_dim; ++i) {
        double eigenvalue = hessian[i * hessian_dim + i];
        // Physical frequencies require positive eigenvalues, use fabs for robustness.
        result_sum += sqrt(fabs(eigenvalue));
    }
}

void cleanup() {
    free(atomic_masses);
    free(hessian);
    free(eigenvectors);
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
    printf("%.6f\n", result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}