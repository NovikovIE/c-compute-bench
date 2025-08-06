#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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

// --- BENCHMARK DATA AND PARAMETERS ---
int num_samples;
int num_features;
int num_components;

// Input data matrix
double** data;

// Intermediate and final matrices/vectors
double** temp_cov_matrix; // Writable copy of covariance matrix for deflation
double** eigenvectors;
double* eigenvalues;

// Final result accumulator
double final_result = 0.0;

// --- UTILITY FUNCTIONS ---
double generate_random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

double** allocate_2d_matrix(int rows, int cols) {
    double** matrix = (double**)malloc(rows * sizeof(double*));
    if (!matrix) return NULL;
    for (int i = 0; i < rows; i++) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        if (!matrix[i]) {
            // Clean up previously allocated rows
            for (int j = 0; j < i; j++) free(matrix[j]);
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

void free_2d_matrix(double** matrix, int rows) {
    if (matrix) {
        for (int i = 0; i < rows; i++) {
            free(matrix[i]);
        }
        free(matrix);
    }
}

// --- BENCHMARK IMPLEMENTATION ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_samples num_features num_components seed\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_features = atoi(argv[2]);
    num_components = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (num_components > num_features) {
        fprintf(stderr, "Error: num_components cannot be greater than num_features.\n");
        exit(1);
    }

    data = allocate_2d_matrix(num_samples, num_features);
    temp_cov_matrix = allocate_2d_matrix(num_features, num_features);
    eigenvectors = allocate_2d_matrix(num_components, num_features);
    eigenvalues = (double*)malloc(num_components * sizeof(double));

    if (!data || !temp_cov_matrix || !eigenvectors || !eigenvalues) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_features; j++) {
            data[i][j] = generate_random_double() * 100.0;
        }
    }
}

void run_computation() {
    // Step 1: Center the data (subtract the mean of each feature)
    double* mean = (double*)calloc(num_features, sizeof(double));
    for (int j = 0; j < num_features; j++) {
        for (int i = 0; i < num_samples; i++) {
            mean[j] += data[i][j];
        }
        mean[j] /= num_samples;
    }

    for (int i = 0; i < num_samples; i++) {
        for (int j = 0; j < num_features; j++) {
            data[i][j] -= mean[j];
        }
    }
    free(mean);

    // Step 2: Compute the covariance matrix
    // Cov(X, Y) = E[(X - E[X])(Y - E[Y])] = (X_centered' * X_centered) / (n-1)
    for (int i = 0; i < num_features; i++) {
        for (int j = i; j < num_features; j++) {
            double cov = 0.0;
            for (int k = 0; k < num_samples; k++) {
                cov += data[k][i] * data[k][j];
            }
            temp_cov_matrix[i][j] = cov / (num_samples - 1);
            temp_cov_matrix[j][i] = temp_cov_matrix[i][j]; // Symmetric matrix
        }
    }

    // Step 3: Find top K eigenvectors/eigenvalues using Power Iteration with Deflation
    double* vec_b = (double*)malloc(num_features * sizeof(double));
    double* vec_Ab = (double*)malloc(num_features * sizeof(double));
    const int POWER_ITERATIONS = 100;

    for (int k = 0; k < num_components; k++) {
        // Initialize a random vector b
        for (int i = 0; i < num_features; i++) {
            vec_b[i] = generate_random_double();
        }

        // Power Iteration to find the dominant eigenvector
        for (int iter = 0; iter < POWER_ITERATIONS; iter++) {
            // Matrix-vector multiplication: Ab = A * b
            for (int i = 0; i < num_features; i++) {
                vec_Ab[i] = 0.0;
                for (int j = 0; j < num_features; j++) {
                    vec_Ab[i] += temp_cov_matrix[i][j] * vec_b[j];
                }
            }
            // Normalize the vector
            double norm = 0.0;
            for (int i = 0; i < num_features; i++) {
                norm += vec_Ab[i] * vec_Ab[i];
            }
            norm = sqrt(norm);
            for (int i = 0; i < num_features; i++) {
                vec_b[i] = vec_Ab[i] / norm;
            }
        }

        // Store the eigenvector
        for (int i = 0; i < num_features; i++) {
            eigenvectors[k][i] = vec_b[i];
        }

        // Calculate eigenvalue: lambda = b' * A * b
        for (int i = 0; i < num_features; i++) {
            vec_Ab[i] = 0.0;
            for (int j = 0; j < num_features; j++) {
                vec_Ab[i] += temp_cov_matrix[i][j] * vec_b[j];
            }
        }
        double lambda = 0.0;
        for (int i = 0; i < num_features; i++) {
            lambda += vec_b[i] * vec_Ab[i];
        }
        eigenvalues[k] = lambda;

        // Deflate the matrix: A = A - lambda * v * v'
        for (int i = 0; i < num_features; i++) {
            for (int j = 0; j < num_features; j++) {
                temp_cov_matrix[i][j] -= lambda * vec_b[i] * vec_b[j];
            }
        }
    }

    free(vec_b);
    free(vec_Ab);

    // Step 4: Accumulate result to prevent dead code elimination
    final_result = 0.0;
    for (int i = 0; i < num_components; i++) {
        final_result += eigenvalues[i];
    }
}

void cleanup() {
    free_2d_matrix(data, num_samples);
    free_2d_matrix(temp_cov_matrix, num_features);
    free_2d_matrix(eigenvectors, num_components);
    free(eigenvalues);
}

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
