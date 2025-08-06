#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Globals ---

// This benchmark simulates symbolic matrix inversion by tracking expression complexity.
// A symbolic expression's complexity is represented by a 64-bit unsigned integer.
// This avoids overflow, as complexity grows factorially with matrix dimension.
typedef uint64_t SymbolicValue;

// Global parameters and data structures
int matrix_dim;
SymbolicValue **matrix; // Input matrix, initialized with base complexity.
SymbolicValue **inverse_matrix; // Resulting inverse matrix after computation.
SymbolicValue final_result; // Accumulated complexity to prevent dead code elimination.

// --- Symbolic Operation Models ---

// Adding two symbolic expressions results in a new expression whose
// complexity is the sum of the operands' complexities.
SymbolicValue symbolic_add(SymbolicValue a, SymbolicValue b) {
    return a + b;
}

// Multiplying two symbolic expressions (e.g., (a+b)*(c+d)) leads to a combinatorial
// increase in terms. We model this by multiplying their complexities.
// This is the primary driver of the 'expression swell' characteristic of this benchmark.
SymbolicValue symbolic_mul(SymbolicValue a, SymbolicValue b) {
    return a * b;
}

// --- Forward Declarations for Core Logic ---
void get_submatrix(SymbolicValue **mat, SymbolicValue **sub, int p, int q, int n);
SymbolicValue symbolic_determinant(SymbolicValue **mat, int n);

/**
 * @brief Sets up benchmark data.
 * 
 * Parses command-line arguments for matrix dimension and seed.
 * Allocates memory for the input and output matrices.
 * Initializes the input matrix with base symbolic variables, each having a
 * complexity of 1, representing simple variables like 'x_ij'.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_dim> <seed>\n", argv[0]);
        exit(1);
    }

    matrix_dim = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (matrix_dim <= 0) {
        fprintf(stderr, "FATAL: matrix_dim must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    matrix = (SymbolicValue **)malloc(matrix_dim * sizeof(SymbolicValue *));
    inverse_matrix = (SymbolicValue **)malloc(matrix_dim * sizeof(SymbolicValue *));
    for (int i = 0; i < matrix_dim; i++) {
        matrix[i] = (SymbolicValue *)malloc(matrix_dim * sizeof(SymbolicValue));
        inverse_matrix[i] = (SymbolicValue *)malloc(matrix_dim * sizeof(SymbolicValue));
    }

    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            // Each element of the initial matrix is a simple variable, so its complexity is 1.
            matrix[i][j] = 1;
            mt_rand(); // Consume a random number.
        }
    }
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    for (int i = 0; i < matrix_dim; i++) {
        free(matrix[i]);
        free(inverse_matrix[i]);
    }
    free(matrix);
    free(inverse_matrix);
}

/**
 * @brief Helper to create a submatrix by excluding a specific row and column.
 * @param mat Source matrix of size n x n.
 * @param sub Destination submatrix of size (n-1) x (n-1).
 * @param p Row to exclude.
 * @param q Column to exclude.
 * @param n Size of the source matrix.
 */
void get_submatrix(SymbolicValue **mat, SymbolicValue **sub, int p, int q, int n) {
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        if (row == p) continue;
        j = 0;
        for (int col = 0; col < n; col++) {
            if (col == q) continue;
            sub[i][j++] = mat[row][col];
        }
        i++;
    }
}

/**
 * @brief Recursively calculates the complexity of a symbolic determinant.
 * 
 * This function's factorial time complexity is the heart of the benchmark.
 * It allocates temporary storage for submatrices at each step of the recursion,
 * simulating the creation of intermediate expressions in a real computer algebra system.
 * @param mat The matrix for which to calculate the determinant complexity.
 * @param n The dimension of the matrix.
 * @return The calculated complexity of the determinant.
 */
SymbolicValue symbolic_determinant(SymbolicValue **mat, int n) {
    if (n == 1) {
        return mat[0][0];
    }

    SymbolicValue det = 0;

    // Allocate temporary submatrix for cofactor expansion.
    SymbolicValue **submat = (SymbolicValue **)malloc((n - 1) * sizeof(SymbolicValue *));
    for (int i = 0; i < n - 1; i++) {
        submat[i] = (SymbolicValue *)malloc((n - 1) * sizeof(SymbolicValue));
    }

    // Expand along the first row to calculate determinant.
    for (int j = 0; j < n; j++) {
        get_submatrix(mat, submat, 0, j, n);
        SymbolicValue sub_det = symbolic_determinant(submat, n - 1);
        
        // The complexity of a term in the expansion is the product of the element's
        // complexity and its cofactor's (sub-determinant) complexity.
        SymbolicValue term_complexity = symbolic_mul(mat[0][j], sub_det);
        det = symbolic_add(det, term_complexity);
    }

    for (int i = 0; i < n - 1; i++) {
        free(submat[i]);
    }
    free(submat);

    return det;
}

/**
 * @brief Runs the core symbolic computation.
 * 
 * This function simulates finding the inverse of a symbolic matrix A, which is
 * given by `A^-1 = adj(A) / det(A)`.
 * 1. It calculates the complexity of each element of the adjugate matrix, `adj(A)`.
 * 2. It calculates the complexity of the determinant, `det(A)`.
 * 3. It models the complexity of the final inverse matrix elements and accumulates
 *    them into a single result.
 */
void run_computation() {
    if (matrix_dim < 1) {
        final_result = 0;
        return;
    }
    if (matrix_dim == 1) {
         // adj[0][0]=1, det=matrix[0][0]=1. inv[0][0] = 1+1=2.
        final_result = symbolic_add(1, matrix[0][0]);
        inverse_matrix[0][0] = final_result;
        return;
    }

    // Allocate temporary storage for adjugate matrix and its subproblems.
    SymbolicValue **adj = (SymbolicValue **)malloc(matrix_dim * sizeof(SymbolicValue *));
    SymbolicValue **temp_submatrix = (SymbolicValue **)malloc((matrix_dim - 1) * sizeof(SymbolicValue *));
    for (int i = 0; i < matrix_dim; i++) {
        adj[i] = (SymbolicValue *)malloc(matrix_dim * sizeof(SymbolicValue));
    }
    for (int i = 0; i < matrix_dim - 1; i++) {
        temp_submatrix[i] = (SymbolicValue *)malloc((matrix_dim - 1) * sizeof(SymbolicValue));
    }

    // 1. Calculate adjugate matrix complexity.
    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            get_submatrix(matrix, temp_submatrix, i, j, matrix_dim);
            adj[j][i] = symbolic_determinant(temp_submatrix, matrix_dim - 1);
        }
    }

    // 2. Calculate determinant complexity.
    SymbolicValue det = symbolic_determinant(matrix, matrix_dim);

    // 3. Calculate final inverse matrix complexity and accumulate result.
    // Complexity(A/B) is modeled as Complexity(A) + Complexity(B).
    final_result = 0;
    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            inverse_matrix[i][j] = symbolic_add(adj[i][j], det);
            final_result = symbolic_add(final_result, inverse_matrix[i][j]);
        }
    }

    // Free temporary matrices allocated within this function.
    for (int i = 0; i < matrix_dim; i++) free(adj[i]);
    free(adj);
    for (int i = 0; i < matrix_dim - 1; i++) free(temp_submatrix[i]);
    free(temp_submatrix);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated complexity to stdout
    printf("%llu\n", (unsigned long long)final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
