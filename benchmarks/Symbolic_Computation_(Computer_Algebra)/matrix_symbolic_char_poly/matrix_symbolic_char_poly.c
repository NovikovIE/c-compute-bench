/*
 * Program: matrix_symbolic_char_poly
 * Theme:   Symbolic Computation (Computer Algebra)
 * Description: This benchmark computes the characteristic polynomial of a matrix symbolically.
 *              The characteristic polynomial of an N x N matrix A is given by det(A - λ*I),
 *              where I is the identity matrix and λ is a symbolic variable.
 *              The computation is performed using recursive cofactor expansion. The complexity
 *              arises from 'expression swell', where the size and coefficients of intermediate
 *              polynomials grow very rapidly.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>


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

// --- Symbolic Computation Data Structures and Globals ---

// Represents a polynomial P(x) = c_n*x^n + ... + c_1*x^1 + c_0
// coeffs[i] stores the coefficient c_i
typedef struct {
    long long* coeffs;
    int degree;
} Polynomial;

// Global variables for the benchmark
int matrix_dim;
int** A; // Input integer matrix
Polynomial*** M; // Matrix of polynomials (A - λ*I)
Polynomial* final_char_poly; // The result of the computation
long long final_result; // A checksum of the final polynomial's coefficients

// --- Polynomial Helper Functions ---

Polynomial* poly_create(int degree) {
    Polynomial* p = (Polynomial*)malloc(sizeof(Polynomial));
    if (!p) exit(1);
    p->degree = degree;
    // calloc initializes coefficients to zero
    p->coeffs = (long long*)calloc(degree + 1, sizeof(long long));
    if (!p->coeffs) exit(1);
    return p;
}

void poly_free(Polynomial* p) {
    if (p) {
        free(p->coeffs);
        free(p);
    }
}

Polynomial* poly_add(const Polynomial* p1, const Polynomial* p2, int sign) {
    int max_degree = (p1->degree > p2->degree) ? p1->degree : p2->degree;
    Polynomial* result = poly_create(max_degree);

    for (int i = 0; i <= p1->degree; i++) {
        result->coeffs[i] += p1->coeffs[i];
    }
    for (int i = 0; i <= p2->degree; i++) {
        result->coeffs[i] += sign * p2->coeffs[i];
    }
    return result;
}

Polynomial* poly_multiply(const Polynomial* p1, const Polynomial* p2) {
    if (!p1 || !p2) return NULL;
    int new_degree = p1->degree + p2->degree;
    Polynomial* result = poly_create(new_degree);

    for (int i = 0; i <= p1->degree; i++) {
        for (int j = 0; j <= p2->degree; j++) {
            result->coeffs[i + j] += p1->coeffs[i] * p2->coeffs[j];
        }
    }
    return result;
}

// --- Core Symbolic Computation ---

// Create a submatrix by excluding a given row and column.
// The returned matrix structure must be freed, but it contains pointers to the original polynomials.
Polynomial*** create_submatrix(Polynomial*** mat, int n, int row_to_skip, int col_to_skip) {
    Polynomial*** sub_mat = (Polynomial***)malloc((n - 1) * sizeof(Polynomial**));
    if (!sub_mat) exit(1);
    for (int i = 0, r = 0; i < n; i++) {
        if (i == row_to_skip) continue;
        sub_mat[r] = (Polynomial**)malloc((n - 1) * sizeof(Polynomial*));
        if (!sub_mat[r]) exit(1);
        for (int j = 0, c = 0; j < n; j++) {
            if (j == col_to_skip) continue;
            sub_mat[r][c] = mat[i][j];
            c++;
        }
        r++;
    }
    return sub_mat;
}

void free_submatrix_struct(Polynomial*** sub_mat, int n_sub) {
    for (int i = 0; i < n_sub; i++) {
        free(sub_mat[i]);
    }
    free(sub_mat);
}

Polynomial* symbolic_determinant(Polynomial*** mat, int n) {
    if (n == 1) {
        // Base case: for a 1x1 matrix, the determinant is the single element.
        // Return a copy, as all returned polynomials must be owned by the caller.
        Polynomial* result = poly_create(mat[0][0]->degree);
        memcpy(result->coeffs, mat[0][0]->coeffs, (mat[0][0]->degree + 1) * sizeof(long long));
        return result;
    }

    Polynomial* total_det = poly_create(0); // Initialize with P(x) = 0
    int sign = 1; 

    for (int j = 0; j < n; j++) {
        Polynomial*** sub_mat = create_submatrix(mat, n, 0, j);
        Polynomial* sub_det = symbolic_determinant(sub_mat, n - 1);
        
        Polynomial* term = poly_multiply(mat[0][j], sub_det);
        
        Polynomial* old_total = total_det;
        total_det = poly_add(old_total, term, sign);
        
        sign = -sign;

        // Cleanup temporary allocations
        poly_free(old_total);
        poly_free(sub_det);
        poly_free(term);
        free_submatrix_struct(sub_mat, n - 1);
    }
    return total_det;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_dim> <seed>\n", argv[0]);
        exit(1);
    }
    matrix_dim = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (matrix_dim <= 0) {
        fprintf(stderr, "Matrix dimension must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // 1. Allocate and fill the input integer matrix A
    A = (int**)malloc(matrix_dim * sizeof(int*));
    for (int i = 0; i < matrix_dim; i++) {
        A[i] = (int*)malloc(matrix_dim * sizeof(int));
        for (int j = 0; j < matrix_dim; j++) {
            // Generate small integers to prevent coefficient overflow
            A[i][j] = (mt_rand() % 21) - 10;
        }
    }

    // 2. Create the symbolic matrix M = A - λ*I
    M = (Polynomial***)malloc(matrix_dim * sizeof(Polynomial**));
    for (int i = 0; i < matrix_dim; i++) {
        M[i] = (Polynomial**)malloc(matrix_dim * sizeof(Polynomial*));
    }

    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            if (i == j) {
                // Diagonal element: A[i][i] - λ
                M[i][j] = poly_create(1); // Degree 1
                M[i][j]->coeffs[0] = A[i][j];
                M[i][j]->coeffs[1] = -1;
            } else {
                // Off-diagonal element: A[i][j]
                M[i][j] = poly_create(0); // Degree 0 (constant)
                M[i][j]->coeffs[0] = A[i][j];
            }
        }
    }
    final_char_poly = NULL;
}

void run_computation() {
    final_char_poly = symbolic_determinant(M, matrix_dim);
    
    // Calculate a checksum of coefficients to prevent dead code elimination
    final_result = 0;
    if (final_char_poly != NULL) {
        for (int i = 0; i <= final_char_poly->degree; i++) {
            final_result += final_char_poly->coeffs[i];
        }
    }
}

void cleanup() {
    if (final_char_poly) {
        poly_free(final_char_poly);
    }

    // M owns all the initial polynomials
    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            poly_free(M[i][j]);
        }
        free(M[i]);
    }
    free(M);

    // Free the integer matrix
    for (int i = 0; i < matrix_dim; i++) {
        free(A[i]);
    }
    free(A);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout as a checksum
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
