#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

/* --- MERSENNE TWISTER (DO NOT MODIFY) --- */
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
/* --- END MERSENNE TWISTER --- */


/* --- BENCHMARK GLOBALS --- */

// Benchmark Parameters
static int NUM_VARIABLES;
static int POLY_DEGREE1;
static int POLY_DEGREE2;
static int COEFF_POLY_TERMS; // Derived from NUM_VARIABLES
static unsigned int SEED;

// Polynomial Data
// P(x1) = sum_{i=0..d1} C_i(x2..xv) * x1^i
// Q(x1) = sum_{j=0..d2} D_j(x2..xv) * x1^j
// We represent C_i and D_j (coefficient polynomials) as dense arrays of doubles.
static double** p_coeffs; // array of d1+1 pointers to coefficient polynomials
static double** q_coeffs; // array of d2+1 pointers to coefficient polynomials
static double* zero_poly; // A pre-allocated zero polynomial for padding

// Workspace for determinant calculation
#define TEMP_POLY_POOL_SIZE 256 // Max recursion depth * 2
static double* temp_poly_pool[TEMP_POLY_POOL_SIZE];
static int temp_poly_pool_idx;
static double*** sub_matrix_scratch; // Scratch space for sub-matrices

// Final result accumulator
static double final_result_sum;


/* --- FORWARD DECLARATIONS --- */
double* determinant(double*** M, int n);

/* --- POLYNOMIAL ARITHMETIC --- */

// Result <- A + B
void poly_add(const double* a, const double* b, double* result) {
    for (int i = 0; i < COEFF_POLY_TERMS; i++) {
        result[i] = a[i] + b[i];
    }
}

// Result <- A - B
void poly_sub(const double* a, const double* b, double* result) {
    for (int i = 0; i < COEFF_POLY_TERMS; i++) {
        result[i] = a[i] - b[i];
    }
}

// Result <- A * B (convolution with wrap-around)
// This simulates multiplication where expression swell is handled by combining
// terms in a system with a fixed number of available term representations.
void poly_mult(const double* a, const double* b, double* result) {
    memset(result, 0, sizeof(double) * COEFF_POLY_TERMS);
    for (int i = 0; i < COEFF_POLY_TERMS; i++) {
        if (a[i] == 0.0) continue;
        for (int j = 0; j < COEFF_POLY_TERMS; j++) {
            // Index wraps around to keep the result polynomial the same size
            int k = (i + j) % COEFF_POLY_TERMS;
            result[k] += a[i] * b[j];
        }
    }
}

/* --- BENCHMARK SETUP --- */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_variables> <poly_degree1> <poly_degree2> <seed>\n", argv[0]);
        exit(1);
    }
    
    NUM_VARIABLES = atoi(argv[1]);
    POLY_DEGREE1 = atoi(argv[2]);
    POLY_DEGREE2 = atoi(argv[3]);
    SEED = (unsigned int)atoi(argv[4]);

    if (NUM_VARIABLES <= 1 || POLY_DEGREE1 < 0 || POLY_DEGREE2 < 0) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }
    
    // The number of terms in the coefficient polynomials is a function of the
    // number of "other" variables. This models the exponential growth.
    // For tuning simplicity, we use a linear factor.
    COEFF_POLY_TERMS = 16 * NUM_VARIABLES;

    mt_seed(SEED);

    // Allocate pointers for coefficient polynomials of P and Q
    p_coeffs = (double**)malloc((POLY_DEGREE1 + 1) * sizeof(double*));
    q_coeffs = (double**)malloc((POLY_DEGREE2 + 1) * sizeof(double*));
    if (!p_coeffs || !q_coeffs) {
        fprintf(stderr, "p_coeffs/q_coeffs malloc failed\n");
        exit(1);
    }
    
    // Allocate and initialize each coefficient polynomial for P
    for (int i = 0; i <= POLY_DEGREE1; i++) {
        p_coeffs[i] = (double*)malloc(COEFF_POLY_TERMS * sizeof(double));
        if (!p_coeffs[i]) exit(1);
        for (int j = 0; j < COEFF_POLY_TERMS; j++) {
            p_coeffs[i][j] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        }
    }

    // Allocate and initialize each coefficient polynomial for Q
    for (int i = 0; i <= POLY_DEGREE2; i++) {
        q_coeffs[i] = (double*)malloc(COEFF_POLY_TERMS * sizeof(double));
        if (!q_coeffs[i]) exit(1);
        for (int j = 0; j < COEFF_POLY_TERMS; j++) {
            q_coeffs[i][j] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        }
    }

    // Allocate helper structures
    zero_poly = (double*)calloc(COEFF_POLY_TERMS, sizeof(double));

    for (int i = 0; i < TEMP_POLY_POOL_SIZE; i++) {
        temp_poly_pool[i] = (double*)malloc(COEFF_POLY_TERMS * sizeof(double));
    }
    
    int max_dim = POLY_DEGREE1 + POLY_DEGREE2;
    if (max_dim > 0) {
        sub_matrix_scratch = (double***)malloc((max_dim - 1) * sizeof(double**));
        for (int i = 0; i < max_dim - 1; i++) {
            sub_matrix_scratch[i] = (double**)malloc((max_dim - 1) * sizeof(double*));
        }
    } else {
        sub_matrix_scratch = NULL;
    }
}

/* --- COMPUTATION --- */
void run_computation() {
    int d1 = POLY_DEGREE1;
    int d2 = POLY_DEGREE2;
    int n = d1 + d2;

    if (n == 0) {
        final_result_sum = 1.0; // P and Q are constants, resultant is 1 if non-zero
        return;
    }
    if (n*2 >= TEMP_POLY_POOL_SIZE) {
        fprintf(stderr, "Matrix size too large for temp pool.\n");
        exit(1);
    }
    
    // Create the Sylvester Matrix S (a matrix of pointers to polynomials)
    double*** S = (double***)malloc(n * sizeof(double**));
    for (int i = 0; i < n; i++) {
        S[i] = (double**)malloc(n * sizeof(double*));
    }

    // Populate the first d2 rows from polynomial P's coefficients
    for (int i = 0; i < d2; i++) {
        for (int j = 0; j < n; j++) {
            int coeff_idx = d1 - (j - i);
            if (j >= i && j <= i + d1) {
                S[i][j] = p_coeffs[coeff_idx];
            } else {
                S[i][j] = zero_poly;
            }
        }
    }

    // Populate the next d1 rows from polynomial Q's coefficients
    for (int i = 0; i < d1; i++) {
        int row = i + d2;
        for (int j = 0; j < n; j++) {
             int coeff_idx = d2 - (j - i);
            if (j >= i && j <= i + d2) {
                S[row][j] = q_coeffs[coeff_idx];
            } else {
                S[row][j] = zero_poly;
            }
        }
    }
    
    // Reset temp pool index and compute determinant
    temp_poly_pool_idx = 0;
    double* resultant_poly = determinant(S, n);

    // To prevent dead code elimination, sum up the resultant polynomial coefficients
    double sum = 0.0;
    for (int i = 0; i < COEFF_POLY_TERMS; i++) {
        sum += resultant_poly[i];
    }
    final_result_sum = sum;

    // Free the Sylvester matrix structure
    for (int i = 0; i < n; i++) {
        free(S[i]);
    }
    free(S);
}

// Recursively compute the determinant of a matrix of polynomials using cofactor expansion.
// Returns a pointer to a polynomial from the temp_poly_pool.
// This function uses a simple stack-like allocation from the pool.
double* determinant(double*** M, int n) {
    // Base case: 1x1 matrix, determinant is the single element.
    // We can return a pointer to the original polynomial, no copy needed.
    if (n == 1) {
        return M[0][0];
    }

    // Get a buffer for the final result of this level's determinant.
    double* final_det = temp_poly_pool[temp_poly_pool_idx++];
    memset(final_det, 0, sizeof(double) * COEFF_POLY_TERMS);

    // Get a temporary buffer for intermediate products (M[0][j] * sub_det).
    double* term_prod = temp_poly_pool[temp_poly_pool_idx++];

    // Cofactor expansion along the first row
    for (int j = 0; j < n; j++) {
        // If the coefficient is a zero polynomial, the whole term is zero.
        if (M[0][j] == zero_poly) {
            continue;
        }

        // Create the submatrix by removing row 0 and column j.
        // We use the pre-allocated scratch space to avoid malloc in the loop.
        int sub_col = 0;
        for (int col = 0; col < n; col++) {
            if (col == j) continue;
            for (int row = 1; row < n; row++) {
                sub_matrix_scratch[row - 1][sub_col] = M[row][col];
            }
            sub_col++;
        }

        double* sub_det = determinant(sub_matrix_scratch, n - 1);
        
        poly_mult(M[0][j], sub_det, term_prod);

        // Accumulate into the final determinant for this level
        if (j % 2 == 1) { // - term
            poly_sub(final_det, term_prod, final_det);
        } else { // + term
            poly_add(final_det, term_prod, final_det);
        }
    }

    // "Free" the temporary product buffer by decrementing the pool index.
    // The final_det buffer is left for the caller.
    temp_poly_pool_idx--;

    // The recursion will naturally unwind the index for sub-determinant results.
    return final_det;
}


/* --- BENCHMARK CLEANUP --- */
void cleanup() {
    int max_dim = POLY_DEGREE1 + POLY_DEGREE2;
    if (sub_matrix_scratch != NULL) {
        for (int i = 0; i < max_dim - 1; i++) {
            free(sub_matrix_scratch[i]);
        }
        free(sub_matrix_scratch);
    }

    for (int i = 0; i < TEMP_POLY_POOL_SIZE; i++) {
        free(temp_poly_pool[i]);
    }
    free(zero_poly);
    
    for (int i = 0; i <= POLY_DEGREE1; i++) {
        free(p_coeffs[i]);
    }
    free(p_coeffs);

    for (int i = 0; i <= POLY_DEGREE2; i++) {
        free(q_coeffs[i]);
    }
    free(q_coeffs);
}


/* --- MAIN FUNCTION --- */
int main(int argc, char *argv[]) {
    struct timespec start, end;
    
    setup_benchmark(argc, argv);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print final accumulated result to stdout
    printf("%.6f\n", final_result_sum);
    
    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);
    
    return 0;
}
