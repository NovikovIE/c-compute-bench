#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA AND CONFIG ---

// A simple struct to represent a symbolic expression by counting its terms
typedef struct {
    unsigned long long num_terms;
} Expression;

int matrix_dim;
Expression** matrix;
unsigned long long final_result;

// --- FORWARD DECLARATIONS ---

Expression determinant(Expression** mat, int n);

// --- BENCHMARK IMPLEMENTATION ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_dim> <seed>\n", argv[0]);
        exit(1);
    }

    matrix_dim = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (matrix_dim <= 0) {
        fprintf(stderr, "FATAL: matrix_dim must be positive.\n");
        exit(1);
    }
    if (matrix_dim > 20) {
        // Factorial overflows unsigned long long beyond 20!
        fprintf(stderr, "FATAL: matrix_dim > 20 will overflow result.\n");
        exit(1);
    }

    matrix = (Expression**)malloc(matrix_dim * sizeof(Expression*));
    if (!matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed for matrix rows.\n");
        exit(1);
    }
    for (int i = 0; i < matrix_dim; i++) {
        matrix[i] = (Expression*)malloc(matrix_dim * sizeof(Expression));
        if (!matrix[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for matrix columns.\n");
            for(int j = 0; j < i; ++j) free(matrix[j]);
            free(matrix);
            exit(1);
        }
    }

    // Initialize matrix. Each entry is a single symbolic variable, hence has 1 term.
    // Using mt_rand() adds work but doesn't change the problem's nature.
    for (int i = 0; i < matrix_dim; i++) {
        for (int j = 0; j < matrix_dim; j++) {
            matrix[i][j].num_terms = 1;
        }
    }
}

// Recursively calculates the number of terms in the determinant of a symbolic matrix
// using cofactor expansion. This method has O(n!) complexity and demonstrates
// "expression swell", where the number of terms grows factorially.
Expression determinant(Expression** mat, int n) {
    if (n == 1) {
        return mat[0][0];
    }
    
    Expression total = {0};

    Expression** sub_matrix = (Expression**)malloc((n - 1) * sizeof(Expression*));
    for (int i = 0; i < n - 1; i++) {
        sub_matrix[i] = (Expression*)malloc((n - 1) * sizeof(Expression));
    }

    // Iterate through first row to calculate cofactors
    for (int p = 0; p < n; p++) {
        // Create the (n-1)x(n-1) submatrix for cofactor M_0p
        int sub_i = 0;
        for (int i = 1; i < n; i++) {
            int sub_j = 0;
            for (int j = 0; j < n; j++) {
                if (j == p) continue; // Skip column p
                sub_matrix[sub_i][sub_j] = mat[i][j];
                sub_j++;
            }
            sub_i++;
        }

        Expression det_sub = determinant(sub_matrix, n - 1);
        
        // Multiply(A, B) results in Terms(A) * Terms(B) terms.
        unsigned long long term_count = mat[0][p].num_terms * det_sub.num_terms;
        
        // In symbolic math with unique variables, Add/Subtract both sum the term counts.
        total.num_terms += term_count;
    }
    
    for (int i = 0; i < n - 1; i++) {
        free(sub_matrix[i]);
    }
    free(sub_matrix);

    return total;
}

void run_computation() {
    Expression result = determinant(matrix, matrix_dim);
    final_result = result.num_terms;
}

void cleanup() {
    for (int i = 0; i < matrix_dim; i++) {
        free(matrix[i]);
    }
    free(matrix);
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
    printf("%llu\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
