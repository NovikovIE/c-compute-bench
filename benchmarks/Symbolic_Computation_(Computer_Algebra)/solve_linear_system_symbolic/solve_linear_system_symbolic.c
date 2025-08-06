#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h> // For llabs

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Data Structures and Globals ---

// Parameters
int NUM_EQUATIONS;
int NUM_VARIABLES;
const int MAX_COEFF = 50;
const int MAX_DENOM = 10;

// Data structure for a rational number (fraction)
typedef struct {
    long long num; // numerator
    long long den; // denominator
} Rational;

// The augmented matrix [A|b]
Rational **augmented_matrix;

// Final result for printing
long long final_result_sum;

// --- Symbolic Arithmetic Helper Functions ---

// Greatest Common Divisor (Euclidean algorithm)
long long gcd(long long a, long long b) {
    a = llabs(a);
    b = llabs(b);
    while (b) {
        long long temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// Simplify a rational number to its canonical form
void simplify(Rational *r) {
    if (r->den == 0) {
        fprintf(stderr, "FATAL: Division by zero in rational number.\n");
        exit(1);
    }
    if (r->num == 0) {
        r->den = 1;
        return;
    }
    long long common_divisor = gcd(r->num, r->den);
    r->num /= common_divisor;
    r->den /= common_divisor;
    // Ensure denominator is positive for a canonical representation
    if (r->den < 0) {
        r->num = -r->num;
        r->den = -r->den;
    }
}

// Subtract two rational numbers
Rational subtract(Rational a, Rational b) {
    Rational result;
    result.num = a.num * b.den - b.num * a.den;
    result.den = a.den * b.den;
    simplify(&result);
    return result;
}

// Multiply two rational numbers
Rational multiply(Rational a, Rational b) {
    Rational result;
    result.num = a.num * b.num;
    result.den = a.den * b.den;
    simplify(&result);
    return result;
}

// Divide two rational numbers
Rational divide(Rational a, Rational b) {
    if (b.num == 0) {
        fprintf(stderr, "FATAL: Symbolic division by zero.\n");
        exit(1);
    }
    Rational result;
    result.num = a.num * b.den;
    result.den = a.den * b.num;
    simplify(&result);
    return result;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_equations> <num_variables> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_EQUATIONS = atoi(argv[1]);
    NUM_VARIABLES = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (NUM_EQUATIONS <= 0 || NUM_VARIABLES <= 0) {
        fprintf(stderr, "FATAL: Number of equations and variables must be positive.\n");
        exit(1);
    }

    // Allocate the augmented matrix: NUM_EQUATIONS rows, NUM_VARIABLES + 1 columns
    int num_cols = NUM_VARIABLES + 1;
    augmented_matrix = (Rational **)malloc(NUM_EQUATIONS * sizeof(Rational *));
    if (!augmented_matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed for matrix rows.\n");
        exit(1);
    }
    for (int i = 0; i < NUM_EQUATIONS; i++) {
        augmented_matrix[i] = (Rational *)malloc(num_cols * sizeof(Rational));
        if (!augmented_matrix[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for matrix columns.\n");
            for (int k = 0; k < i; ++k) free(augmented_matrix[k]);
            free(augmented_matrix);
            exit(1);
        }
    }

    // Populate the matrix with random rational numbers
    for (int i = 0; i < NUM_EQUATIONS; i++) {
        for (int j = 0; j < num_cols; j++) {
            long long num = (mt_rand() % (2 * MAX_COEFF)) - MAX_COEFF;
            long long den = (mt_rand() % (MAX_DENOM - 1)) + 1;
            augmented_matrix[i][j].num = num;
            augmented_matrix[i][j].den = den;
            simplify(&augmented_matrix[i][j]);
        }
    }

    final_result_sum = 0;
}

void run_computation() {
    int pivot_row = 0;
    int pivot_col = 0;
    int num_cols = NUM_VARIABLES + 1;

    // Forward elimination (Gaussian elimination) to get row echelon form
    while (pivot_row < NUM_EQUATIONS && pivot_col < NUM_VARIABLES) {
        // Find pivot: row with largest absolute value in the current column
        int max_row = pivot_row;
        for (int i = pivot_row + 1; i < NUM_EQUATIONS; i++) {
            double val_current = (double)llabs(augmented_matrix[i][pivot_col].num) / augmented_matrix[i][pivot_col].den;
            double val_max = (double)llabs(augmented_matrix[max_row][pivot_col].num) / augmented_matrix[max_row][pivot_col].den;
            if (val_current > val_max) {
                max_row = i;
            }
        }

        // Swap pivot row to bring the best pivot up
        if (max_row != pivot_row) {
            Rational *temp_row = augmented_matrix[pivot_row];
            augmented_matrix[pivot_row] = augmented_matrix[max_row];
            augmented_matrix[max_row] = temp_row;
        }

        // If pivot is zero, this column has no pivot, move to next column
        if (augmented_matrix[pivot_row][pivot_col].num == 0) {
            pivot_col++;
            continue;
        }

        // Eliminate entries *below* the pivot in this column
        for (int i = pivot_row + 1; i < NUM_EQUATIONS; i++) {
            if (augmented_matrix[i][pivot_col].num != 0) {
                Rational factor = divide(augmented_matrix[i][pivot_col], augmented_matrix[pivot_row][pivot_col]);
                for (int j = pivot_col; j < num_cols; j++) {
                    Rational term = multiply(factor, augmented_matrix[pivot_row][j]);
                    augmented_matrix[i][j] = subtract(augmented_matrix[i][j], term);
                }
            }
        }
        pivot_row++;
        pivot_col++;
    }

    // Accumulate a result to prevent dead code elimination.
    // We sum the numerators of the final augmented column vector.
    for (int i = 0; i < NUM_EQUATIONS; i++) {
        final_result_sum += augmented_matrix[i][NUM_VARIABLES].num;
    }
}

void cleanup() {
    if (augmented_matrix) {
        for (int i = 0; i < NUM_EQUATIONS; i++) {
            if (augmented_matrix[i]) {
                free(augmented_matrix[i]);
            }
        }
        free(augmented_matrix);
    }
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
    printf("%lld\n", final_result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
