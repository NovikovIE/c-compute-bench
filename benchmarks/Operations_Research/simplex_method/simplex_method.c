#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <float.h>

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
// --- END Mersenne Twister ---

// --- Global Benchmark Data ---
int NUM_VARIABLES;
int NUM_CONSTRAINTS;
double** tableau;
long long final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_variables> <num_constraints> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_VARIABLES = atoi(argv[1]);
    NUM_CONSTRAINTS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_VARIABLES <= 0 || NUM_CONSTRAINTS <= 0) {
        fprintf(stderr, "FATAL: Number of variables and constraints must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    int rows = NUM_CONSTRAINTS + 1;
    int cols = NUM_VARIABLES + NUM_CONSTRAINTS + 1;

    tableau = (double**)malloc(rows * sizeof(double*));
    if (tableau == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for tableau rows.\n");
        exit(1);
    }
    for (int i = 0; i < rows; i++) {
        tableau[i] = (double*)malloc(cols * sizeof(double));
        if (tableau[i] == NULL) {
            fprintf(stderr, "FATAL: Memory allocation failed for tableau columns.\n");
            exit(1);
        }
    }

    // Generate a standard form LP problem: max c'x s.t. Ax <= b, x >= 0
    // Tableau structure:
    // [  A | I | b ]
    // [ -c | 0 | 0 ]

    // Generate A (constraint coefficients)
    for (int i = 0; i < NUM_CONSTRAINTS; i++) {
        for (int j = 0; j < NUM_VARIABLES; j++) {
            tableau[i][j] = ((double)mt_rand() / UINT32_MAX) * 9.0 + 1.0;
        }
    }

    // Set I (identity matrix for slack variables)
    for (int i = 0; i < NUM_CONSTRAINTS; i++) {
        for (int j = 0; j < NUM_CONSTRAINTS; j++) {
            tableau[i][NUM_VARIABLES + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Generate b (constraint bounds), ensuring they are non-negative for initial feasibility
    for (int i = 0; i < NUM_CONSTRAINTS; i++) {
        tableau[i][cols - 1] = ((double)mt_rand() / UINT32_MAX) * 100.0 + 50.0;
    }

    // Generate -c (objective function), ensuring some are negative to start the algorithm
    for (int j = 0; j < NUM_VARIABLES; j++) {
        tableau[rows - 1][j] = -(((double)mt_rand() / UINT32_MAX) * 9.0 + 1.0);
    }
    // Set rest of objective row to 0
    for (int j = NUM_VARIABLES; j < cols; j++) {
        tableau[rows - 1][j] = 0.0;
    }
}

void run_computation() {
    int rows = NUM_CONSTRAINTS + 1;
    int cols = NUM_VARIABLES + NUM_CONSTRAINTS + 1;
    
    int max_iterations = 3 * (NUM_VARIABLES + NUM_CONSTRAINTS);
    int iterations = 0;

    while (iterations < max_iterations) {
        // Find pivot column: most negative entry in objective row
        int pivot_col = -1;
        double min_obj_val = -1e-9; // Tolerance for floating point zero
        
        for (int j = 0; j < cols - 1; j++) {
            if (tableau[rows - 1][j] < min_obj_val) {
                min_obj_val = tableau[rows - 1][j];
                pivot_col = j;
            }
        }

        // If no negative value found, the current solution is optimal
        if (pivot_col == -1) {
            break; 
        }

        // Find pivot row: minimum ratio test
        int pivot_row = -1;
        double min_ratio = DBL_MAX;

        for (int i = 0; i < rows - 1; i++) {
            if (tableau[i][pivot_col] > 1e-9) { // Check for positive divisor
                double ratio = tableau[i][cols - 1] / tableau[i][pivot_col];
                if (ratio >= 0 && ratio < min_ratio) {
                    min_ratio = ratio;
                    pivot_row = i;
                }
            }
        }

        // If no suitable pivot row, the problem is unbounded
        if (pivot_row == -1) {
            break;
        }

        // Perform pivot operation
        double pivot_element = tableau[pivot_row][pivot_col];

        // 1. Normalize the pivot row
        for (int j = 0; j < cols; j++) {
            tableau[pivot_row][j] /= pivot_element;
        }

        // 2. Clear out the pivot column in all other rows
        for (int i = 0; i < rows; i++) {
            if (i != pivot_row) {
                double factor = tableau[i][pivot_col];
                for (int j = 0; j < cols; j++) {
                    tableau[i][j] -= factor * tableau[pivot_row][j];
                }
            }
        }
        
        iterations++;
    }

    final_result = iterations;
}

void cleanup() {
    int rows = NUM_CONSTRAINTS + 1;
    for (int i = 0; i < rows; i++) {
        free(tableau[i]);
    }
    free(tableau);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
