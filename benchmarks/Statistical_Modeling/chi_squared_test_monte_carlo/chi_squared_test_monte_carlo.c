#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Verbatim Mersenne Twister ---
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
// --- End Verbatim Mersenne Twister ---

// --- Benchmark Globals ---
int num_rows;
int num_cols;
int num_simulations;
int num_swaps_per_simulation;

// Data
int** observed_table;
int** simulation_table;
double** expected_values;
double* row_totals;
double* col_totals;
double grand_total;

// Result
int p_value_counter;


// --- Helper Functions ---
void* safe_malloc(size_t size) {
    void* ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "Fatal: malloc failed.\n");
        exit(1);
    }
    return ptr;
}

int** allocate_int_matrix(int rows, int cols) {
    int** matrix = (int**)safe_malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; ++i) {
        matrix[i] = (int*)safe_malloc(cols * sizeof(int));
    }
    return matrix;
}

double** allocate_double_matrix(int rows, int cols) {
    double** matrix = (double**)safe_malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; ++i) {
        matrix[i] = (double*)safe_malloc(cols * sizeof(double));
    }
    return matrix;
}

void free_matrix(void** matrix, int rows) {
    if (!matrix) return;
    for (int i = 0; i < rows; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}

double calculate_chi_squared(int** table) {
    double chi_sq_stat = 0.0;
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            // All expected values are pre-calculated to be > 0
            double diff = table[i][j] - expected_values[i][j];
            chi_sq_stat += (diff * diff) / expected_values[i][j];
        }
    }
    return chi_sq_stat;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_rows num_cols num_simulations num_swaps_per_simulation seed\n", argv[0]);
        exit(1);
    }

    num_rows = atoi(argv[1]);
    num_cols = atoi(argv[2]);
    num_simulations = atoi(argv[3]);
    num_swaps_per_simulation = atoi(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);
    
    mt_seed(seed);

    // Allocate memory
    observed_table = allocate_int_matrix(num_rows, num_cols);
    simulation_table = allocate_int_matrix(num_rows, num_cols);
    expected_values = allocate_double_matrix(num_rows, num_cols);
    row_totals = (double*)safe_malloc(num_rows * sizeof(double));
    col_totals = (double*)safe_malloc(num_cols * sizeof(double));

    // Initialize row/column totals to zero
    for(int i=0; i<num_rows; ++i) row_totals[i] = 0.0;
    for(int j=0; j<num_cols; ++j) col_totals[j] = 0.0;
    grand_total = 0.0;
    
    // Populate observed table with random data and calculate totals
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            // +1 to ensure all cell values, and thus row/col/grand totals are positive.
            observed_table[i][j] = (mt_rand() % 10) + 1;
            row_totals[i] += observed_table[i][j];
            col_totals[j] += observed_table[i][j];
        }
    }
    
    for(int i=0; i<num_rows; ++i) {
        grand_total += row_totals[i];
    }
    
    // Pre-calculate expected values
    if (grand_total == 0) {
        fprintf(stderr, "Fatal: grand_total is zero, cannot compute expected values.\n");
        exit(1);
    }
    for (int i = 0; i < num_rows; ++i) {
        for (int j = 0; j < num_cols; ++j) {
            expected_values[i][j] = (row_totals[i] * col_totals[j]) / grand_total;
        }
    }
}

void run_computation() {
    double observed_chi_sq = calculate_chi_squared(observed_table);
    p_value_counter = 0;

    for (int s = 0; s < num_simulations; ++s) {
        // Copy observed to simulation table for shuffling
        for (int i = 0; i < num_rows; ++i) {
            for (int j = 0; j < num_cols; ++j) {
                simulation_table[i][j] = observed_table[i][j];
            }
        }

        // Shuffle the table using swaps that preserve row/col sums
        for (int k = 0; k < num_swaps_per_simulation; ++k) {
            uint32_t r1 = mt_rand() % num_rows;
            uint32_t r2 = mt_rand() % num_rows;
            uint32_t c1 = mt_rand() % num_cols;
            uint32_t c2 = mt_rand() % num_cols;
            
            if (r1 == r2 || c1 == c2) continue;
            
            // Perform a "Diaconis-Gangolli" move
            if (simulation_table[r1][c1] > 0 && simulation_table[r2][c2] > 0) {
                simulation_table[r1][c1]--;
                simulation_table[r1][c2]++;
                simulation_table[r2][c1]++;
                simulation_table[r2][c2]--;
            }
        }

        double simulated_chi_sq = calculate_chi_squared(simulation_table);
        if (simulated_chi_sq >= observed_chi_sq) {
            p_value_counter++;
        }
    }
}

void cleanup() {
    free_matrix((void**)observed_table, num_rows);
    free_matrix((void**)simulation_table, num_rows);
    free_matrix((void**)expected_values, num_rows);
    free(row_totals);
    free(col_totals);
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
    printf("%d\n", p_value_counter);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
