#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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
// --- End of Mersenne Twister ---

// Global structure to hold benchmark data
typedef struct {
    int board_size;
    unsigned long long solution_count;
    int* queens; // Array to store column of queen for each row
} benchmark_data_t;

static benchmark_data_t g_data;

// Helper function for the recursive N-Queens solver
static void solve_queens_recursive(int row) {
    // Base case: If all queens are placed, we found a solution
    if (row == g_data.board_size) {
        g_data.solution_count++;
        return;
    }

    // Try placing a queen in each column of the current row
    for (int col = 0; col < g_data.board_size; col++) {
        int is_safe = 1;

        // Check for conflicts with queens in previous rows
        for (int prev_row = 0; prev_row < row; prev_row++) {
            // Check if the new queen is in the same column or on the same diagonal
            if (g_data.queens[prev_row] == col || 
                abs(g_data.queens[prev_row] - col) == (row - prev_row)) {
                is_safe = 0;
                break;
            }
        }

        // If the position is safe, place the queen and recurse
        if (is_safe) {
            g_data.queens[row] = col;
            solve_queens_recursive(row + 1);
        }
    }
}


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <board_size> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.board_size = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if(g_data.board_size <= 0) {
        fprintf(stderr, "Error: board_size must be a positive integer.\n");
        exit(1);
    }
    
    // Seed the random number generator (as required, though not used in this deterministic benchmark)
    mt_seed(seed);

    // Initialize benchmark data
    g_data.solution_count = 0;
    g_data.queens = (int*)malloc(g_data.board_size * sizeof(int));
    if (g_data.queens == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for queens array.\n");
        exit(1);
    }
}

void run_computation() {
    // Start the recursive backtracking from the first row (row 0)
    solve_queens_recursive(0);
}

void cleanup() {
    free(g_data.queens);
    g_data.queens = NULL;
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print the final result (total number of solutions) to stdout
    printf("%llu\n", g_data.solution_count);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
