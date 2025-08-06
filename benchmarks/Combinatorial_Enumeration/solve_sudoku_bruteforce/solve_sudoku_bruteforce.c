#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

// --- Mersenne Twister (MT19937) start ---
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
// --- Mersenne Twister end ---

// --- Global Benchmark Data ---
typedef struct {
    int board_size;
    int sr_board_size; // Square root of board_size
    int num_prefilled_cells;
    int **board;
    unsigned long long solution_count;
} SudokuBenchmark;

SudokuBenchmark G_Args;

// --- Sudoku Helper Functions ---

bool is_safe(int row, int col, int num) {
    int N = G_Args.board_size;
    int SRN = G_Args.sr_board_size;

    // Check row
    for (int x = 0; x < N; x++) {
        if (G_Args.board[row][x] == num) return false;
    }

    // Check column
    for (int x = 0; x < N; x++) {
        if (G_Args.board[x][col] == num) return false;
    }

    // Check 3x3 box
    int startRow = row - row % SRN;
    int startCol = col - col % SRN;
    for (int i = 0; i < SRN; i++) {
        for (int j = 0; j < SRN; j++) {
            if (G_Args.board[i + startRow][j + startCol] == num) return false;
        }
    }

    return true;
}

bool generate_full_board(int row, int col) {
    int N = G_Args.board_size;
    if (row == N - 1 && col == N) return true;

    if (col == N) {
        row++;
        col = 0;
    }

    if (G_Args.board[row][col] != 0) return generate_full_board(row, col + 1);

    int nums[N];
    for(int i=0; i<N; i++) nums[i] = i+1;
    // Shuffle numbers to get a random valid board
    for (int i = N - 1; i > 0; i--) {
        int j = mt_rand() % (i + 1);
        int temp = nums[i];
        nums[i] = nums[j];
        nums[j] = temp;
    }

    for (int i = 0; i < N; i++) {
        int num = nums[i];
        if (is_safe(row, col, num)) {
            G_Args.board[row][col] = num;
            if (generate_full_board(row, col + 1)) return true;
        }
    }

    G_Args.board[row][col] = 0; // Backtrack
    return false;
}

void count_solutions_recursive(int row, int col) {
    int N = G_Args.board_size;
    
    if (row == N - 1 && col == N) {
        G_Args.solution_count++;
        return;
    }

    if (col == N) {
        row++;
        col = 0;
    }

    if (G_Args.board[row][col] != 0) {
        count_solutions_recursive(row, col + 1);
        return;
    }

    for (int num = 1; num <= N; num++) {
        if (is_safe(row, col, num)) {
            G_Args.board[row][col] = num;
            count_solutions_recursive(row, col + 1);
        }
    }

    G_Args.board[row][col] = 0; // Backtrack
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <board_size> <num_prefilled_cells> <seed>\n", argv[0]);
        exit(1);
    }

    G_Args.board_size = atoi(argv[1]);
    G_Args.num_prefilled_cells = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    G_Args.sr_board_size = (int)sqrt((double)G_Args.board_size);
    if (G_Args.sr_board_size * G_Args.sr_board_size != G_Args.board_size) {
        fprintf(stderr, "FATAL: board_size must be a perfect square.\n");
        exit(1);
    }

    int N = G_Args.board_size;
    int total_cells = N * N;
    if (G_Args.num_prefilled_cells < 0 || G_Args.num_prefilled_cells >= total_cells) {
        fprintf(stderr, "FATAL: num_prefilled_cells must be between 0 and %d.\n", total_cells - 1);
        exit(1);
    }

    // Allocate and initialize board to 0
    G_Args.board = (int **)malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        G_Args.board[i] = (int *)calloc(N, sizeof(int));
    }

    // 1. Generate a fully solved Sudoku board
    if (!generate_full_board(0, 0)) {
        fprintf(stderr, "FATAL: Could not generate a solved Sudoku board.\n");
        exit(1);
    }

    // 2. Poke holes in the board to create the puzzle
    int *indices = (int *)malloc(total_cells * sizeof(int));
    for (int i = 0; i < total_cells; i++) {
        indices[i] = i;
    }
    
    // Fisher-Yates shuffle
    for (int i = total_cells - 1; i > 0; i--) {
        int j = mt_rand() % (i + 1);
        int temp = indices[i];
        indices[i] = indices[j];
        indices[j] = temp;
    }
    
    int cells_to_remove = total_cells - G_Args.num_prefilled_cells;
    for (int i = 0; i < cells_to_remove; i++) {
        int cell_idx = indices[i];
        int row = cell_idx / N;
        int col = cell_idx % N;
        G_Args.board[row][col] = 0;
    }

    free(indices);
    G_Args.solution_count = 0;
}

void run_computation() {
    count_solutions_recursive(0, 0);
}

void cleanup() {
    for (int i = 0; i < G_Args.board_size; i++) {
        free(G_Args.board[i]);
    }
    free(G_Args.board);
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
    printf("%llu\n", G_Args.solution_count);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
