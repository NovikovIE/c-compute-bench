#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK STATE ---
typedef struct {
    int grid_dim_x;
    int grid_dim_y;
    int num_time_steps;

    double **grid_prev;
    double **grid_curr;
    double **grid_new;

    // Keep original pointers for cleanup
    double *__data_block1;
    double *__data_block2;
    double *__data_block3;
    double **__grid1;
    double **__grid2;
    double **__grid3;

    double final_checksum;
} BenchmarkState;

BenchmarkState state;

// --- HELPER FUNCTIONS ---
double** allocate_grid(int dim_x, int dim_y, double** data_block_ptr) {
    *data_block_ptr = (double*)malloc(dim_x * dim_y * sizeof(double));
    if (*data_block_ptr == NULL) {
        fprintf(stderr, "Failed to allocate data block\n");
        exit(1);
    }

    double** grid = (double**)malloc(dim_y * sizeof(double*));
    if (grid == NULL) {
        fprintf(stderr, "Failed to allocate grid pointers\n");
        free(*data_block_ptr);
        exit(1);
    }

    for (int i = 0; i < dim_y; ++i) {
        grid[i] = &((*data_block_ptr)[i * dim_x]);
    }
    return grid;
}

void initialize_grid(double** grid, int dim_x, int dim_y) {
    for (int i = 0; i < dim_y; ++i) {
        for (int j = 0; j < dim_x; ++j) {
             // Set a small random initial displacement, but keep boundaries zero
            if (i == 0 || i == dim_y - 1 || j == 0 || j == dim_x - 1) {
                grid[i][j] = 0.0;
            } else {
                grid[i][j] = (mt_rand() / 4294967295.0 - 0.5) * 0.01;
            }
        }
    }
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_dim_x grid_dim_y num_time_steps seed\n", argv[0]);
        exit(1);
    }

    state.grid_dim_x = atoi(argv[1]);
    state.grid_dim_y = atoi(argv[2]);
    state.num_time_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate three grids for previous, current, and next time steps
    state.__grid1 = allocate_grid(state.grid_dim_x, state.grid_dim_y, &state.__data_block1);
    state.__grid2 = allocate_grid(state.grid_dim_x, state.grid_dim_y, &state.__data_block2);
    state.__grid3 = allocate_grid(state.grid_dim_x, state.grid_dim_y, &state.__data_block3);
    
    state.grid_prev = state.__grid1;
    state.grid_curr = state.__grid2;
    state.grid_new = state.__grid3;

    // Initialize grids with initial conditions
    initialize_grid(state.grid_prev, state.grid_dim_x, state.grid_dim_y);
    initialize_grid(state.grid_curr, state.grid_dim_x, state.grid_dim_y);
    
    state.final_checksum = 0.0;
}

void run_computation() {
    const double c = 1.0;  // Wave speed
    const double dt = 0.1; // Time step
    const double dx = 1.0; // Grid spacing in x
    const double dy = 1.0; // Grid spacing in y
    const double C_x = (c * dt / dx) * (c * dt / dx);
    const double C_y = (c * dt / dy) * (c * dt / dy);
    const double stencil_C0 = 2.0 - 2.0 * C_x - 2.0 * C_y;

    for (int t = 0; t < state.num_time_steps; ++t) {
        for (int i = 1; i < state.grid_dim_y - 1; ++i) {
            for (int j = 1; j < state.grid_dim_x - 1; ++j) {
                double laplacian = C_x * (state.grid_curr[i][j+1] + state.grid_curr[i][j-1]) +
                                   C_y * (state.grid_curr[i+1][j] + state.grid_curr[i-1][j]);

                state.grid_new[i][j] = stencil_C0 * state.grid_curr[i][j] 
                                     + laplacian 
                                     - state.grid_prev[i][j];
            }
        }

        // Swap pointers for the next time step
        double **temp = state.grid_prev;
        state.grid_prev = state.grid_curr;
        state.grid_curr = state.grid_new;
        state.grid_new = temp;
    }

    // Calculate a checksum to prevent dead code elimination
    double checksum = 0.0;
    for (int i = 0; i < state.grid_dim_y; ++i) {
        for (int j = 0; j < state.grid_dim_x; ++j) {
            checksum += state.grid_curr[i][j];
        }
    }
    state.final_checksum = checksum;
}

void cleanup() {
    free(state.__data_block1);
    free(state.__grid1);
    free(state.__data_block2);
    free(state.__grid2);
    free(state.__data_block3);
    free(state.__grid3);
}

// --- MAIN --- 
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", state.final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
