#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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
// --- End of MT19937 ---

// --- Benchmark Data and Globals ---
typedef struct {
    int grid_size;
    int num_states;
    int num_steps;
    uint8_t **grid_a;
    uint8_t **grid_b;
    long long final_result;
} BenchmarkData;

BenchmarkData bench_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_size num_states num_steps seed\n", argv[0]);
        exit(1);
    }

    bench_data.grid_size = atoi(argv[1]);
    bench_data.num_states = atoi(argv[2]);
    bench_data.num_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    if (bench_data.grid_size <= 0 || bench_data.num_states <= 1 || bench_data.num_steps <= 0) {
        fprintf(stderr, "FATAL: grid_size, num_states, and num_steps must be positive integers.\n");
        exit(1);
    }
    if (bench_data.num_states > 256) {
        fprintf(stderr, "FATAL: num_states cannot exceed 256 for uint8_t grid.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate two grids for double-buffering
    bench_data.grid_a = (uint8_t **)malloc(bench_data.grid_size * sizeof(uint8_t *));
    bench_data.grid_b = (uint8_t **)malloc(bench_data.grid_size * sizeof(uint8_t *));
    if (!bench_data.grid_a || !bench_data.grid_b) {
        fprintf(stderr, "FATAL: Memory allocation failed for grid pointers.\n");
        exit(1);
    }

    for (int i = 0; i < bench_data.grid_size; ++i) {
        bench_data.grid_a[i] = (uint8_t *)malloc(bench_data.grid_size * sizeof(uint8_t));
        bench_data.grid_b[i] = (uint8_t *)malloc(bench_data.grid_size * sizeof(uint8_t));
        if (!bench_data.grid_a[i] || !bench_data.grid_b[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for grid rows.\n");
            exit(1);
        }
    }

    // Initialize the primary grid with random states
    for (int i = 0; i < bench_data.grid_size; ++i) {
        for (int j = 0; j < bench_data.grid_size; ++j) {
            bench_data.grid_a[i][j] = mt_rand() % bench_data.num_states;
        }
    }
}

void run_computation() {
    uint8_t **current_grid = bench_data.grid_a;
    uint8_t **next_grid = bench_data.grid_b;
    
    int size = bench_data.grid_size;
    int n_states = bench_data.num_states;

    // Moore neighborhood (8 neighbors)
    int dx[] = {-1, -1, -1, 0, 0, 1, 1, 1};
    int dy[] = {-1, 0, 1, -1, 1, -1, 0, 1};

    for (int step = 0; step < bench_data.num_steps; ++step) {
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                uint8_t current_state = current_grid[i][j];
                uint8_t next_state_if_changed = (current_state + 1) % n_states;
                int found_neighbor = 0;

                for (int k = 0; k < 8; ++k) {
                    // Toroidal (wrapping) boundary conditions
                    int ni = (i + dy[k] + size) % size;
                    int nj = (j + dx[k] + size) % size;
                    if (current_grid[ni][nj] == next_state_if_changed) {
                        found_neighbor = 1;
                        break;
                    }
                }

                if (found_neighbor) {
                    next_grid[i][j] = next_state_if_changed;
                } else {
                    next_grid[i][j] = current_state;
                }
            }
        }

        // Swap grids
        uint8_t **temp = current_grid;
        current_grid = next_grid;
        next_grid = temp;
    }
    
    // Accumulate a final result to prevent dead code elimination
    long long sum = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            sum += current_grid[i][j];
        }
    }
    bench_data.final_result = sum;
}

void cleanup() {
    for (int i = 0; i < bench_data.grid_size; ++i) {
        free(bench_data.grid_a[i]);
        free(bench_data.grid_b[i]);
    }
    free(bench_data.grid_a);
    free(bench_data.grid_b);
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
    printf("%lld\n", bench_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
