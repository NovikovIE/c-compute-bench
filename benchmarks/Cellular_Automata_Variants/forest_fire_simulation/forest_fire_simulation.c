#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (Provided) ---
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

// --- Benchmark Specifics ---

// Cell states
#define EMPTY 0
#define TREE 1
#define BURNING 2

// Global state struct
typedef struct {
    int grid_size;
    int num_steps;
    double fire_probability;
    double growth_probability;
    unsigned int final_result;

    unsigned char *grid_a;
    unsigned char *grid_b;
} BenchmarkData;

BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <grid_size> <num_steps> <fire_probability> <growth_probability> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.grid_size = atoi(argv[1]);
    g_data.num_steps = atoi(argv[2]);
    g_data.fire_probability = atof(argv[3]);
    g_data.growth_probability = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if (g_data.grid_size <= 0 || g_data.num_steps <= 0) {
        fprintf(stderr, "FATAL: grid_size and num_steps must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    size_t grid_bytes = (size_t)g_data.grid_size * g_data.grid_size * sizeof(unsigned char);
    g_data.grid_a = (unsigned char *)malloc(grid_bytes);
    g_data.grid_b = (unsigned char *)malloc(grid_bytes);

    if (!g_data.grid_a || !g_data.grid_b) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize grid_a with a mix of empty and tree cells
    for (int i = 0; i < g_data.grid_size * g_data.grid_size; ++i) {
        if (((double)mt_rand() / (double)UINT32_MAX) < 0.5) {
            g_data.grid_a[i] = TREE;
        } else {
            g_data.grid_a[i] = EMPTY;
        }
    }
}

void run_computation() {
    int size = g_data.grid_size;
    
    unsigned char *current_grid = g_data.grid_a;
    unsigned char *next_grid = g_data.grid_b;

    uint32_t fire_thresh = (uint32_t)(g_data.fire_probability * UINT32_MAX);
    uint32_t growth_thresh = (uint32_t)(g_data.growth_probability * UINT32_MAX);

    for (int step = 0; step < g_data.num_steps; ++step) {
        for (int r = 0; r < size; ++r) {
            for (int c = 0; c < size; ++c) {
                int current_idx = r * size + c;
                unsigned char current_state = current_grid[current_idx];
                unsigned char next_state = current_state;

                switch (current_state) {
                    case BURNING:
                        next_state = EMPTY;
                        break;
                    case TREE: {
                        int neighbor_on_fire = 0;
                        // Check 8 neighbors with toroidal boundary conditions
                        for (int dr = -1; dr <= 1; ++dr) {
                            for (int dc = -1; dc <= 1; ++dc) {
                                if (dr == 0 && dc == 0) continue;
                                int nr = (r + dr + size) % size;
                                int nc = (c + dc + size) % size;
                                if (current_grid[nr * size + nc] == BURNING) {
                                    neighbor_on_fire = 1;
                                    break;
                                }
                            }
                            if (neighbor_on_fire) break;
                        }

                        if (neighbor_on_fire) {
                            next_state = BURNING;
                        } else {
                            if (mt_rand() < fire_thresh) {
                                next_state = BURNING;
                            }
                        }
                        break;
                    }
                    case EMPTY:
                        if (mt_rand() < growth_thresh) {
                            next_state = TREE;
                        }
                        break;
                }
                next_grid[current_idx] = next_state;
            }
        }

        unsigned char *temp = current_grid;
        current_grid = next_grid;
        next_grid = temp;
    }

    unsigned int tree_count = 0;
    for (int i = 0; i < size * size; ++i) {
        if (current_grid[i] == TREE) {
            tree_count++;
        }
    }
    g_data.final_result = tree_count;
}

void cleanup() {
    free(g_data.grid_a);
    free(g_data.grid_b);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%u\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
