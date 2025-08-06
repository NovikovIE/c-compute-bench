#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// DO NOT MODIFY - Mersenne Twister (MT19937) Generator
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
// End of Mersenne Twister

// --- Benchmark Data ---
typedef struct {
    char *grid;
    long grid_size;
    long num_steps;
    int ant_x;
    int ant_y;
    int ant_dir; // 0: Up, 1: Right, 2: Down, 3: Left
    long final_result;
} BenchmarkData;

BenchmarkData data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <grid_size> <num_steps> <seed>\n", argv[0]);
        exit(1);
    }

    data.grid_size = atol(argv[1]);
    data.num_steps = atol(argv[2]);
    uint32_t seed = atoi(argv[3]);
    
    mt_seed(seed);

    if (data.grid_size <= 0 || data.num_steps <= 0) {
        fprintf(stderr, "Invalid arguments: grid_size and num_steps must be positive.\n");
        exit(1);
    }

    data.grid = (char*)malloc(data.grid_size * data.grid_size * sizeof(char));
    if (!data.grid) {
        fprintf(stderr, "Failed to allocate memory for the grid.\n");
        exit(1);
    }
    memset(data.grid, 0, data.grid_size * data.grid_size * sizeof(char));

    data.ant_x = data.grid_size / 2;
    data.ant_y = data.grid_size / 2;
    data.ant_dir = 0; // Start facing Up
    data.final_result = 0;
}

void run_computation() {
    // LUTs for direction changes based on ant_dir: 0:U, 1:R, 2:D, 3:L
    const int dx[] = {0, 1, 0, -1};
    const int dy[] = {-1, 0, 1, 0};

    for (long i = 0; i < data.num_steps; ++i) {
        long current_pos_idx = (long)data.ant_y * data.grid_size + data.ant_x;
        
        if (data.grid[current_pos_idx] == 0) { // White square
            data.ant_dir = (data.ant_dir + 1) % 4; // Turn right 90 degrees
            data.grid[current_pos_idx] = 1;       // Flip color to black
        } else { // Black square
            data.ant_dir = (data.ant_dir + 3) % 4; // Turn left 90 degrees
            data.grid[current_pos_idx] = 0;       // Flip color to white
        }

        data.ant_x += dx[data.ant_dir];
        data.ant_y += dy[data.ant_dir];

        // Wrap around grid boundaries (toroidal grid)
        if (data.ant_x < 0) data.ant_x = data.grid_size - 1;
        if (data.ant_x >= data.grid_size) data.ant_x = 0;
        if (data.ant_y < 0) data.ant_y = data.grid_size - 1;
        if (data.ant_y >= data.grid_size) data.ant_y = 0;
    }
    
    long black_cells = 0;
    long grid_area = data.grid_size * data.grid_size;
    for (long i = 0; i < grid_area; ++i) {
        if (data.grid[i] == 1) {
            black_cells++;
        }
    }
    data.final_result = black_cells;
}

void cleanup() {
    free(data.grid);
    data.grid = NULL;
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%ld\n", data.final_result);

    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
