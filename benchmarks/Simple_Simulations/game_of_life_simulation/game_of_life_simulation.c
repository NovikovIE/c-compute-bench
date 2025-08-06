#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Globals ---
int GRID_SIZE;
int NUM_GENERATIONS;
char* grid_current;
char* grid_next;
int final_live_cells;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <grid_size> <num_generations> <seed>\n", argv[0]);
        exit(1);
    }

    GRID_SIZE = atoi(argv[1]);
    NUM_GENERATIONS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (GRID_SIZE <= 0 || NUM_GENERATIONS <= 0) {
        fprintf(stderr, "Invalid grid_size or num_generations.\n");
        exit(1);
    }
    
    mt_seed(seed);

    size_t grid_bytes = (size_t)GRID_SIZE * GRID_SIZE * sizeof(char);
    grid_current = (char*)malloc(grid_bytes);
    grid_next = (char*)malloc(grid_bytes);

    if (!grid_current || !grid_next) {
        fprintf(stderr, "Failed to allocate memory for grids.\n");
        exit(1);
    }

    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        grid_current[i] = mt_rand() % 2;
    }
}

void run_computation() {
    for (int gen = 0; gen < NUM_GENERATIONS; gen++) {
        for (int y = 0; y < GRID_SIZE; y++) {
            for (int x = 0; x < GRID_SIZE; x++) {
                int live_neighbors = 0;
                for (int dy = -1; dy <= 1; dy++) {
                    for (int dx = -1; dx <= 1; dx++) {
                        if (dx == 0 && dy == 0) continue;

                        // Toroidal (wrapping) boundary conditions
                        int nx = (x + dx + GRID_SIZE) % GRID_SIZE;
                        int ny = (y + dy + GRID_SIZE) % GRID_SIZE;

                        if (grid_current[ny * GRID_SIZE + nx]) {
                            live_neighbors++;
                        }
                    }
                }

                int current_cell_index = y * GRID_SIZE + x;
                char current_state = grid_current[current_cell_index];
                char next_state = 0;

                if (current_state == 1) { // Live cell
                    if (live_neighbors == 2 || live_neighbors == 3) {
                        next_state = 1; // Survives
                    }
                } else { // Dead cell
                    if (live_neighbors == 3) {
                        next_state = 1; // Reproduction
                    }
                }
                grid_next[current_cell_index] = next_state;
            }
        }
        
        // Swap grids for the next generation
        char* temp = grid_current;
        grid_current = grid_next;
        grid_next = temp;
    }

    // Calculate final result to prevent dead code elimination
    int total_live_cells = 0;
    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        if (grid_current[i]) {
            total_live_cells++;
        }
    }
    final_live_cells = total_live_cells;
}

void cleanup() {
    free(grid_current);
    free(grid_next);
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

    // Print result to stdout
    printf("%d\n", final_live_cells);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}