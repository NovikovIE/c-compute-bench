#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// Benchmark state
// Wireworld cell states:
// 0: Empty, 1: Electron Head, 2: Electron Tail, 3: Conductor
int *grid_current;
int *grid_next;
int grid_width;
int grid_height;
int num_steps;
long long final_result; // Use long long for checksum to avoid overflow

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <grid_width> <grid_height> <num_steps> <seed>\n", argv[0]);
        exit(1);
    }

    grid_width = atoi(argv[1]);
    grid_height = atoi(argv[2]);
    num_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    if (grid_width <= 0 || grid_height <= 0 || num_steps <= 0) {
        fprintf(stderr, "Error: grid dimensions and steps must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    size_t grid_size = (size_t)grid_width * grid_height;
    grid_current = (int *)malloc(grid_size * sizeof(int));
    grid_next = (int *)malloc(grid_size * sizeof(int));

    if (!grid_current || !grid_next) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize grid with a random a mix of conductors and a few electron heads
    // to ensure interesting simulation activity.
    for (size_t i = 0; i < grid_size; ++i) {
        uint32_t r = mt_rand() % 100;
        if (r < 2) { // 2% chance of being an electron head
            grid_current[i] = 1;
        } else if (r < 50) { // 48% chance of being a conductor
            grid_current[i] = 3;
        } else { // 50% chance of being empty
            grid_current[i] = 0;
        }
    }
}

void run_computation() {
    for (int step = 0; step < num_steps; ++step) {
        for (int y = 0; y < grid_height; ++y) {
            for (int x = 0; x < grid_width; ++x) {
                int current_pos = y * grid_width + x;
                int current_state = grid_current[current_pos];
                int next_state = current_state;

                switch (current_state) {
                    case 0: // Empty
                        break;
                    case 1: // Electron Head
                        next_state = 2; // Becomes Tail
                        break;
                    case 2: // Electron Tail
                        next_state = 3; // Becomes Conductor
                        break;
                    case 3: // Conductor
                        {
                            int head_neighbors = 0;
                            for (int dy = -1; dy <= 1; ++dy) {
                                for (int dx = -1; dx <= 1; ++dx) {
                                    if (dx == 0 && dy == 0) continue;
                                    int nx = x + dx;
                                    int ny = y + dy;
                                    if (nx >= 0 && nx < grid_width && ny >= 0 && ny < grid_height) {
                                        if (grid_current[ny * grid_width + nx] == 1) {
                                            head_neighbors++;
                                        }
                                    }
                                }
                            }
                            if (head_neighbors == 1 || head_neighbors == 2) {
                                next_state = 1; // Becomes Head
                            }
                        }
                        break;
                }
                grid_next[current_pos] = next_state;
            }
        }

        // Swap grids for the next iteration
        int *temp = grid_current;
        grid_current = grid_next;
        grid_next = temp;
    }

    // Calculate a checksum to prevent dead code elimination
    long long checksum = 0;
    size_t grid_size = (size_t)grid_width * grid_height;
    for (size_t i = 0; i < grid_size; ++i) {
        checksum += grid_current[i];
    }
    final_result = checksum;
}

void cleanup() {
    free(grid_current);
    free(grid_next);
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
