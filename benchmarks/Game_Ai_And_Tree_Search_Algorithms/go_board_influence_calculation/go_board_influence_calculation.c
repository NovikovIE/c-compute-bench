#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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
// --- END MERSENNE TWISTER ---

// Benchmark parameters
static int board_size;
static int num_stones;

// Data structures
static int* board; // 0=empty, 1=black, 2=white
static float* influence_map;
static int* stone_x;
static int* stone_y;
static int* stone_color;

// Final result accumulator
static long long final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <board_size> <num_stones> <seed>\n", argv[0]);
        exit(1);
    }

    board_size = atoi(argv[1]);
    num_stones = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (board_size <= 0 || num_stones <= 0) {
        fprintf(stderr, "ERROR: board_size and num_stones must be positive integers.\n");
        exit(1);
    }

    if (num_stones > board_size * board_size) {
        fprintf(stderr, "ERROR: num_stones cannot be greater than board_size * board_size.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    unsigned long long board_area = (unsigned long long)board_size * board_size;
    board = (int*)calloc(board_area, sizeof(int));
    influence_map = (float*)malloc(board_area * sizeof(float));
    stone_x = (int*)malloc(num_stones * sizeof(int));
    stone_y = (int*)malloc(num_stones * sizeof(int));
    stone_color = (int*)malloc(num_stones * sizeof(int));

    if (!board || !influence_map || !stone_x || !stone_y || !stone_color) {
        fprintf(stderr, "ERROR: Memory allocation failed.\n");
        exit(1);
    }

    // Place stones randomly on the board, avoiding collisions
    int stones_placed = 0;
    while (stones_placed < num_stones) {
        int x = mt_rand() % board_size;
        int y = mt_rand() % board_size;
        if (board[y * board_size + x] == 0) {
            int color = (mt_rand() % 2) + 1; // 1 for black, 2 for white
            board[y * board_size + x] = color;
            stone_x[stones_placed] = x;
            stone_y[stones_placed] = y;
            stone_color[stones_placed] = color;
            stones_placed++;
        }
    }
}

void run_computation() {
    unsigned long long board_area = (unsigned long long)board_size * board_size;

    // 1. Initialize influence map to zero
    for (unsigned long long i = 0; i < board_area; ++i) {
        influence_map[i] = 0.0f;
    }

    // 2. Calculate influence from each stone
    const float INITIAL_INFLUENCE = 100.0f;
    for (int s = 0; s < num_stones; ++s) {
        int sx = stone_x[s];
        int sy = stone_y[s];
        float stone_influence_sign = (stone_color[s] == 1) ? 1.0f : -1.0f;

        for (int y = 0; y < board_size; ++y) {
            for (int x = 0; x < board_size; ++x) {
                // Manhattan distance
                int dist = abs(x - sx) + abs(y - sy);
                if (dist > 0) {
                    // Influence decays with the square of the distance
                    influence_map[y * board_size + x] +=
                        (stone_influence_sign * INITIAL_INFLUENCE) / (float)(dist * dist);
                }
            }
        }
    }

    // 3. Accumulate a final result to prevent dead code elimination
    long long total_influence = 0;
    for (unsigned long long i = 0; i < board_area; ++i) {
        // Scale and cast to long long for accumulation
        total_influence += (long long)(influence_map[i] * 100.0f);
    }
    final_result = total_influence;
}

void cleanup() {
    free(board);
    free(influence_map);
    free(stone_x);
    free(stone_y);
    free(stone_color);
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
