#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <inttypes.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

#define NUM_PIECE_TYPES 2

// --- Global Benchmark Data ---
typedef struct {
    int x;
    int y;
    int piece_type;
} Move;

static int g_board_size;
static long g_num_positions;

static uint64_t* zobrist_table = NULL;
static Move* moves = NULL;
static uint64_t final_hash;

// Helper to generate a 64-bit random number from the MT
static inline uint64_t mt_rand64(void) {
    uint64_t r1 = mt_rand();
    uint64_t r2 = mt_rand();
    return (r1 << 32) | r2;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <board_size> <num_positions> <seed>\n", argv[0]);
        exit(1);
    }

    g_board_size = atoi(argv[1]);
    g_num_positions = atol(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_board_size <= 0 || g_num_positions <= 0) {
        fprintf(stderr, "Invalid arguments: board_size and num_positions must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate Zobrist table
    size_t table_size = (size_t)g_board_size * g_board_size * NUM_PIECE_TYPES;
    zobrist_table = (uint64_t*)malloc(table_size * sizeof(uint64_t));
    if (!zobrist_table) {
        fprintf(stderr, "Failed to allocate Zobrist table.\n");
        exit(1);
    }

    // Fill Zobrist table with random 64-bit numbers
    for (size_t i = 0; i < table_size; ++i) {
        zobrist_table[i] = mt_rand64();
    }

    // Allocate moves array
    moves = (Move*)malloc((size_t)g_num_positions * sizeof(Move));
    if (!moves) {
        fprintf(stderr, "Failed to allocate moves array.\n");
        free(zobrist_table);
        exit(1);
    }

    // Generate random moves
    for (long i = 0; i < g_num_positions; ++i) {
        moves[i].x = mt_rand() % g_board_size;
        moves[i].y = mt_rand() % g_board_size;
        moves[i].piece_type = mt_rand() % NUM_PIECE_TYPES;
    }
}

void run_computation() {
    uint64_t current_hash = 0;
    for (long i = 0; i < g_num_positions; ++i) {
        Move m = moves[i];
        size_t index = ((size_t)m.y * g_board_size + m.x) * NUM_PIECE_TYPES + m.piece_type;
        current_hash ^= zobrist_table[index];
    }
    final_hash = current_hash;
}

void cleanup() {
    free(zobrist_table);
    free(moves);
}


// --- Main ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%" PRIu64 "\n", final_hash);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
