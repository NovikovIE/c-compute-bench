#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (Do Not Modify) ---
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

// --- Benchmark Globals ---
typedef struct {
    int row_length;
    int num_rows;
    char* current_row;
    char* next_row;
    unsigned long long total_ones;
} BenchmarkData;

BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <row_length> <num_rows> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.row_length = atoi(argv[1]);
    g_data.num_rows = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.row_length <= 1 || g_data.num_rows <= 0) {
        fprintf(stderr, "FATAL: row_length must be > 1 and num_rows > 0.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.current_row = (char*)malloc(g_data.row_length * sizeof(char));
    g_data.next_row = (char*)malloc(g_data.row_length * sizeof(char));
    g_data.total_ones = 0;

    if (!g_data.current_row || !g_data.next_row) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize the first row randomly and count its 'ones'.
    for (int i = 0; i < g_data.row_length; i++) {
        char cell_state = mt_rand() % 2;
        g_data.current_row[i] = cell_state;
        g_data.total_ones += cell_state;
    }
}

void run_computation() {
    for (int r = 1; r < g_data.num_rows; r++) {
        // Handle boundaries by assuming cells outside the row are 0.
        // Cell 0:
        char left0 = 0;
        char center0 = g_data.current_row[0];
        char right0 = g_data.current_row[1];
        char next_state0 = left0 ^ (center0 | right0);
        g_data.next_row[0] = next_state0;
        g_data.total_ones += next_state0;

        // Internal cells:
        for (int i = 1; i < g_data.row_length - 1; i++) {
            char left = g_data.current_row[i - 1];
            char center = g_data.current_row[i];
            char right = g_data.current_row[i + 1];
            char next_state = left ^ (center | right);
            g_data.next_row[i] = next_state;
            g_data.total_ones += next_state;
        }

        // Cell ROW_LENGTH - 1:
        char leftN = g_data.current_row[g_data.row_length - 2];
        char centerN = g_data.current_row[g_data.row_length - 1];
        char rightN = 0;
        char next_stateN = leftN ^ (centerN | rightN);
        g_data.next_row[g_data.row_length - 1] = next_stateN;
        g_data.total_ones += next_stateN;

        // Swap buffers for the next iteration.
        char* temp = g_data.current_row;
        g_data.current_row = g_data.next_row;
        g_data.next_row = temp;
    }
}

void cleanup() {
    free(g_data.current_row);
    free(g_data.next_row);
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
    printf("%llu\n", g_data.total_ones);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
