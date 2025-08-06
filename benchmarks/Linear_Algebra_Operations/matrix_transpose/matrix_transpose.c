#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Globals ---
typedef double DTYPE;

int MATRIX_ROWS, MATRIX_COLS;
DTYPE *matrix_A;
DTYPE *matrix_B; // Transposed matrix
double verification_sum = 0.0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_rows> <matrix_cols> <seed>\n", argv[0]);
        exit(1);
    }

    MATRIX_ROWS = atoi(argv[1]);
    MATRIX_COLS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (MATRIX_ROWS <= 0 || MATRIX_COLS <= 0) {
        fprintf(stderr, "Error: Matrix dimensions must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    matrix_A = (DTYPE *)malloc((size_t)MATRIX_ROWS * MATRIX_COLS * sizeof(DTYPE));
    matrix_B = (DTYPE *)malloc((size_t)MATRIX_COLS * MATRIX_ROWS * sizeof(DTYPE));

    if (matrix_A == NULL || matrix_B == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < MATRIX_ROWS; ++i) {
        for (int j = 0; j < MATRIX_COLS; ++j) {
            matrix_A[i * MATRIX_COLS + j] = (DTYPE)mt_rand() / (DTYPE)UINT32_MAX;
        }
    }
}

void run_computation() {
    // Transpose matrix_A into matrix_B
    for (int i = 0; i < MATRIX_ROWS; ++i) {
        for (int j = 0; j < MATRIX_COLS; ++j) {
            matrix_B[j * MATRIX_ROWS + i] = matrix_A[i * MATRIX_COLS + j];
        }
    }

    // Calculate a sum of the transposed matrix to prevent dead code elimination.
    // This is a simple checksum that ensures the computation has an observable effect.
    double sum = 0.0;
    for (int i = 0; i < MATRIX_COLS * MATRIX_ROWS; ++i) {
        sum += matrix_B[i];
    }
    verification_sum = sum;
}

void cleanup() {
    free(matrix_A);
    free(matrix_B);
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
    printf("%f\n", verification_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
