#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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

// --- Benchmark Globals ---
static int matrix_size;
static double* matrix;
static double final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_size> <seed>\n", argv[0]);
        exit(1);
    }

    matrix_size = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (matrix_size <= 0) {
        fprintf(stderr, "FATAL: matrix_size must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    size_t num_elements = (size_t)matrix_size * matrix_size;
    matrix = (double*)malloc(num_elements * sizeof(double));
    if (matrix == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for matrix.\n");
        exit(1);
    }

    for (size_t i = 0; i < num_elements; ++i) {
        matrix[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

void run_computation() {
    // The trace an O(N) operation on an O(N^2) data structure. To create
    // a workload that lasts a meaningful amount of time without requiring
    // enormous amounts of memory, we calculate the trace repeatedly.
    // The total number of additions is kept roughly constant to make runtimes
    // comparable across different matrix sizes.
    const long long total_ops_target = 2000000000LL;
    long repetitions = total_ops_target / matrix_size;
    if (repetitions < 1) {
        repetitions = 1;
    }

    final_result = 0.0;
    for (long k = 0; k < repetitions; ++k) {
        double local_trace = 0.0;
        for (int i = 0; i < matrix_size; ++i) {
            // Accessing diagonal element A[i][i] in a 1D-mapped array
            local_trace += matrix[(size_t)i * matrix_size + i];
        }
        // Accumulate result to prevent compiler from optimizing away the loops
        final_result += local_trace;
    }
}

void cleanup() {
    free(matrix);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
