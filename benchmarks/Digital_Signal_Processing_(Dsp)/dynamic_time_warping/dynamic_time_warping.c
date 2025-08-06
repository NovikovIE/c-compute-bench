#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Globals and Helpers ---

// Benchmark parameters
int N; // sequence1_length
int M; // sequence2_length
int W; // window_size

// Data pointers
double* seq1;
double* seq2;
double* dtw_row1;
double* dtw_row2;

// Result
double final_cost;

// Utility macros
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

// Distance function (squared Euclidean)
static inline double dist(double a, double b) {
    double d = a - b;
    return d * d;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <sequence1_length> <sequence2_length> <window_size> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    M = atoi(argv[2]);
    W = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    if (N <= 0 || M <= 0 || W < 0) {
        fprintf(stderr, "FATAL: lengths must be > 0 and window_size must be >= 0.\n");
        exit(1);
    }

    if (abs(N - M) > W) {
       // A valid path may not be found; proceed, but result might be inf/-1.
    }

    mt_seed(seed);

    seq1 = (double*)malloc(N * sizeof(double));
    seq2 = (double*)malloc(M * sizeof(double));
    dtw_row1 = (double*)malloc(M * sizeof(double));
    dtw_row2 = (double*)malloc(M * sizeof(double));

    if (!seq1 || !seq2 || !dtw_row1 || !dtw_row2) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < N; i++) {
        seq1[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
    for (int i = 0; i < M; i++) {
        seq2[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

void run_computation() {
    // Two rows of the cost matrix are stored to save space.
    // prev_row is cost[i-1], curr_row is cost[i].
    double *prev_row = dtw_row1;
    double *curr_row = dtw_row2;

    // --- Initialize first row (i=0) based on calculation in main loop ---
    // Effectively, prev_row represents the state before the loop starts (row i=-1).
    // This initialization makes the main loop cleaner.
    for (int j = 0; j < M; j++) {
        prev_row[j] = DBL_MAX;
    }
    prev_row[0] = 0; // The path must start from a zero-cost origin before (0,0)

    // --- Main DP loop ---
    for (int i = 0; i < N; i++) {
        int j_start = MAX(0, i - W);
        int j_end = MIN(M - 1, i + W);

        // Initialize first element of curr_row separately
        if (j_start == 0) {
            curr_row[0] = prev_row[0] + dist(seq1[i], seq2[0]);
        } else {
          curr_row[0] = DBL_MAX;
        }

        for (int j = MAX(1, j_start); j <= j_end; j++) {
            double min_prev_cost = fmin(curr_row[j-1], fmin(prev_row[j], prev_row[j-1]));
            if (min_prev_cost == DBL_MAX) {
                curr_row[j] = DBL_MAX;
            } else {
                curr_row[j] = min_prev_cost + dist(seq1[i], seq2[j]);
            }
        }

        // Swap rows for the next iteration
        double *temp = prev_row;
        prev_row = curr_row;
        curr_row = temp;
    }

    final_cost = prev_row[M - 1];
    if (final_cost == DBL_MAX) {
        final_cost = -1.0; // Indicate path not found
    }
}

void cleanup() {
    free(seq1);
    free(seq2);
    free(dtw_row1);
    free(dtw_row2);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lf\n", final_cost);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
