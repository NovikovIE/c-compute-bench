#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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


// --- Benchmark Globals ---
int num_dimensions;
int num_iterations;
double *x; // Solution vector
double *Q; // Quadratic form matrix (symmetric, positive-definite)
double *b; // Linear form vector
double final_result;

// --- Benchmark Functions ---

// Generates a random double between -1.0 and 1.0
double rand_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_dimensions> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    num_dimensions = atoi(argv[1]);
    num_iterations = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    
    mt_seed(seed);

    if (num_dimensions <= 0 || num_iterations <= 0) {
        fprintf(stderr, "FATAL: num_dimensions and num_iterations must be positive.\n");
        exit(1);
    }
    
    // Allocate memory
    x = (double*)malloc(num_dimensions * sizeof(double));
    b = (double*)malloc(num_dimensions * sizeof(double));
    Q = (double*)malloc((size_t)num_dimensions * num_dimensions * sizeof(double));

    if (!x || !b || !Q) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize data
    // Initialize solution vector x to zeros
    for (int i = 0; i < num_dimensions; i++) {
        x[i] = 0.0;
    }

    // Initialize vector b with random values
    for (int i = 0; i < num_dimensions; i++) {
        b[i] = rand_double();
    }

    // Initialize matrix Q to be symmetric and diagonally dominant.
    // This ensures there's a unique minimum and the algorithm is stable.
    for (int i = 0; i < num_dimensions; i++) {
        for (int j = i; j < num_dimensions; j++) {
            if (i == j) {
                // Large positive diagonal elements
                Q[i * num_dimensions + j] = (double)num_dimensions;
            } else {
                double val = rand_double();
                Q[i * num_dimensions + j] = val;
                Q[j * num_dimensions + i] = val; // Symmetry
            }
        }
    }
}

void run_computation() {
    // Coordinate descent for minimizing f(x) = 0.5 * x'Qx - b'x
    // The update rule for a single coordinate x_i is:
    // x_i = (b_i - sum_{j!=i} Q_{ij} * x_j) / Q_{ii}
    for (int iter = 0; iter < num_iterations; iter++) {
        for (int i = 0; i < num_dimensions; i++) {
            double sigma = 0.0;
            for (int j = 0; j < num_dimensions; j++) {
                if (i != j) {
                    sigma += Q[i * num_dimensions + j] * x[j];
                }
            }
            x[i] = (b[i] - sigma) / Q[i * num_dimensions + i];
        }
    }
    
    // Accumulate a result to prevent dead code elimination
    double sum = 0.0;
    for (int i = 0; i < num_dimensions; i++) {
        sum += x[i];
    }
    final_result = sum;
}

void cleanup() {
    free(x);
    free(b);
    free(Q);
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

    // Print final result to stdout
    printf("%f\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
