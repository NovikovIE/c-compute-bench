#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// Mersenne Twister (MT19937) PRNG
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

// --- BENCHMARK DATA AND PARAMETERS ---

int N; // Matrix size
double** A; // Matrix
double final_result; // Accumulated result to prevent dead-code elimination

// --- UTILITY FUNCTIONS ---

double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <matrix_size> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);

    if (N <= 0) {
        fprintf(stderr, "ERROR: matrix_size must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    A = (double**)malloc(N * sizeof(double*));
    if (A == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for matrix rows.\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        A[i] = (double*)malloc(N * sizeof(double));
        if (A[i] == NULL) {
            fprintf(stderr, "FATAL: Memory allocation failed for matrix columns.\n");
            exit(1);
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = random_double() * 200.0 - 100.0; // Random values between -100 and 100
        }
    }

    final_result = 0.0;
}

void run_computation() {
    for (int j = 0; j < N; j++) {
        for (int i = j + 1; i < N; i++) {
            double a = A[j][j];
            double b = A[i][j];

            if (fabs(b) > 1e-9) { // No rotation needed if element is already zero
                double r, c, s;

                r = sqrt(a * a + b * b);
                c = a / r;
                s = -b / r;

                // Apply the rotation to the remaining columns of rows j and i
                for (int k = j; k < N; k++) {
                    double ajk = A[j][k];
                    double aik = A[i][k];
                    A[j][k] = c * ajk - s * aik;
                    A[i][k] = s * ajk + c * aik;
                }
            }
        }
    }

    // Accumulate the sum of diagonal elements to produce a result
    // This prevents the compiler from optimizing away the entire computation
    for (int i = 0; i < N; i++) {
        final_result += A[i][i];
    }
}

void cleanup() {
    if (A) {
        for (int i = 0; i < N; i++) {
            if (A[i]) {
                free(A[i]);
            }
        }
        free(A);
    }
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
