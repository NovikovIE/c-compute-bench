#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

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
int N; // matrix_size
int P; // power
double *A_base;    // The matrix to be powered
double *R_result;  // The final result matrix (holds A^P)
double *M_power;   // A^(2^k) in the exponentiation by squaring
double *T_temp;    // Temporary for matrix multiplication products
double final_sum;  // Sum of elements in R_result for output

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_size> <power> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    P = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (N <= 0 || P < 0) {
        fprintf(stderr, "ERROR: matrix_size must be positive and power must be non-negative.\n");
        exit(1);
    }
    
    mt_seed(seed);

    size_t matrix_bytes = (size_t)N * N * sizeof(double);
    A_base   = (double*)malloc(matrix_bytes);
    R_result = (double*)malloc(matrix_bytes);
    M_power  = (double*)malloc(matrix_bytes);
    T_temp   = (double*)malloc(matrix_bytes);

    if (!A_base || !R_result || !M_power || !T_temp) {
        fprintf(stderr, "ERROR: Memory allocation failed.\n");
        exit(1);
    }
    
    for (int i = 0; i < N * N; ++i) {
        // Generate a small random double between 0 and 1 for numerical stability
        A_base[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

void cleanup() {
    free(A_base);
    free(R_result);
    free(M_power);
    free(T_temp);
}

// Helper to multiply dense matrices C = A * B
// All matrices are n x n and stored in row-major order
void multiply(double *C, const double *A, const double *B, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
}

void run_computation() {
    size_t matrix_bytes = (size_t)N * N * sizeof(double);

    if (P == 0) {
        // If power is 0, result is the identity matrix
        memset(R_result, 0, matrix_bytes);
        for(int i = 0; i < N; ++i) {
            R_result[i * N + i] = 1.0;
        }
    } else {
        // Initialize R_result to identity matrix
        memset(R_result, 0, matrix_bytes);
        for(int i = 0; i < N; ++i) {
            R_result[i * N + i] = 1.0;
        }

        // Initialize M_power to A_base
        memcpy(M_power, A_base, matrix_bytes);

        int p_local = P;
        // Exponentiation by squaring
        while (p_local > 0) {
            // If p is odd, R_result = R_result * M_power
            if (p_local & 1) {
                multiply(T_temp, R_result, M_power, N);
                memcpy(R_result, T_temp, matrix_bytes);
            }
            // M_power = M_power * M_power
            multiply(T_temp, M_power, M_power, N);
            memcpy(M_power, T_temp, matrix_bytes);
            
            p_local >>= 1; // p = p / 2
        }
    }

    // Calculate sum of elements to prevent dead code elimination
    final_sum = 0.0;
    for (int i = 0; i < N * N; ++i) {
        final_sum += R_result[i];
    }
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
    printf("%f\n", final_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
