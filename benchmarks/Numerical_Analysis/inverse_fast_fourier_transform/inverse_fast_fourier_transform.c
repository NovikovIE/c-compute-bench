#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
// --- End Mersenne Twister ---

// --- Benchmark Data and Globals ---
typedef struct {
    int array_size;
    int dimensions;
    unsigned long long total_elements;
    double complex* data;
    double result_checksum;
} BenchmarkData;

static BenchmarkData B;

// --- Helper Functions ---

// Performs a 1D iterative Cooley-Tukey IFFT. The result is unscaled.
static void ifft_1d(double complex* x, int n) {
    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            double complex temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }

    // Cooley-Tukey iterative algorithm
    for (int len = 2; len <= n; len <<= 1) {
        double angle = 2.0 * M_PI / len;
        double complex wlen = cexp(I * angle); // Positive sign for IFFT
        for (int i = 0; i < n; i += len) {
            double complex w = 1.0;
            for (int j = 0; j < len / 2; j++) {
                double complex u = x[i + j];
                double complex v = x[i + j + len / 2] * w;
                x[i + j] = u + v;
                x[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s array_size dimensions seed\n", argv[0]);
        exit(1);
    }

    B.array_size = atoi(argv[1]);
    B.dimensions = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (B.array_size <= 0 || (B.array_size & (B.array_size - 1)) != 0) {
        fprintf(stderr, "FATAL: array_size must be a positive power of 2.\n");
        exit(1);
    }
    if (B.dimensions <= 0) {
        fprintf(stderr, "FATAL: dimensions must be positive.\n");
        exit(1);
    }

    B.total_elements = 1;
    for(int i = 0; i < B.dimensions; ++i) {
        unsigned long long prev_total = B.total_elements;
        B.total_elements *= B.array_size;
        if (B.total_elements / B.array_size != prev_total) { // Overflow check
            fprintf(stderr, "FATAL: Total element count overflows unsigned long long.\n");
            exit(1);
        }
    }

    B.data = (double complex*)malloc(B.total_elements * sizeof(double complex));
    if (B.data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for data array.\n");
        exit(1);
    }

    mt_seed(seed);
    for (unsigned long long i = 0; i < B.total_elements; i++) {
        double real_part = (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0;
        double imag_part = (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0;
        B.data[i] = real_part + imag_part * I;
    }

    B.result_checksum = 0.0;
}

void run_computation() {
    double complex* temp_slice = (double complex*)malloc(B.array_size * sizeof(double complex));
    if (!temp_slice) {
        fprintf(stderr, "FATAL: Memory allocation failed for temp slice.\n");
        exit(1);
    }

    long long N = B.array_size;

    // Perform 1D IFFT along each dimension
    for (int d = 0; d < B.dimensions; d++) {
        long long stride = 1;
        for (int k = 0; k < d; k++) stride *= N;
        long long num_chunks = B.total_elements / (stride * N);

        for (long long chunk = 0; chunk < num_chunks; chunk++) {
            for (long long offset = 0; offset < stride; offset++) {
                long long base_idx = chunk * stride * N + offset;

                // Gather the 1D slice
                for (long long k = 0; k < N; k++) {
                    temp_slice[k] = B.data[base_idx + k * stride];
                }

                // In-place 1D IFFT on the slice
                ifft_1d(temp_slice, N);

                // Scatter the result back
                for (long long k = 0; k < N; k++) {
                    B.data[base_idx + k * stride] = temp_slice[k];
                }
            }
        }
    } 

    free(temp_slice);

    // Final scaling and checksum calculation
    double scale_factor = 1.0 / (double)B.total_elements;
    double checksum = 0.0;
    for (unsigned long long i = 0; i < B.total_elements; i++) {
        B.data[i] *= scale_factor;
        checksum += cabs(B.data[i]);
    }
    B.result_checksum = checksum;
}

void cleanup() {
    if (B.data) {
        free(B.data);
        B.data = NULL;
    }
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
    printf("%f\n", B.result_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
