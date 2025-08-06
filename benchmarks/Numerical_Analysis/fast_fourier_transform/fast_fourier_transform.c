#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Define PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

// --- Benchmark Globals ---
typedef struct {
    double real;
    double imag;
} Complex;

int ARRAY_SIZE;
int DIMENSIONS;
unsigned long long TOTAL_ELEMENTS;
Complex *data;
double final_result_sum;

// --- Helper Functions ---
int is_power_of_two(int n) {
    if (n <= 0) return 0;
    return (n & (n - 1)) == 0;
}

// In-place 1D Cooley-Tukey FFT
void fft_1d(Complex *buf, int n) {
    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            Complex temp = buf[i];
            buf[i] = buf[j];
            buf[j] = temp;
        }
    }

    // Cooley-Tukey passes
    for (int len = 2; len <= n; len <<= 1) {
        double angle = -2.0 * M_PI / len;
        Complex wlen = {cos(angle), sin(angle)};
        for (int i = 0; i < n; i += len) {
            Complex w = {1.0, 0.0};
            for (int j = 0; j < len / 2; j++) {
                Complex u = buf[i + j];
                Complex v;
                v.real = buf[i + j + len / 2].real * w.real - buf[i + j + len / 2].imag * w.imag;
                v.imag = buf[i + j + len / 2].real * w.imag + buf[i + j + len / 2].imag * w.real;
                
                buf[i + j].real = u.real + v.real;
                buf[i + j].imag = u.imag + v.imag;

                buf[i + j + len / 2].real = u.real - v.real;
                buf[i + j + len / 2].imag = u.imag - v.imag;
                
                double w_temp_real = w.real * wlen.real - w.imag * wlen.imag;
                w.imag = w.real * wlen.imag + w.imag * wlen.real;
                w.real = w_temp_real;
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

    ARRAY_SIZE = atoi(argv[1]);
    DIMENSIONS = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (!is_power_of_two(ARRAY_SIZE)) {
        fprintf(stderr, "Error: array_size must be a power of two.\n");
        exit(1);
    }
    if (DIMENSIONS < 1 || DIMENSIONS > 2) {
        fprintf(stderr, "Error: This implementation supports 1 or 2 dimensions only.\n");
        exit(1);
    }

    TOTAL_ELEMENTS = 1;
    for(int i = 0; i < DIMENSIONS; ++i) {
        TOTAL_ELEMENTS *= ARRAY_SIZE;
    }

    data = (Complex *)malloc(TOTAL_ELEMENTS * sizeof(Complex));
    if (data == NULL) {
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        exit(1);
    }

    mt_seed(seed);
    for (unsigned long long i = 0; i < TOTAL_ELEMENTS; i++) {
        data[i].real = (double)mt_rand() / (double)UINT32_MAX;
        data[i].imag = (double)mt_rand() / (double)UINT32_MAX;
    }

    final_result_sum = 0.0;
}

void run_computation() {
    if (DIMENSIONS == 1) {
        fft_1d(data, ARRAY_SIZE);
    } else if (DIMENSIONS == 2) {
        Complex *temp_col = (Complex *)malloc(ARRAY_SIZE * sizeof(Complex));
        if (temp_col == NULL) {
            fprintf(stderr, "Failed to allocate temp column buffer.\n");
            exit(1);
        }

        // Apply 1D FFT to each row
        for (int i = 0; i < ARRAY_SIZE; i++) {
            fft_1d(data + i * ARRAY_SIZE, ARRAY_SIZE);
        }

        // Apply 1D FFT to each column
        for (int j = 0; j < ARRAY_SIZE; j++) {
            // Extract column j to temp buffer
            for (int i = 0; i < ARRAY_SIZE; i++) {
                temp_col[i] = data[i * ARRAY_SIZE + j];
            }
            // Transform the column
            fft_1d(temp_col, ARRAY_SIZE);
            // Put transformed column back
            for (int i = 0; i < ARRAY_SIZE; i++) {
                data[i * ARRAY_SIZE + j] = temp_col[i];
            }
        }
        free(temp_col);
    }

    // Calculate a checksum to prevent dead code elimination
    for (unsigned long long i = 0; i < TOTAL_ELEMENTS; i++) {
        final_result_sum += sqrt(data[i].real * data[i].real + data[i].imag * data[i].imag);
    }
}

void cleanup() {
    free(data);
    data = NULL;
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
    printf("%f\n", final_result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
