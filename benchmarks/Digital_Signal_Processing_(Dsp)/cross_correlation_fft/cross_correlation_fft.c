#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- START: Mersenne Twister (MT19937) implementation ---
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
// --- END: Mersenne Twister (MT19937) implementation ---

// Benchmark parameters
static int signal_a_length;
static int signal_b_length;
static int fft_size;

// Data arrays
static double complex *signal_a_fft; // Holds padded signal A, then its FFT, then correlation product
static double complex *signal_b_fft; // Holds padded signal B, then its FFT

static double final_result;

// Helper to find the next power of 2, required for FFT
static int next_power_of_2(int n) {
    int power = 1;
    while (power < n) {
        power <<= 1;
    }
    return power;
}

// In-place Cooley-Tukey Radix-2 FFT/IFFT
// is_inverse = -1 for forward FFT, +1 for inverse FFT
static void fft_inplace(double complex *buf, int n, int is_inverse) {
    if (n <= 1) return;

    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            double complex temp = buf[i];
            buf[i] = buf[j];
            buf[j] = temp;
        }
    }

    // Iterative FFT
    for (int len = 2; len <= n; len <<= 1) {
        double angle = is_inverse * 2 * M_PI / len;
        double complex wlen = cexp(I * angle);
        for (int i = 0; i < n; i += len) {
            double complex w = 1.0;
            for (int j = 0; j < len / 2; j++) {
                double complex u = buf[i + j];
                double complex v = buf[i + j + len / 2] * w;
                buf[i + j] = u + v;
                buf[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s signal_a_length signal_b_length seed\n", argv[0]);
        exit(1);
    }

    signal_a_length = atoi(argv[1]);
    signal_b_length = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // For linear cross-correlation, pad to length N1 + N2 - 1
    int min_fft_size = signal_a_length + signal_b_length - 1;
    fft_size = next_power_of_2(min_fft_size);

    // Allocate memory for the signals
    signal_a_fft = (double complex *)calloc(fft_size, sizeof(double complex));
    signal_b_fft = (double complex *)calloc(fft_size, sizeof(double complex));

    if (!signal_a_fft || !signal_b_fft) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate random signal A and pad with zeros
    for (int i = 0; i < signal_a_length; ++i) {
        signal_a_fft[i] = (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0;
    }

    // Generate random signal B and pad with zeros
    for (int i = 0; i < signal_b_length; ++i) {
        signal_b_fft[i] = (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0;
    }
}

void run_computation() {
    // 1. Compute FFT of signal A
    fft_inplace(signal_a_fft, fft_size, -1);

    // 2. Compute FFT of signal B
    fft_inplace(signal_b_fft, fft_size, 1); // Use +1 for conjugate trick

    // 3. Point-wise product of FFT(A) and conj(FFT(B))
    // We computed FFT(B) with +1 (IFFT), which is equivalent to conj(FFT(B))
    // up to a scaling factor, so we combine it here.
    for (int i = 0; i < fft_size; ++i) {
        signal_a_fft[i] *= signal_b_fft[i];
    }

    // 4. Compute Inverse FFT of the product
    fft_inplace(signal_a_fft, fft_size, 1);

    // 5. Scale and accumulate result
    double sum = 0.0;
    double scale = (double)fft_size;
    for (int i = 0; i < fft_size; ++i) {
        // IFFT of FFT(A)*IFFT(B) = IFFT(A) * B. We need IFFT(FFT(A)*conj(FFT(B))).
        // The sign trick for FFT(B) above gave conj(FFT(B))/N. So overall result needs *N scaling.
        // And the IFFT added another 1/N scaling. So total scaling is neutral.
        // Correct scaling for IFFT of FFT()*conj(FFT()) is 1/N for each transform.
        // Our IFFT doesn't scale. We scale at the end.
        double val = creal(signal_a_fft[i]) / scale;
        sum += val;
    }
    final_result = sum;
}

void cleanup() {
    free(signal_a_fft);
    free(signal_b_fft);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}