/**
 * @file convolution_direct.c
 * @brief Benchmark for direct 1D convolution, a fundamental DSP operation.
 *
 * This program measures the performance of calculating the convolution of a large
 * input signal with a smaller kernel using the direct summation method. Convolution
 * is a mathematical operation on two functions (or signals) that produces a third
 * function expressing how the shape of one is modified by the other.
 * The formula for discrete convolution is: y[n] = sum_{k=-inf}^{inf} x[k] * h[n-k]
 * For finite signals, this becomes y[n] = sum_{k=0}^{M-1} h[k] * x[n-k].
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) a_prng --- (DO NOT MODIFY)
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
// --- End of Mersenne Twister ---

// Benchmark parameters
long input_signal_length;
long kernel_length;
long output_signal_length;

// Data arrays
double *input_signal;
double *kernel;
double *output_signal;

// Final result to prevent dead code elimination
double final_checksum = 0.0;

/**
 * @brief Sets up the benchmark by parsing arguments, allocating memory, and
 *        generating pseudo-random input data.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_signal_length> <kernel_length> <seed>\n", argv[0]);
        exit(1);
    }

    input_signal_length = atol(argv[1]);
    kernel_length = atol(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (input_signal_length <= 0 || kernel_length <= 0) {
        fprintf(stderr, "FATAL: Signal and kernel lengths must be positive.\n");
        exit(1);
    }
    if (kernel_length > input_signal_length) {
        fprintf(stderr, "FATAL: Kernel length cannot be greater than signal length.\n");
        exit(1);
    }

    mt_seed(seed);

    output_signal_length = input_signal_length + kernel_length - 1;

    // Allocate memory
    input_signal = (double*)malloc(input_signal_length * sizeof(double));
    kernel = (double*)malloc(kernel_length * sizeof(double));
    output_signal = (double*)malloc(output_signal_length * sizeof(double));

    if (!input_signal || !kernel || !output_signal) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize input signal and kernel with random values between -1.0 and 1.0
    for (long i = 0; i < input_signal_length; ++i) {
        input_signal[i] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }
    for (long i = 0; i < kernel_length; ++i) {
        kernel[i] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }
}

/**
 * @brief Runs the core computation of the benchmark: direct 1D convolution.
 */
void run_computation() {
    // Perform direct convolution
    for (long n = 0; n < output_signal_length; ++n) {
        double sum = 0.0;
        for (long k = 0; k < kernel_length; ++k) {
            long input_index = n - k;
            // Check boundaries: equivalent to zero-padding the input signal
            if (input_index >= 0 && input_index < input_signal_length) {
                sum += input_signal[input_index] * kernel[k];
            }
        }
        output_signal[n] = sum;
    }

    // Calculate a checksum to ensure computation is not optimized away
    double checksum = 0.0;
    for (long i = 0; i < output_signal_length; ++i) {
        checksum += output_signal[i];
    }
    final_checksum = checksum;
}

/**
 * @brief Frees all memory allocated during setup.
 */
void cleanup() {
    free(input_signal);
    free(kernel);
    free(output_signal);
}

/**
 * @brief Main function to orchestrate the benchmark.
 */
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final checksum to stdout
    printf("%f\n", final_checksum);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
