/**
 * @file auto_correlation_direct.c
 * @brief Benchmark for direct autocorrelation calculation in DSP.
 * 
 * This program computes the autocorrelation of a signal with itself at
 * different time offsets (lags). Autocorrelation is a mathematical tool used
 * for finding repeating patterns in a signal, such as determining the presence
 * of a periodic signal which is buried under noise.
 * 
 * The direct method calculates the sum of products for each lag, which is
 * computationally intensive (O(N*M) where N is signal length and M is max lag).
 * The benchmark measures the performance of this core computational task.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) BEGIN ---
// (Provided implementation, do not modify)
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
// --- Mersenne Twister (MT19937) END ---

// --- Benchmark Globals ---
int SIGNAL_LENGTH;
int MAX_LAG;

double* signal_data;
double* autocorr_result;
double final_result; // For preventing dead code elimination

/**
 * @brief Generates a random double between -1.0 and 1.0.
 */
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0;
}

/**
 * @brief Sets up benchmark data.
 * 
 * Parses command-line arguments, allocates memory for the signal and result
 * arrays, and populates the signal array with random data using the MT19937
 * generator.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <signal_length> <max_lag> <seed>\n", argv[0]);
        exit(1);
    }

    SIGNAL_LENGTH = atoi(argv[1]);
    MAX_LAG = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (SIGNAL_LENGTH <= 0 || MAX_LAG < 0) {
        fprintf(stderr, "Error: signal_length must be > 0 and max_lag must be >= 0.\n");
        exit(1);
    }
    if (MAX_LAG >= SIGNAL_LENGTH) {
        fprintf(stderr, "Error: max_lag must be less than signal_length.\n");
        exit(1);
    }

    mt_seed(seed);

    signal_data = (double*)malloc(SIGNAL_LENGTH * sizeof(double));
    if (!signal_data) {
        fprintf(stderr, "Failed to allocate memory for signal_data.\n");
        exit(1);
    }

    // MAX_LAG + 1 to include lag 0
    autocorr_result = (double*)malloc((MAX_LAG + 1) * sizeof(double));
    if (!autocorr_result) {
        fprintf(stderr, "Failed to allocate memory for autocorr_result.\n");
        free(signal_data);
        exit(1);
    }

    for (int i = 0; i < SIGNAL_LENGTH; ++i) {
        signal_data[i] = random_double();
    }
}

/**
 * @brief Runs the core computation.
 *
 * Calculates the autocorrelation for each lag from 0 to MAX_LAG.
 * The autocorrelation R(k) for a lag k is given by:
 * R(k) = sum_{n=0}^{N-k-1} x[n] * x[n+k]
 * where x is the signal and N is its length.
 * The final result is a sum of all autocorrelation values to prevent optimization.
 */
void run_computation() {
    for (int k = 0; k <= MAX_LAG; ++k) {
        double sum = 0.0;
        for (int n = 0; n < SIGNAL_LENGTH - k; ++n) {
            sum += signal_data[n] * signal_data[n + k];
        }
        autocorr_result[k] = sum;
    }

    final_result = 0.0;
    for (int k = 0; k <= MAX_LAG; ++k) {
        final_result += autocorr_result[k];
    }
}

/**
 * @brief Frees all allocated memory.
 */
void cleanup() {
    free(signal_data);
    free(autocorr_result);
    signal_data = NULL;
    autocorr_result = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%f\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
