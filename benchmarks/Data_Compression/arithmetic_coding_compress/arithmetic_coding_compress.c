#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// Mersenne Twister (verbatim)
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

// Benchmark-specific global variables
static uint8_t *input_data;
static size_t input_size;
static int model_precision_bits;

// Frequency model for arithmetic coding
static uint64_t *freqs;
static uint64_t *cum_freqs;

// Result of computation (must be checksum of some sort)
static uint64_t final_result;

// Arithmetic coder state constants
static uint64_t top_value;
static uint64_t half;
static uint64_t quarter;
static uint64_t three_quarters;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <model_precision_bits> <seed>\n", argv[0]);
        exit(1);
    }

    size_t input_size_mb = atoi(argv[1]);
    model_precision_bits = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);
    
    input_size = input_size_mb * 1024 * 1024;

    if (model_precision_bits < 16 || model_precision_bits > 62) {
        fprintf(stderr, "Model precision bits must be between 16 and 62.\n");
        exit(1);
    }
    
    // Allocate memory
    input_data = (uint8_t *)malloc(input_size * sizeof(uint8_t));
    freqs = (uint64_t *)calloc(256, sizeof(uint64_t));
    cum_freqs = (uint64_t *)calloc(257, sizeof(uint64_t)); // 257 for cumulative property

    if (!input_data || !freqs || !cum_freqs) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate pseudo-compressible random data by creating a non-uniform distribution
    for (size_t i = 0; i < input_size; i++) {
        // Products of small random numbers create a skewed distribution where smaller values are more frequent
        uint8_t val1 = mt_rand() % 16;
        uint8_t val2 = mt_rand() % 16;
        input_data[i] = val1 * val2;
    }

    // Build frequency model from input data
    for (size_t i = 0; i < input_size; i++) {
        freqs[input_data[i]]++;
    }

    // Build cumulative frequency table. Give every possible symbol a minimum count of 1
    // to avoid zero-width intervals, which would break the arithmetic coder.
    cum_freqs[0] = 0;
    for (int i = 0; i < 256; i++) {
        if (freqs[i] == 0) {
            freqs[i] = 1;
        }
        cum_freqs[i + 1] = cum_freqs[i] + freqs[i];
    }

    // Set up constants for the arithmetic coder based on precision
    top_value = (1ULL << model_precision_bits) - 1;
    half = (top_value >> 1) + 1;
    quarter = (top_value >> 2) + 1;
    three_quarters = quarter * 3;
}

void run_computation() {
    uint64_t low = 0;
    uint64_t high = top_value;
    uint64_t underflow_bits = 0;
    const uint64_t total_freqs = cum_freqs[256];

    for (size_t i = 0; i < input_size; ++i) {
        uint8_t symbol = input_data[i];
        
        // 1. Update the range based on the symbol's frequency
        const uint64_t range = high - low + 1;
        const uint64_t symbol_low_freq = cum_freqs[symbol];
        const uint64_t symbol_high_freq = cum_freqs[symbol + 1];

        high = low + (range * symbol_high_freq / total_freqs) - 1;
        low = low + (range * symbol_low_freq / total_freqs);

        // 2. Renormalize the range (E1, E2, E3 conditions) to keep precision
        for (;;) {
            if (high < half) {
                // E1 condition: range is in the lower half. Shift out MSB 0.
                low <<= 1;
                high = (high << 1) | 1;
                underflow_bits = 0; // In a real coder, output 0 + pending 1s
            } else if (low >= half) {
                // E2 condition: range is in the upper half. Shift out MSB 1.
                low = (low - half) << 1;
                high = (high - half) << 1 | 1;
                 underflow_bits = 0; // In a real coder, output 1 + pending 0s
            } else if (low >= quarter && high < three_quarters) {
                // E3 condition: range contains the midpoint. Underflow is possible.
                low = (low - quarter) << 1;
                high = (high - quarter) << 1 | 1;
                underflow_bits++;
            } else {
                // Range is large enough, no more renormalization needed for this symbol.
                break;
            }
        }
    }
    
    // Use the final state to produce a result, preventing dead code elimination.
    final_result = low + underflow_bits;
}

void cleanup() {
    free(input_data);
    input_data = NULL;
    free(freqs);
    freqs = NULL;
    free(cum_freqs);
    cum_freqs = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final computational result to stdout
    printf("%llu\n", (unsigned long long)final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
