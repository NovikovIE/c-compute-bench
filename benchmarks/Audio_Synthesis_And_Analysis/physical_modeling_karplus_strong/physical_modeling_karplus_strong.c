#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND PARAMETERS ---
typedef struct {
    long duration_seconds;
    long sample_rate_hz;
    long delay_line_length;
    long num_samples;
    float* output_buffer;
    float* delay_line;
    double result_accumulator;
} BenchmarkData;

BenchmarkData g_data;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s duration_seconds sample_rate_hz delay_line_length seed\n", argv[0]);
        exit(1);
    }

    g_data.duration_seconds = atol(argv[1]);
    g_data.sample_rate_hz = atol(argv[2]);
    g_data.delay_line_length = atol(argv[3]);
    uint32_t seed = (uint32_t)atol(argv[4]);
    mt_seed(seed);

    if (g_data.duration_seconds <= 0 || g_data.sample_rate_hz <= 0 || g_data.delay_line_length <= 1) {
        fprintf(stderr, "FATAL: Invalid parameters. Values must be positive and delay line length > 1.\n");
        exit(1);
    }

    g_data.num_samples = g_data.duration_seconds * g_data.sample_rate_hz;
    
    g_data.output_buffer = (float*)malloc(g_data.num_samples * sizeof(float));
    if (!g_data.output_buffer) {
        fprintf(stderr, "FATAL: Failed to allocate memory for output buffer.\n");
        exit(1);
    }

    g_data.delay_line = (float*)malloc(g_data.delay_line_length * sizeof(float));
    if (!g_data.delay_line) {
        fprintf(stderr, "FATAL: Failed to allocate memory for delay line.\n");
        free(g_data.output_buffer);
        exit(1);
    }

    // Initialize the delay line with random noise (-0.5 to 0.5)
    // This is the initial "pluck" of the string.
    for (long i = 0; i < g_data.delay_line_length; ++i) {
        g_data.delay_line[i] = ((float)mt_rand() / (float)UINT32_MAX) - 0.5f;
    }

    g_data.result_accumulator = 0.0;
}

void run_computation() {
    long cursor = 0;
    long next_cursor;
    double accumulator = 0.0;

    for (long i = 0; i < g_data.num_samples; ++i) {
        // Determine the next sample in the circular buffer
        next_cursor = (cursor + 1) % g_data.delay_line_length;

        // Get the current sample from the delay line
        float current_sample = g_data.delay_line[cursor];

        // Store the current sample as the output
        g_data.output_buffer[i] = current_sample;

        // Create the new sample by averaging the current and next samples.
        // This is a simple low-pass filter that causes the sound to decay.
        float new_sample = (current_sample + g_data.delay_line[next_cursor]) * 0.5f;

        // Write the new, filtered sample back into the delay line
        g_data.delay_line[cursor] = new_sample;
        
        // Move the cursor forward
        cursor = next_cursor;
    }

    // Accumulate a result to prevent dead code elimination
    for (long i = 0; i < g_data.num_samples; ++i) {
        accumulator += g_data.output_buffer[i];
    }
    g_data.result_accumulator = accumulator;
}

void cleanup() {
    free(g_data.output_buffer);
    free(g_data.delay_line);
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

    // Print result to stdout
    printf("%f\n", g_data.result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
