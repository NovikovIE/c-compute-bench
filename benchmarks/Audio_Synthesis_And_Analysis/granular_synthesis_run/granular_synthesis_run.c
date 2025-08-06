#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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

// --- Benchmark-specific code ---

// Global parameters and data structures
const int SAMPLE_RATE = 44100;

float *input_audio_buffer;
float *output_audio_buffer;
float *grain_window;

long p_output_duration_seconds;
long p_input_audio_length;
long p_num_grains_per_second;
long p_grain_size_ms;

long output_audio_length;
long grain_size_samples;
long total_grains;

double final_result; // Use double for accumulator to avoid precision loss

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s output_duration_seconds input_audio_length num_grains_per_second grain_size_ms seed\n", argv[0]);
        exit(1);
    }

    p_output_duration_seconds = atol(argv[1]);
    p_input_audio_length = atol(argv[2]);
    p_num_grains_per_second = atol(argv[3]);
    p_grain_size_ms = atol(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);
    mt_seed(seed);

    // Derived parameters
    output_audio_length = p_output_duration_seconds * SAMPLE_RATE;
    grain_size_samples = p_grain_size_ms * SAMPLE_RATE / 1000;
    total_grains = p_output_duration_seconds * p_num_grains_per_second;

    if (grain_size_samples > p_input_audio_length) {
        fprintf(stderr, "Error: grain size cannot be larger than input audio length.\n");
        exit(1);
    }

    // Allocate memory
    input_audio_buffer = (float *)malloc(p_input_audio_length * sizeof(float));
    output_audio_buffer = (float *)malloc(output_audio_length * sizeof(float));
    grain_window = (float *)malloc(grain_size_samples * sizeof(float));

    if (!input_audio_buffer || !output_audio_buffer || !grain_window) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Initialize input buffer with random data (like white noise)
    for (long i = 0; i < p_input_audio_length; ++i) {
        input_audio_buffer[i] = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }

    // Initialize output buffer to silence
    for (long i = 0; i < output_audio_length; ++i) {
        output_audio_buffer[i] = 0.0f;
    }

    // Pre-calculate a triangular window for the grains to prevent clicking
    long half_grain = grain_size_samples / 2;
    for (long i = 0; i < half_grain; ++i) {
        grain_window[i] = (float)i / (float)half_grain;
    }
    for (long i = half_grain; i < grain_size_samples; ++i) {
        grain_window[i] = 1.0f - ((float)(i - half_grain) / (float)half_grain);
    }
}

void run_computation() {
    long max_source_start = p_input_audio_length - grain_size_samples;
    long max_dest_start = output_audio_length - grain_size_samples;

    for (long i = 0; i < total_grains; ++i) {
        long source_start_index = mt_rand() % (max_source_start + 1);
        long dest_start_index = mt_rand() % (max_dest_start + 1);

        // Additive synthesis: overlay windowed grains
        for (long j = 0; j < grain_size_samples; ++j) {
            output_audio_buffer[dest_start_index + j] += 
                input_audio_buffer[source_start_index + j] * grain_window[j];
        }
    }

    // Calculate a checksum (sum of squares) to prevent dead code elimination
    final_result = 0.0;
    for (long i = 0; i < output_audio_length; ++i) {
        final_result += (double)output_audio_buffer[i] * (double)output_audio_buffer[i];
    }
}

void cleanup() {
    free(input_audio_buffer);
    free(output_audio_buffer);
    free(grain_window);
}

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

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
