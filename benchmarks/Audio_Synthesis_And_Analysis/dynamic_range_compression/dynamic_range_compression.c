#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
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

// Benchmark-specific data structures and globals
typedef struct {
    // Parameters from command line
    int audio_length_seconds;
    int sample_rate_hz;
    float threshold_db;
    float ratio;
    float attack_ms;
    float release_ms;
    float makeup_gain_db;
    uint32_t seed;

    // Derived data
    size_t num_samples;
    float attack_coeff;
    float release_coeff;
    float makeup_gain_linear;

    // Data arrays
    float* input_signal;
    float* output_signal;

    // Result accumulator
    double result_accumulator;
} BenchmarkData;

static BenchmarkData g_data;

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 9) {
        fprintf(stderr, "Usage: %s audio_length_seconds sample_rate_hz threshold_db ratio attack_ms release_ms makeup_gain_db seed\n", argv[0]);
        exit(1);
    }

    // Parse arguments
    g_data.audio_length_seconds = atoi(argv[1]);
    g_data.sample_rate_hz = atoi(argv[2]);
    g_data.threshold_db = atof(argv[3]);
    g_data.ratio = atof(argv[4]);
    g_data.attack_ms = atof(argv[5]);
    g_data.release_ms = atof(argv[6]);
    g_data.makeup_gain_db = atof(argv[7]);
    g_data.seed = (uint32_t)strtoul(argv[8], NULL, 10);

    // Seed the random number generator
    mt_seed(g_data.seed);

    // Calculate derived parameters
    g_data.num_samples = (size_t)g_data.audio_length_seconds * g_data.sample_rate_hz;
    g_data.makeup_gain_linear = powf(10.0f, g_data.makeup_gain_db / 20.0f);
    
    // Attack and release coefficients are calculated for a one-pole filter.
    // A value of 0.0 results in an infinite time, so we add a small epsilon.
    g_data.attack_coeff = expf(-1.0f / (0.001f * g_data.attack_ms * g_data.sample_rate_hz + 1e-6f));
    g_data.release_coeff = expf(-1.0f / (0.001f * g_data.release_ms * g_data.sample_rate_hz + 1e-6f));

    // Allocate memory
    g_data.input_signal = (float*)malloc(g_data.num_samples * sizeof(float));
    g_data.output_signal = (float*)malloc(g_data.num_samples * sizeof(float));
    if (!g_data.input_signal || !g_data.output_signal) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random input audio signal (-1.0 to 1.0)
    for (size_t i = 0; i < g_data.num_samples; ++i) {
        g_data.input_signal[i] = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }

    g_data.result_accumulator = 0.0;
}

void run_computation() {
    float envelope_db = 0.0f; // This tracks the current gain reduction in dB

    for (size_t i = 0; i < g_data.num_samples; ++i) {
        float input_sample = g_data.input_signal[i];

        // 1. Level Detection: Get the absolute value of the input signal.
        // A small epsilon prevents taking the log of zero.
        float input_level = fabsf(input_sample) + 1e-9f;

        // 2. Gain Computation: Determine the desired gain reduction.
        float input_level_db = 20.0f * log10f(input_level);
        float over_threshold_db = input_level_db - g_data.threshold_db;

        float target_gain_db = 0.0f;
        if (over_threshold_db > 0.0f) {
            // The classic compressor formula
            target_gain_db = -over_threshold_db * (1.0f - 1.0f / g_data.ratio);
        }

        // 3. Ballistics: Smooth the gain change using attack and release times.
        // The envelope moves towards the target gain reduction.
        if (target_gain_db < envelope_db) { // Attack phase
            envelope_db = g_data.attack_coeff * envelope_db + (1.0f - g_data.attack_coeff) * target_gain_db;
        } else { // Release phase
            envelope_db = g_data.release_coeff * envelope_db + (1.0f - g_data.release_coeff) * target_gain_db;
        }

        // 4. Gain Application: Apply the computed gain.
        // Convert the gain from dB to a linear multiplier.
        float gain_linear = powf(10.0f, envelope_db / 20.0f);
        g_data.output_signal[i] = input_sample * gain_linear * g_data.makeup_gain_linear;

        // Accumulate the result to prevent dead code elimination
        g_data.result_accumulator += g_data.output_signal[i];
    }
}

void cleanup() {
    free(g_data.input_signal);
    free(g_data.output_signal);
}

int main(int argc, char* argv[]) {
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
