#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (MT19937) Generator --- DO NOT MODIFY
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

// Global structure to hold benchmark data
typedef struct {
    // Parameters
    int duration_seconds;
    int sample_rate_hz;
    int wavetable_size;
    int num_oscillators; // Mapped from num_interpolations argument

    // Data arrays
    float* wavetable;
    float* frequencies;
    float* phases;
    float* output_buffer;

    long total_samples;
    double accumulator;
} BenchmarkData;

static BenchmarkData g_data;

// Setup: Parse args, allocate memory, initialize data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s duration_seconds sample_rate_hz wavetable_size num_interpolations seed\n", argv[0]);
        exit(1);
    }

    g_data.duration_seconds = atoi(argv[1]);
    g_data.sample_rate_hz = atoi(argv[2]);
    g_data.wavetable_size = atoi(argv[3]);
    g_data.num_oscillators = atoi(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);

    mt_seed(seed);

    g_data.total_samples = (long)g_data.duration_seconds * g_data.sample_rate_hz;
    g_data.accumulator = 0.0;

    // Allocate memory
    g_data.wavetable = (float*)malloc(g_data.wavetable_size * sizeof(float));
    g_data.frequencies = (float*)malloc(g_data.num_oscillators * sizeof(float));
    g_data.phases = (float*)malloc(g_data.num_oscillators * sizeof(float));
    g_data.output_buffer = (float*)malloc(g_data.total_samples * sizeof(float));

    if (!g_data.wavetable || !g_data.frequencies || !g_data.phases || !g_data.output_buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize wavetable with a sine wave
    for (int i = 0; i < g_data.wavetable_size; ++i) {
        g_data.wavetable[i] = sinf(2.0f * (float)M_PI * i / g_data.wavetable_size);
    }

    // Initialize oscillators with random frequencies and zero phase
    for (int i = 0; i < g_data.num_oscillators; ++i) {
        // Generate frequency in a musically relevant range (e.g., 20Hz - 5000Hz)
        g_data.frequencies[i] = 20.0f + ((float)mt_rand() / (float)UINT32_MAX) * 4980.0f;
        g_data.phases[i] = 0.0f;
    }
    
    // Initialize output buffer to zero
    for (long i = 0; i < g_data.total_samples; ++i) {
        g_data.output_buffer[i] = 0.0f;
    }
}

// Computation: Synthesize audio using wavetable oscillators
void run_computation() {
    for (long i = 0; i < g_data.total_samples; ++i) {
        float sample_sum = 0.0f;

        for (int j = 0; j < g_data.num_oscillators; ++j) {
            // Get integer and fractional parts of the phase for interpolation
            float current_phase = g_data.phases[j];
            int index1 = (int)current_phase;
            float fraction = current_phase - index1;

            // Get the two table values for interpolation
            float value1 = g_data.wavetable[index1];
            int index2 = (index1 + 1) % g_data.wavetable_size; // Wrap around
            float value2 = g_data.wavetable[index2];

            // Linear interpolation
            float interpolated_sample = value1 + fraction * (value2 - value1);
            sample_sum += interpolated_sample;

            // Update phase for the next sample
            float phase_increment = (g_data.frequencies[j] * g_data.wavetable_size) / g_data.sample_rate_hz;
            g_data.phases[j] += phase_increment;
            if (g_data.phases[j] >= g_data.wavetable_size) {
                g_data.phases[j] -= g_data.wavetable_size;
            }
        }
        // Store the mixed sample (simple average to avoid clipping)
        g_data.output_buffer[i] = sample_sum / g_data.num_oscillators;
    }

    // Accumulate the final result to prevent dead code elimination
    for (long i = 0; i < g_data.total_samples; ++i) {
        g_data.accumulator += g_data.output_buffer[i];
    }
}

// Cleanup: Free all allocated memory
void cleanup() {
    free(g_data.wavetable);
    free(g_data.frequencies);
    free(g_data.phases);
    free(g_data.output_buffer);
    g_data.wavetable = NULL;
    g_data.frequencies = NULL;
    g_data.phases = NULL;
    g_data.output_buffer = NULL;
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
    printf("%f\n", g_data.accumulator);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
