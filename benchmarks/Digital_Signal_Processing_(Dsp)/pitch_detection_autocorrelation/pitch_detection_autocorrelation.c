#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
// --- END MERSENNE TWISTER ---

// Global structure to hold benchmark data
typedef struct {
    size_t num_samples;
    int sample_rate;
    float fundamental_frequency;
    float *audio_buffer;
    float *autocorrelation_result;
    float detected_frequency; // To store the final result
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <audio_buffer_size_ms> <sample_rate_hz> <fundamental_frequency_hz> <seed>\n", argv[0]);
        exit(1);
    }

    int audio_buffer_size_ms = atoi(argv[1]);
    g_data.sample_rate = atoi(argv[2]);
    g_data.fundamental_frequency = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    g_data.num_samples = (size_t)(((long)audio_buffer_size_ms * g_data.sample_rate) / 1000L);

    // Allocate memory on the heap
    g_data.audio_buffer = (float *)malloc(g_data.num_samples * sizeof(float));
    if (!g_data.audio_buffer) {
        fprintf(stderr, "Failed to allocate memory for audio_buffer\n");
        exit(1);
    }

    g_data.autocorrelation_result = (float *)malloc(g_data.num_samples * sizeof(float));
    if (!g_data.autocorrelation_result) {
        fprintf(stderr, "Failed to allocate memory for autocorrelation_result\n");
        free(g_data.audio_buffer);
        exit(1);
    }

    // Generate a sine wave with some noise
    for (size_t i = 0; i < g_data.num_samples; ++i) {
        double time = (double)i / g_data.sample_rate;
        float sine_wave = sin(2.0 * M_PI * g_data.fundamental_frequency * time);
        float noise = ((mt_rand() / (float)UINT32_MAX) - 0.5f) * 0.2f; // Low amplitude noise
        g_data.audio_buffer[i] = sine_wave + noise;
    }

    g_data.detected_frequency = 0.0f; // Initialize result
}

void run_computation() {
    size_t n_samples = g_data.num_samples;
    float *buffer = g_data.audio_buffer;
    float *autocorr = g_data.autocorrelation_result;

    // Compute autocorrelation for all possible lags
    for (size_t lag = 0; lag < n_samples; ++lag) {
        double sum = 0.0;
        for (size_t i = 0; i < n_samples - lag; ++i) {
            sum += buffer[i] * buffer[i + lag];
        }
        autocorr[lag] = (float)sum;
    }

    // Find the peak in the autocorrelation function to determine the fundamental period.
    // We must skip the peak at lag 0.
    // We define a minimum lag to start searching from, corresponding to a max detectable frequency.
    size_t min_lag = (size_t)(g_data.sample_rate / 2000.0f); // Max freq = 2000Hz
    if (min_lag >= n_samples) min_lag = 1; // Sanity check

    size_t max_lag = min_lag;
    float max_corr = 0.0f;

    for (size_t lag = min_lag; lag < n_samples; ++lag) {
        if (autocorr[lag] > max_corr) {
            max_corr = autocorr[lag];
            max_lag = lag;
        }
    }

    // Calculate frequency from the lag (period)
    if (max_lag > 0) {
        g_data.detected_frequency = (float)g_data.sample_rate / (float)max_lag;
    }
}

void cleanup() {
    free(g_data.audio_buffer);
    free(g_data.autocorrelation_result);
    g_data.audio_buffer = NULL;
    g_data.autocorrelation_result = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", g_data.detected_frequency);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
