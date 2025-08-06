#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- START: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// --- Benchmark Data and Globals ---
typedef struct {
    // Parameters
    int audio_length_seconds;
    int sample_rate_hz;
    int fft_size;
    int hop_size;
    long total_samples;
    int num_frames;

    // Data arrays
    float* audio_signal;
    float* window;
    float* frame;
    float* fft_real;
    float* fft_imag;
    float* magnitude_spectrum;

    // Result
    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s audio_length_seconds sample_rate_hz fft_size hop_size seed\n", argv[0]);
        exit(1);
    }

    g_data.audio_length_seconds = atoi(argv[1]);
    g_data.sample_rate_hz = atoi(argv[2]);
    g_data.fft_size = atoi(argv[3]);
    g_data.hop_size = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    g_data.total_samples = (long)g_data.audio_length_seconds * g_data.sample_rate_hz;
    if (g_data.total_samples < g_data.fft_size) {
        fprintf(stderr, "Error: Audio length must be at least one FFT size.\n");
        exit(1);
    }
    g_data.num_frames = (g_data.total_samples - g_data.fft_size) / g_data.hop_size + 1;

    // Seed the random number generator
    mt_seed(seed);

    // Allocate all memory on the heap
    g_data.audio_signal = (float*)malloc(g_data.total_samples * sizeof(float));
    g_data.window = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.frame = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.fft_real = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.fft_imag = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.magnitude_spectrum = (float*)malloc((g_data.fft_size / 2 + 1) * sizeof(float));

    if (!g_data.audio_signal || !g_data.window || !g_data.frame || !g_data.fft_real || !g_data.fft_imag || !g_data.magnitude_spectrum) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random audio signal (-1.0 to 1.0)
    for (long i = 0; i < g_data.total_samples; ++i) {
        g_data.audio_signal[i] = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }

    // Pre-calculate Hanning window
    for (int i = 0; i < g_data.fft_size; ++i) {
        g_data.window[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (g_data.fft_size - 1.0f)));
    }

    g_data.final_result = 0.0;
}

void run_computation() {
    double total_centroid_sum = 0.0;
    int fft_half = g_data.fft_size / 2 + 1;

    for (int i = 0; i < g_data.num_frames; ++i) {
        long current_pos = (long)i * g_data.hop_size;

        // 1. Apply window to the current frame
        for (int j = 0; j < g_data.fft_size; ++j) {
            g_data.frame[j] = g_data.audio_signal[current_pos + j] * g_data.window[j];
        }

        // 2. Compute DFT (slow, for benchmark purposes)
        for (int k = 0; k < fft_half; ++k) {
            double sum_real = 0.0;
            double sum_imag = 0.0;
            for (int n = 0; n < g_data.fft_size; ++n) {
                double angle = 2.0 * M_PI * k * n / g_data.fft_size;
                sum_real += g_data.frame[n] * cos(angle);
                sum_imag += g_data.frame[n] * -sin(angle);
            }
            g_data.fft_real[k] = (float)sum_real;
            g_data.fft_imag[k] = (float)sum_imag;
        }

        // 3. Compute magnitude spectrum
        for (int k = 0; k < fft_half; ++k) {
            g_data.magnitude_spectrum[k] = sqrtf(g_data.fft_real[k] * g_data.fft_real[k] + g_data.fft_imag[k] * g_data.fft_imag[k]);
        }

        // 4. Calculate spectral centroid for the frame
        double weighted_sum_freq = 0.0;
        double sum_magnitudes = 0.0;
        for (int k = 0; k < fft_half; ++k) {
            double freq_of_bin = (double)k * g_data.sample_rate_hz / g_data.fft_size;
            weighted_sum_freq += freq_of_bin * g_data.magnitude_spectrum[k];
            sum_magnitudes += g_data.magnitude_spectrum[k];
        }

        if (sum_magnitudes > 1e-9) {
            total_centroid_sum += weighted_sum_freq / sum_magnitudes;
        }
    }

    g_data.final_result = total_centroid_sum;
}

void cleanup() {
    free(g_data.audio_signal);
    free(g_data.window);
    free(g_data.frame);
    free(g_data.fft_real);
    free(g_data.fft_imag);
    free(g_data.magnitude_spectrum);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%f\n", g_data.final_result);

    // Print timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
