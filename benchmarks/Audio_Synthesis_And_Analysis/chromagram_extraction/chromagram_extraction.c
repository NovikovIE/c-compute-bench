#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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


// --- BENCHMARK DATA AND FUNCTIONS --- 

#define NUM_CHROMA_BINS 12

typedef struct {
    // Parameters
    int audio_length_seconds;
    int sample_rate_hz;
    int fft_size;
    int hop_size;

    // Derived values
    long num_samples;
    int num_frames;

    // Data arrays
    float *audio_signal;
    float *window_function;
    float *frame_buffer;
    float *fft_real;
    float *fft_imag;
    float *fft_magnitude;
    float *chroma_vector;
    float **chromagram;
    int *bin_to_chroma_map;

    // Final result
    double final_checksum;
} BenchmarkData;

static BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    printf("%f\n", g_data.final_checksum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s audio_length_seconds sample_rate_hz fft_size hop_size seed\n", argv[0]);
        exit(1);
    }

    g_data.audio_length_seconds = atoi(argv[1]);
    g_data.sample_rate_hz = atoi(argv[2]);
    g_data.fft_size = atoi(argv[3]);
    g_data.hop_size = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Calculate derived parameters
    g_data.num_samples = (long)g_data.audio_length_seconds * g_data.sample_rate_hz;
    if (g_data.num_samples < g_data.fft_size) {
        fprintf(stderr, "FATAL: Audio length is shorter than FFT size.\n");
        exit(1);
    }
    g_data.num_frames = 1 + (g_data.num_samples - g_data.fft_size) / g_data.hop_size;
    int fft_half_size = g_data.fft_size / 2 + 1;

    // Allocate memory
    g_data.audio_signal = (float*)malloc(g_data.num_samples * sizeof(float));
    g_data.window_function = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.frame_buffer = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.fft_real = (float*)malloc(fft_half_size * sizeof(float));
    g_data.fft_imag = (float*)malloc(fft_half_size * sizeof(float));
    g_data.fft_magnitude = (float*)malloc(fft_half_size * sizeof(float));
    g_data.chroma_vector = (float*)malloc(NUM_CHROMA_BINS * sizeof(float));
    g_data.bin_to_chroma_map = (int*)malloc(fft_half_size * sizeof(int));
    g_data.chromagram = (float**)malloc(g_data.num_frames * sizeof(float*));
    for (int i = 0; i < g_data.num_frames; ++i) {
        g_data.chromagram[i] = (float*)malloc(NUM_CHROMA_BINS * sizeof(float));
    }

    // Initialize data
    // Generate a random audio signal
    for (long i = 0; i < g_data.num_samples; ++i) {
        g_data.audio_signal[i] = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }

    // Create a Hann window
    for (int i = 0; i < g_data.fft_size; ++i) {
        g_data.window_function[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (g_data.fft_size - 1)));
    }

    // Pre-calculate the mapping from FFT bin to chroma bin
    const float C0_FREQ = 16.3516f; // Reference frequency for note C0
    for (int k = 0; k < fft_half_size; ++k) {
        float freq = (float)k * g_data.sample_rate_hz / g_data.fft_size;
        if (freq <= 0.0f) {
            g_data.bin_to_chroma_map[k] = -1; // Ignore DC and invalid frequencies
        } else {
            int chroma_index = (int)roundf(12.0f * log2f(freq / C0_FREQ));
            g_data.bin_to_chroma_map[k] = (chroma_index % NUM_CHROMA_BINS + NUM_CHROMA_BINS) % NUM_CHROMA_BINS;
        }
    }
}

void run_computation() {
    int fft_half_size = g_data.fft_size / 2 + 1;
    g_data.final_checksum = 0.0;

    for (int i = 0; i < g_data.num_frames; ++i) {
        long frame_start = (long)i * g_data.hop_size;

        // 1. Copy to buffer and apply window function
        for (int n = 0; n < g_data.fft_size; ++n) {
            g_data.frame_buffer[n] = g_data.audio_signal[frame_start + n] * g_data.window_function[n];
        }

        // 2. Perform Discrete Fourier Transform (DFT) - this is the bottleneck
        for (int k = 0; k < fft_half_size; ++k) {
            g_data.fft_real[k] = 0.0f;
            g_data.fft_imag[k] = 0.0f;
            for (int n = 0; n < g_data.fft_size; ++n) {
                float angle = -2.0f * M_PI * k * n / g_data.fft_size;
                g_data.fft_real[k] += g_data.frame_buffer[n] * cosf(angle);
                g_data.fft_imag[k] += g_data.frame_buffer[n] * sinf(angle);
            }
            g_data.fft_magnitude[k] = sqrtf(g_data.fft_real[k] * g_data.fft_real[k] + g_data.fft_imag[k] * g_data.fft_imag[k]);
        }

        // 3. Bin FFT magnitudes into chroma bins
        for(int b = 0; b < NUM_CHROMA_BINS; ++b) {
            g_data.chroma_vector[b] = 0.0f;
        }
        for (int k = 1; k < fft_half_size; ++k) { // Start at 1 to ignore DC component
            int chroma_bin = g_data.bin_to_chroma_map[k];
            if (chroma_bin != -1) {
                g_data.chroma_vector[chroma_bin] += g_data.fft_magnitude[k];
            }
        }

        // 4. Store the resulting chroma vector in the final chromagram
        for(int b = 0; b < NUM_CHROMA_BINS; ++b) {
             g_data.chromagram[i][b] = g_data.chroma_vector[b];
        }
    }

    // 5. Calculate a final checksum to prevent dead code elimination
    double checksum = 0.0;
    for (int i = 0; i < g_data.num_frames; ++i) {
        for (int b = 0; b < NUM_CHROMA_BINS; ++b) {
            checksum += g_data.chromagram[i][b];
        }
    }
    g_data.final_checksum = checksum;
}

void cleanup() {
    for (int i = 0; i < g_data.num_frames; ++i) {
        free(g_data.chromagram[i]);
    }
    free(g_data.chromagram);
    free(g_data.bin_to_chroma_map);
    free(g_data.chroma_vector);
    free(g_data.fft_magnitude);
    free(g_data.fft_imag);
    free(g_data.fft_real);
    free(g_data.frame_buffer);
    free(g_data.window_function);
    free(g_data.audio_signal);
}
