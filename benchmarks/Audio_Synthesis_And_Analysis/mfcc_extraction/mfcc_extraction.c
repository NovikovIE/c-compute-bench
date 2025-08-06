#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Mersenne Twister (provided)
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
// End of Mersenne Twister

// --- Global data structure ---
typedef struct {
    // Parameters
    int audio_length_seconds;
    int sample_rate_hz;
    int num_mfccs;
    int fft_size;
    int hop_size;
    uint32_t seed;

    // Derived values
    int num_samples;
    int num_frames;
    int spectrum_size;

    // Data arrays
    float *audio_signal;
    float *frame;
    float *window;
    float *spectrum;
    float *mel_filterbank;
    float *mel_energies;
    float *dct_matrix;
    float *mfccs; 
    
    // Result
    double result_accumulator;
} BenchmarkData;

BenchmarkData g_data;

// --- Helper functions ---
float hz_to_mel(float hz) {
    return 2595.0f * log10f(1.0f + hz / 700.0f);
}

float mel_to_hz(float mel) {
    return 700.0f * (powf(10.0f, mel / 2595.0f) - 1.0f);
}

// --- Benchmark functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s audio_length_seconds sample_rate_hz num_mfccs fft_size hop_size seed\n", argv[0]);
        exit(1);
    }

    g_data.audio_length_seconds = atoi(argv[1]);
    g_data.sample_rate_hz = atoi(argv[2]);
    g_data.num_mfccs = atoi(argv[3]);
    g_data.fft_size = atoi(argv[4]);
    g_data.hop_size = atoi(argv[5]);
    g_data.seed = (uint32_t)atoi(argv[6]);

    mt_seed(g_data.seed);

    // Calculate derived parameters
    g_data.num_samples = g_data.audio_length_seconds * g_data.sample_rate_hz;
    g_data.num_frames = (g_data.num_samples - g_data.fft_size) / g_data.hop_size + 1;
    g_data.spectrum_size = g_data.fft_size / 2 + 1;

    // Allocate memory
    g_data.audio_signal = (float*)malloc(g_data.num_samples * sizeof(float));
    g_data.frame = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.window = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.spectrum = (float*)malloc(g_data.spectrum_size * sizeof(float));
    g_data.mel_filterbank = (float*)malloc(g_data.num_mfccs * g_data.spectrum_size * sizeof(float));
    g_data.mel_energies = (float*)malloc(g_data.num_mfccs * sizeof(float));
    g_data.dct_matrix = (float*)malloc(g_data.num_mfccs * g_data.num_mfccs * sizeof(float));
    g_data.mfccs = (float*)malloc(g_data.num_frames * g_data.num_mfccs * sizeof(float));

    if (!g_data.audio_signal || !g_data.frame || !g_data.window || !g_data.spectrum ||
        !g_data.mel_filterbank || !g_data.mel_energies || !g_data.dct_matrix || !g_data.mfccs) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize data
    // 1. Audio signal (random noise)
    for (int i = 0; i < g_data.num_samples; ++i) {
        g_data.audio_signal[i] = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }

    // 2. Hamming window
    for (int i = 0; i < g_data.fft_size; ++i) {
        g_data.window[i] = 0.54f - 0.46f * cosf(2.0f * M_PI * i / (g_data.fft_size - 1));
    }
    
    // 3. Mel filterbank
    float low_freq_mel = hz_to_mel(0);
    float high_freq_mel = hz_to_mel(g_data.sample_rate_hz / 2.0f);
    float* mel_points = (float*)malloc((g_data.num_mfccs + 2) * sizeof(float));
    float* hz_points = (float*)malloc((g_data.num_mfccs + 2) * sizeof(float));
    int* fft_bin = (int*)malloc((g_data.num_mfccs + 2) * sizeof(int));

    for (int i = 0; i < g_data.num_mfccs + 2; ++i) {
        mel_points[i] = low_freq_mel + i * (high_freq_mel - low_freq_mel) / (g_data.num_mfccs + 1);
        hz_points[i] = mel_to_hz(mel_points[i]);
        fft_bin[i] = floor((g_data.fft_size + 1) * hz_points[i] / g_data.sample_rate_hz);
    }

    for (int i = 0; i < g_data.num_mfccs * g_data.spectrum_size; ++i) g_data.mel_filterbank[i] = 0.0f;

    for (int j = 0; j < g_data.num_mfccs; ++j) {
        for (int i = fft_bin[j]; i < fft_bin[j + 1]; ++i) {
            if (fft_bin[j + 1] > fft_bin[j])
                g_data.mel_filterbank[j * g_data.spectrum_size + i] = (float)(i - fft_bin[j]) / (float)(fft_bin[j + 1] - fft_bin[j]);
        }
        for (int i = fft_bin[j + 1]; i < fft_bin[j + 2]; ++i) {
             if (fft_bin[j + 2] > fft_bin[j + 1])
                g_data.mel_filterbank[j * g_data.spectrum_size + i] = (float)(fft_bin[j + 2] - i) / (float)(fft_bin[j + 2] - fft_bin[j + 1]);
        }
    }
    free(mel_points);
    free(hz_points);
    free(fft_bin);

    // 4. DCT matrix (DCT-II)
    for (int k = 0; k < g_data.num_mfccs; ++k) {
        for (int n = 0; n < g_data.num_mfccs; ++n) {
            float scale = (k == 0) ? sqrtf(1.0f / g_data.num_mfccs) : sqrtf(2.0f / g_data.num_mfccs);
            g_data.dct_matrix[k * g_data.num_mfccs + n] = scale * cosf(M_PI * k * (2 * n + 1) / (2.0f * g_data.num_mfccs));
        }
    }
}

void run_computation() {
    g_data.result_accumulator = 0.0;
    
    for (int i = 0; i < g_data.num_frames; ++i) {
        int frame_start = i * g_data.hop_size;

        // 1. Framing and Windowing
        for (int j = 0; j < g_data.fft_size; ++j) {
            g_data.frame[j] = g_data.audio_signal[frame_start + j] * g_data.window[j];
        }

        // 2. DFT Magnitude (explicit calculation for workload)
        for (int k = 0; k < g_data.spectrum_size; ++k) {
            float real = 0.0f, imag = 0.0f;
            for (int n = 0; n < g_data.fft_size; ++n) {
                float angle = -2.0f * M_PI * k * n / g_data.fft_size;
                real += g_data.frame[n] * cosf(angle);
                imag += g_data.frame[n] * sinf(angle);
            }
            g_data.spectrum[k] = sqrtf(real * real + imag * imag);
        }

        // 3. Apply Mel Filterbank & Get Log Power
        for (int j = 0; j < g_data.num_mfccs; ++j) {
            float energy = 0.0f;
            for (int k = 0; k < g_data.spectrum_size; ++k) {
                energy += g_data.spectrum[k] * g_data.mel_filterbank[j * g_data.spectrum_size + k];
            }
            g_data.mel_energies[j] = logf(energy + 1e-6f);
        }

        // 4. DCT
        for (int k = 0; k < g_data.num_mfccs; ++k) {
            float mfcc_val = 0.0f;
            for (int j = 0; j < g_data.num_mfccs; ++j) {
                mfcc_val += g_data.mel_energies[j] * g_data.dct_matrix[k * g_data.num_mfccs + j];
            }
            g_data.mfccs[i * g_data.num_mfccs + k] = mfcc_val;
        }
    }

    // Accumulate result to prevent dead code elimination
    for (int i = 0; i < g_data.num_frames * g_data.num_mfccs; ++i) {
        g_data.result_accumulator += g_data.mfccs[i];
    }
}

void cleanup() {
    free(g_data.audio_signal);
    free(g_data.frame);
    free(g_data.window);
    free(g_data.spectrum);
    free(g_data.mel_filterbank);
    free(g_data.mel_energies);
    free(g_data.dct_matrix);
    free(g_data.mfccs);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", g_data.result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
