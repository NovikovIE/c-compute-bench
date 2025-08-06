#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

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

// Global benchmark data structure
typedef struct {
    // Parameters
    int audio_length_seconds;
    float stretch_factor;
    int fft_size;
    int hop_size;
    uint32_t seed;

    // Derived values
    int num_samples;
    int sample_rate;
    int num_bins;

    // Data arrays
    float *input_signal;
    float *window;
    float *prev_phase;
    float *current_phase;
    float *output_phase;
    float *magnitude;

    // Result accumulator
    double result_accumulator;
} BenchmarkData;

static BenchmarkData g_data;

// Function prototypes
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s audio_length_seconds stretch_factor fft_size hop_size seed\n", argv[0]);
        exit(1);
    }

    g_data.audio_length_seconds = atoi(argv[1]);
    g_data.stretch_factor = atof(argv[2]);
    g_data.fft_size = atoi(argv[3]);
    g_data.hop_size = atoi(argv[4]);
    g_data.seed = (uint32_t)atoi(argv[5]);

    mt_seed(g_data.seed);

    g_data.sample_rate = 44100;
    g_data.num_samples = g_data.audio_length_seconds * g_data.sample_rate;
    g_data.num_bins = g_data.fft_size / 2 + 1;
    g_data.result_accumulator = 0.0;
    
    // Allocate memory
    g_data.input_signal = (float*)malloc(g_data.num_samples * sizeof(float));
    g_data.window = (float*)malloc(g_data.fft_size * sizeof(float));
    g_data.prev_phase = (float*)malloc(g_data.num_bins * sizeof(float));
    g_data.current_phase = (float*)malloc(g_data.num_bins * sizeof(float));
    g_data.output_phase = (float*)malloc(g_data.num_bins * sizeof(float));
    g_data.magnitude = (float*)malloc(g_data.num_bins * sizeof(float));

    if (!g_data.input_signal || !g_data.window || !g_data.prev_phase || !g_data.current_phase || !g_data.output_phase || !g_data.magnitude) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random audio signal
    for (int i = 0; i < g_data.num_samples; i++) {
        g_data.input_signal[i] = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }

    // Generate Hann window
    for (int i = 0; i < g_data.fft_size; i++) {
        g_data.window[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (g_data.fft_size - 1)));
    }

    // Initialize phase arrays to zero
    memset(g_data.prev_phase, 0, g_data.num_bins * sizeof(float));
    memset(g_data.current_phase, 0, g_data.num_bins * sizeof(float));
    memset(g_data.output_phase, 0, g_data.num_bins * sizeof(float));
    memset(g_data.magnitude, 0, g_data.num_bins * sizeof(float));
}

void run_computation() {
    long num_frames = (g_data.num_samples - g_data.fft_size) / g_data.hop_size + 1;
    float *frame = (float*)malloc(g_data.fft_size * sizeof(float));
    if (!frame) {
        fprintf(stderr, "FATAL: Frame allocation failed in computation.\n");
        exit(1);
    }

    double synthesis_hop = g_data.hop_size * g_data.stretch_factor;

    for (long i = 0; i < num_frames; i++) {
        long frame_start = i * g_data.hop_size;

        // Apply window to frame
        for (int n = 0; n < g_data.fft_size; n++) {
            frame[n] = g_data.input_signal[frame_start + n] * g_data.window[n];
        }

        // --- Simulated FFT ---
        for (int k = 0; k < g_data.num_bins; k++) {
            double real = 0.0, imag = 0.0;
            for (int n = 0; n < 16; n++) {
                int sample_index = (n * (g_data.fft_size / 16)) % g_data.fft_size;
                double angle = 2.0 * M_PI * k * sample_index / g_data.fft_size;
                real += frame[sample_index] * cos(angle);
                imag -= frame[sample_index] * sin(angle);
            }
            g_data.magnitude[k] = sqrt(real * real + imag * imag);
            g_data.current_phase[k] = atan2(imag, real);
        }

        // --- Phase Vocoder Processing ---
        for (int k = 0; k < g_data.num_bins; k++) {
            double phase_dev = g_data.current_phase[k] - g_data.prev_phase[k];
            double expected_phase = (double)k * 2.0 * M_PI * g_data.hop_size / g_data.fft_size;
            double freq_dev = phase_dev - expected_phase;
            freq_dev = fmod(freq_dev + M_PI, 2.0 * M_PI) - M_PI;
            double true_freq_rad_per_sample = (expected_phase + freq_dev) / g_data.hop_size;
            
            g_data.output_phase[k] += synthesis_hop * true_freq_rad_per_sample;
            g_data.prev_phase[k] = g_data.current_phase[k];
        }

        // --- Simulated Inverse FFT and Accumulation ---
        for (int k = 0; k < g_data.num_bins; k++) {
             g_data.result_accumulator += g_data.magnitude[k] * (cosf(g_data.output_phase[k]) + sinf(g_data.output_phase[k]));
        }
    }
    free(frame);
}

void cleanup() {
    free(g_data.input_signal);
    free(g_data.window);
    free(g_data.prev_phase);
    free(g_data.current_phase);
    free(g_data.output_phase);
    free(g_data.magnitude);
}

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
