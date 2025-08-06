#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark: onset_detection_spectral_flux ---

typedef struct {
    double real;
    double imag;
} Complex;

// --- Global Benchmark Data ---

// Parameters
int g_audio_length_seconds;
int g_sample_rate_hz;
int g_fft_size;
int g_hop_size;
long g_num_samples;
long g_num_frames;

// Data Buffers
float *g_audio_signal = NULL;
float *g_window = NULL;
float *g_prev_spectrum = NULL;
float *g_curr_spectrum = NULL;
Complex *g_fft_buffer = NULL;

// Final Result
double g_final_result = 0.0;

// --- Forward Declarations ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
void fft(Complex *buf, int n);
static inline float random_float(float min, float max);

// --- Main Execution --- 
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    printf("%f\n", g_final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}

// --- Benchmark Implementation ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s audio_length_seconds sample_rate_hz fft_size hop_size seed\n", argv[0]);
        exit(1);
    }

    g_audio_length_seconds = atoi(argv[1]);
    g_sample_rate_hz = atoi(argv[2]);
    g_fft_size = atoi(argv[3]);
    g_hop_size = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);
    mt_seed(seed);

    // Validate parameters
    assert(g_audio_length_seconds > 0);
    assert(g_sample_rate_hz > 0);
    assert(g_fft_size > 0);
    assert((g_fft_size & (g_fft_size - 1)) == 0); // Must be a power of 2
    assert(g_hop_size > 0);

    g_num_samples = (long)g_audio_length_seconds * g_sample_rate_hz;
    assert(g_num_samples >= g_fft_size);
    g_num_frames = 1 + (g_num_samples - g_fft_size) / g_hop_size;

    // Allocate memory
    g_audio_signal = (float*)malloc(g_num_samples * sizeof(float));
    g_window = (float*)malloc(g_fft_size * sizeof(float));
    int spectrum_size = g_fft_size / 2 + 1;
    g_prev_spectrum = (float*)malloc(spectrum_size * sizeof(float));
    g_curr_spectrum = (float*)malloc(spectrum_size * sizeof(float));
    g_fft_buffer = (Complex*)malloc(g_fft_size * sizeof(Complex));
    
    if (!g_audio_signal || !g_window || !g_prev_spectrum || !g_curr_spectrum || !g_fft_buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate random audio signal
    for (long i = 0; i < g_num_samples; i++) {
        g_audio_signal[i] = random_float(-1.0f, 1.0f);
    }

    // Generate Hann window
    for (int i = 0; i < g_fft_size; i++) {
        g_window[i] = 0.5f * (1.0f - cos(2.0 * M_PI * i / (g_fft_size - 1.0)));
    }
}

void run_computation() {
    int spectrum_size = g_fft_size / 2 + 1;
    g_final_result = 0.0;

    memset(g_prev_spectrum, 0, spectrum_size * sizeof(float));

    for (long frame_idx = 0; frame_idx < g_num_frames; ++frame_idx) {
        long sample_start = frame_idx * g_hop_size;

        for (int i = 0; i < g_fft_size; ++i) {
            g_fft_buffer[i].real = g_audio_signal[sample_start + i] * g_window[i];
            g_fft_buffer[i].imag = 0.0;
        }

        fft(g_fft_buffer, g_fft_size);

        for (int i = 0; i < spectrum_size; ++i) {
            g_curr_spectrum[i] = sqrt(g_fft_buffer[i].real * g_fft_buffer[i].real + g_fft_buffer[i].imag * g_fft_buffer[i].imag);
        }

        double frame_flux = 0.0;
        for (int i = 0; i < spectrum_size; ++i) {
            double diff = g_curr_spectrum[i] - g_prev_spectrum[i];
            if (diff > 0) {
                frame_flux += diff;
            }
        }

        g_final_result += frame_flux;

        memcpy(g_prev_spectrum, g_curr_spectrum, spectrum_size * sizeof(float));
    }
}

void cleanup() {
    free(g_audio_signal);
    free(g_window);
    free(g_prev_spectrum);
    free(g_curr_spectrum);
    free(g_fft_buffer);
    g_audio_signal = NULL;
    g_window = NULL;
    g_prev_spectrum = NULL;
    g_curr_spectrum = NULL;
    g_fft_buffer = NULL;
}


// --- Utility Functions ---

static inline float random_float(float min, float max) {
    return min + ((float)mt_rand() / (float)UINT32_MAX) * (max - min);
}

void fft(Complex *buf, int n) {
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            Complex temp = buf[i];
            buf[i] = buf[j];
            buf[j] = temp;
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * M_PI / len;
        Complex wlen = {cos(ang), sin(ang)};
        for (int i = 0; i < n; i += len) {
            Complex w = {1.0, 0.0};
            for (int j = 0; j < len / 2; j++) {
                Complex u = buf[i + j];
                Complex v = {buf[i + j + len / 2].real * w.real - buf[i + j + len / 2].imag * w.imag,
                             buf[i + j + len / 2].real * w.imag + buf[i + j + len / 2].imag * w.real};
                
                buf[i + j].real = u.real + v.real;
                buf[i + j].imag = u.imag + v.imag;
                buf[i + j + len / 2].real = u.real - v.real;
                buf[i + j + len / 2].imag = u.imag - v.imag;

                double w_real_temp = w.real * wlen.real - w.imag * wlen.imag;
                w.imag = w.real * wlen.imag + w.imag * wlen.real;
                w.real = w_real_temp;
            }
        }
    }
}
