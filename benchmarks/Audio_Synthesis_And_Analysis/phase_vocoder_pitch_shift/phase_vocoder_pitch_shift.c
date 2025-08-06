#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

// For M_PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (Provided) ---
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
// --- end Mersenne Twister ---


typedef struct {
    float real;
    float imag;
} Complex;

// --- Global Benchmark Data ---
int g_audio_length_seconds;
float g_pitch_shift_factor;
int g_fft_size;
int g_hop_size;
const int g_sample_rate = 44100;

int g_total_samples;
int g_num_frames;

// Signal buffers
float *g_input_signal;
float *g_output_signal;
float *g_window;

// Phase vocoder state
float *g_last_input_phase;
float *g_last_output_phase;

// Final Result
float g_final_result;

// --- Helper Functions ---

// Recursive Radix-2 Cooley-Tukey FFT
void fft_recursive(Complex *buf, Complex *out, int n, int step) {
    if (n == 1) {
        out[0] = buf[0];
        return;
    }

    // Recurse on even and odd parts
    fft_recursive(buf, out, n / 2, step * 2);
    fft_recursive(buf + step, out + n / 2, n / 2, step * 2);

    for (int k = 0; k < n / 2; k++) {
        float angle = -2.0f * M_PI * k / n;
        Complex t = {cosf(angle), sinf(angle)};
        Complex odd = out[k + n / 2];

        Complex temp_prod = {t.real * odd.real - t.imag * odd.imag, t.real * odd.imag + t.imag * odd.real};

        out[k + n / 2] = (Complex){out[k].real - temp_prod.real, out[k].imag - temp_prod.imag};
        out[k] = (Complex){out[k].real + temp_prod.real, out[k].imag + temp_prod.imag};
    }
}

void fft(Complex *buf, int n) {
    Complex* out = (Complex*)malloc(n * sizeof(Complex));
    if (!out) { fprintf(stderr, "FFT memory allocation failed\n"); exit(1); }
    fft_recursive(buf, out, n, 1);
    memcpy(buf, out, n * sizeof(Complex));
    free(out);
}

void ifft(Complex *buf, int n) {
    for (int i = 0; i < n; i++) {
        buf[i].imag = -buf[i].imag;
    }

    fft(buf, n);

    for (int i = 0; i < n; i++) {
        buf[i].imag = -buf[i].imag;
        buf[i].real /= n;
        buf[i].imag /= n;
    }
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s audio_length_seconds pitch_shift_factor fft_size hop_size seed\n", argv[0]);
        exit(1);
    }

    g_audio_length_seconds = atoi(argv[1]);
    g_pitch_shift_factor = atof(argv[2]);
    g_fft_size = atoi(argv[3]);
    g_hop_size = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    g_total_samples = g_audio_length_seconds * g_sample_rate;
    g_num_frames = (g_total_samples - g_fft_size) / g_hop_size + 1;

    g_input_signal = (float*)malloc(g_total_samples * sizeof(float));
    g_output_signal = (float*)malloc(g_total_samples * sizeof(float));
    g_window = (float*)malloc(g_fft_size * sizeof(float));
    g_last_input_phase = (float*)malloc(g_fft_size * sizeof(float));
    g_last_output_phase = (float*)malloc(g_fft_size * sizeof(float));

    if (!g_input_signal || !g_output_signal || !g_window || !g_last_input_phase || !g_last_output_phase) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    memset(g_output_signal, 0, g_total_samples * sizeof(float));
    memset(g_last_input_phase, 0, g_fft_size * sizeof(float));
    memset(g_last_output_phase, 0, g_fft_size * sizeof(float));

    for (int i = 0; i < g_fft_size; i++) {
        g_window[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (g_fft_size - 1)));
    }

    for (int i = 0; i < g_total_samples; i++) {
        float t = (float)i / g_sample_rate;
        float random_component = ((float)mt_rand() / (float)UINT32_MAX - 0.5f) * 0.01f;
        g_input_signal[i] = 0.5f * sinf(2.0f * M_PI * 220.0f * t) + 0.3f * sinf(2.0f * M_PI * 330.0f * t) + random_component;
    }
}

void run_computation() {
    Complex *frame_complex = (Complex*)malloc(g_fft_size * sizeof(Complex));
    float *magnitudes = (float*)malloc(g_fft_size * sizeof(float));
    float *frequencies = (float*)malloc(g_fft_size * sizeof(float));
    float *current_phase = (float*)malloc(g_fft_size * sizeof(float));

    float freq_per_bin = (float)g_sample_rate / g_fft_size;
    float expected_phase_advance = 2.0f * M_PI * g_hop_size / g_fft_size;

    for (int frame_idx = 0; frame_idx < g_num_frames; ++frame_idx) {
        int frame_start = frame_idx * g_hop_size;

        for (int i = 0; i < g_fft_size; ++i) {
            frame_complex[i].real = g_input_signal[frame_start + i] * g_window[i];
            frame_complex[i].imag = 0.0f;
        }

        fft(frame_complex, g_fft_size);

        for (int k = 0; k < g_fft_size; ++k) {
            magnitudes[k] = sqrtf(frame_complex[k].real * frame_complex[k].real + frame_complex[k].imag * frame_complex[k].imag);
            current_phase[k] = atan2f(frame_complex[k].imag, frame_complex[k].real);

            float phase_delta = current_phase[k] - g_last_input_phase[k];
            float phase_dev = phase_delta - k * expected_phase_advance;
            phase_dev = remainderf(phase_dev, 2.0f * M_PI);

            float true_freq = k * freq_per_bin + phase_dev * g_sample_rate / (2.0f * M_PI * g_hop_size);
            frequencies[k] = true_freq;
            g_last_input_phase[k] = current_phase[k];
        }

        for (int k = 0; k < g_fft_size; ++k) {
            float new_freq = frequencies[k] * g_pitch_shift_factor;
            g_last_output_phase[k] += 2.0f * M_PI * new_freq * g_hop_size / g_sample_rate;

            frame_complex[k].real = magnitudes[k] * cosf(g_last_output_phase[k]);
            frame_complex[k].imag = magnitudes[k] * sinf(g_last_output_phase[k]);
        }

        ifft(frame_complex, g_fft_size);
        
        for (int i = 0; i < g_fft_size; ++i) {
            if (frame_start + i < g_total_samples) {
                g_output_signal[frame_start + i] += frame_complex[i].real * g_window[i];
            }
        }
    }
    
    float sum = 0.0f;
    for (int i = 0; i < g_total_samples; i++) {
       sum += g_output_signal[i];
    }
    g_final_result = sum;

    free(frame_complex);
    free(magnitudes);
    free(frequencies);
    free(current_phase);
}

void cleanup() {
    free(g_input_signal);
    free(g_output_signal);
    free(g_window);
    free(g_last_input_phase);
    free(g_last_output_phase);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
