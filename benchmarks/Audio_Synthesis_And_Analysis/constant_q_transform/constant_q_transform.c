/*
 * constant_q_transform: A benchmark for audio analysis.
 * 
 * This program implements a Constant-Q Transform (CQT), a time-frequency
 * analysis tool used in music and audio processing. Unlike the standard
 * Fourier Transform, CQT uses geometrically spaced frequency bins, which
 * mirrors the human perception of pitch.
 *
 * The core computation involves convolving the input audio signal with a set of
 * pre-computed, windowed complex sinusoids (kernels), one for each frequency
 * bin. The length of each kernel is inversely proportional to its frequency,
 * ensuring a constant "quality factor" (Q) across all bins. This makes the
 * algorithm computationally demanding, especially for low frequencies which
 * require long analysis windows.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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

// --- Benchmark Globals ---

// Parameters
int audio_length_seconds;
int sample_rate_hz;
int bins_per_octave;

// Derived sizes and constants
int num_samples;
int num_bins;
const double F_MIN = 80.0; // Minimum frequency for CQT, e.g., ~E2

// Data arrays
float* audio_signal;
float** cqt_real;
float** cqt_imag;

// CQT Kernel data
float** kernel_real;
float** kernel_imag;
int* kernel_lengths;

// Result accumulator
double final_result = 0.0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <audio_length_seconds> <sample_rate_hz> <bins_per_octave> <seed>\n", argv[0]);
        exit(1);
    }

    audio_length_seconds = atoi(argv[1]);
    sample_rate_hz = atoi(argv[2]);
    bins_per_octave = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // 1. Calculate derived parameters
    num_samples = audio_length_seconds * sample_rate_hz;
    double f_max = sample_rate_hz / 2.0;
    num_bins = (int)floor(bins_per_octave * log2(f_max / F_MIN));
    double quality_factor = 1.0 / (pow(2.0, 1.0 / bins_per_octave) - 1.0);

    // 2. Allocate and initialize input audio signal
    audio_signal = (float*)malloc(num_samples * sizeof(float));
    if (!audio_signal) { perror("malloc audio_signal"); exit(1); }
    for (int i = 0; i < num_samples; ++i) {
        audio_signal[i] = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }

    // 3. Allocate and pre-compute CQT kernels
    kernel_lengths = (int*)malloc(num_bins * sizeof(int));
    kernel_real = (float**)malloc(num_bins * sizeof(float*));
    kernel_imag = (float**)malloc(num_bins * sizeof(float*));
    if (!kernel_lengths || !kernel_real || !kernel_imag) { perror("malloc kernels"); exit(1); }

    for (int k = 0; k < num_bins; ++k) {
        double freq_k = F_MIN * pow(2.0, (double)k / bins_per_octave);
        int N_k = (int)floor(quality_factor * sample_rate_hz / freq_k);
        if (N_k == 0) N_k = 1;

        kernel_lengths[k] = N_k;
        kernel_real[k] = (float*)malloc(N_k * sizeof(float));
        kernel_imag[k] = (float*)malloc(N_k * sizeof(float));
        if (!kernel_real[k] || !kernel_imag[k]) { perror("malloc kernel bin"); exit(1); }

        for (int n = 0; n < N_k; ++n) {
            double hamming_window = 0.54 - 0.46 * cos(2.0 * M_PI * n / (N_k - 1));
            double angle = -2.0 * M_PI * freq_k * n / sample_rate_hz;
            kernel_real[k][n] = (float)(hamming_window * cos(angle) / N_k);
            kernel_imag[k][n] = (float)(hamming_window * sin(angle) / N_k);
        }
    }

    // 4. Allocate CQT output matrices
    cqt_real = (float**)malloc(num_bins * sizeof(float*));
    cqt_imag = (float**)malloc(num_bins * sizeof(float*));
    if (!cqt_real || !cqt_imag) { perror("malloc cqt_matrix"); exit(1); }
    for (int k = 0; k < num_bins; k++) {
        cqt_real[k] = (float*)malloc(num_samples * sizeof(float));
        cqt_imag[k] = (float*)malloc(num_samples * sizeof(float));
        if (!cqt_real[k] || !cqt_imag[k]) { perror("malloc cqt_matrix bin"); exit(1); }
    }
}

void run_computation() {
    for (int k = 0; k < num_bins; k++) {
        int N_k = kernel_lengths[k];
        int half_win = N_k / 2;
        float* k_real = kernel_real[k];
        float* k_imag = kernel_imag[k];

        for (int t = 0; t < num_samples; t++) {
            double sum_real = 0.0;
            double sum_imag = 0.0;

            for (int n = 0; n < N_k; n++) {
                int signal_idx = t + n - half_win;
                if (signal_idx >= 0 && signal_idx < num_samples) {
                    float signal_val = audio_signal[signal_idx];
                    sum_real += signal_val * k_real[n];
                    sum_imag += signal_val * k_imag[n];
                }
            }
            cqt_real[k][t] = (float)sum_real;
            cqt_imag[k][t] = (float)sum_imag;
        }
    }

    // Accumulate a final result to prevent dead code elimination
    double mag_sum = 0.0;
    for (int k = 0; k < num_bins; k += 4) { // Stride to reduce accumulation overhead
        for (int t = 0; t < num_samples; t += 4) {
            float r = cqt_real[k][t];
            float i = cqt_imag[k][t];
            mag_sum += sqrt(r * r + i * i);
        }
    }
    final_result = mag_sum;
}

void cleanup() {
    for (int k = 0; k < num_bins; k++) {
        free(cqt_real[k]);
        free(cqt_imag[k]);
        free(kernel_real[k]);
        free(kernel_imag[k]);
    }
    free(cqt_real);
    free(cqt_imag);
    free(kernel_real);
    free(kernel_imag);
    free(kernel_lengths);
    free(audio_signal);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
