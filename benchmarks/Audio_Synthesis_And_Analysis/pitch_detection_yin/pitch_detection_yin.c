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

// --- BENCHMARK DATA AND PARAMETERS ---
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Parameters
static int g_audio_length_seconds;
static int g_sample_rate_hz;
static int g_buffer_size;

// Data structures
static long g_total_samples;
static float *g_audio_signal = NULL;
static float *g_yin_buffer = NULL;

// Final result to prevent dead code elimination
static double g_final_result = 0.0;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <audio_length_seconds> <sample_rate_hz> <buffer_size> <seed>\n", argv[0]);
        exit(1);
    }

    g_audio_length_seconds = atoi(argv[1]);
    g_sample_rate_hz = atoi(argv[2]);
    g_buffer_size = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    g_total_samples = (long)g_audio_length_seconds * g_sample_rate_hz;

    g_audio_signal = (float *)malloc(g_total_samples * sizeof(float));
    if (!g_audio_signal) {
        fprintf(stderr, "Failed to allocate memory for audio signal.\n");
        exit(1);
    }

    g_yin_buffer = (float *)malloc((g_buffer_size / 2) * sizeof(float));
    if (!g_yin_buffer) {
        fprintf(stderr, "Failed to allocate memory for YIN buffer.\n");
        free(g_audio_signal);
        exit(1);
    }

    // Generate a synthetic audio signal: a sine wave with sweeping frequency and noise
    double start_freq = 110.0; // A2 note
    double end_freq = 440.0;   // A4 note
    double phase = 0.0;

    for (long i = 0; i < g_total_samples; ++i) {
        double progress = (double)i / (double)g_total_samples;
        double current_freq = start_freq + (end_freq - start_freq) * progress * progress;
        phase += 2.0 * M_PI * current_freq / g_sample_rate_hz;
        if (phase > 2.0 * M_PI) {
            phase -= 2.0 * M_PI;
        }
        float signal = sinf(phase);
        float noise = ((float)mt_rand() / (float)UINT32_MAX - 0.5f) * 0.1f; // Add some noise
        g_audio_signal[i] = signal + noise;
    }
}

void run_computation() {
    const float yin_threshold = 0.15f;
    int yin_buffer_size = g_buffer_size / 2;
    double pitch_sum = 0.0;
    int windows_processed = 0;

    int hop_size = g_buffer_size / 2; // 50% overlap

    for (long i = 0; i <= g_total_samples - g_buffer_size; i += hop_size) {
        float *window = &g_audio_signal[i];

        // 1. Difference function
        for (int tau = 0; tau < yin_buffer_size; ++tau) {
            float sum_sq_diff = 0.0f;
            for (int j = 0; j < yin_buffer_size; ++j) {
                float delta = window[j] - window[j + tau];
                sum_sq_diff += delta * delta;
            }
            g_yin_buffer[tau] = sum_sq_diff;
        }

        // 2. Cumulative mean normalized difference function
        float running_sum = 0.0f;
        g_yin_buffer[0] = 1.0f;
        for (int tau = 1; tau < yin_buffer_size; ++tau) {
            running_sum += g_yin_buffer[tau];
            if (running_sum < 1e-7f) {
                 g_yin_buffer[tau] = 1.0f;
            } else {
                 g_yin_buffer[tau] *= (float)tau / running_sum;
            }
        }

        // 3. Find fundamental period
        int period = -1;
        for (int tau = 2; tau < yin_buffer_size; ++tau) {
            if (g_yin_buffer[tau] < yin_threshold) {
                // Find the first local minimum below the threshold
                if (g_yin_buffer[tau] < g_yin_buffer[tau - 1] && g_yin_buffer[tau] < g_yin_buffer[tau + 1]) {
                     period = tau;
                     break;
                }
            }
        }

        // 4. If no period found, find global minimum (fallback)
        if (period == -1) {
            float min_val = 2.0f;
            for (int tau = 2; tau < yin_buffer_size; ++tau) {
                if (g_yin_buffer[tau] < min_val) {
                    min_val = g_yin_buffer[tau];
                    period = tau;
                }
            }
        }

        if (period > 0) {
            pitch_sum += (double)g_sample_rate_hz / period;
        }
        windows_processed++;
    }

    if (windows_processed > 0) {
        g_final_result = pitch_sum / windows_processed;
    } else {
        g_final_result = 0.0;
    }
}

void cleanup() {
    if (g_audio_signal) {
        free(g_audio_signal);
        g_audio_signal = NULL;
    }
    if (g_yin_buffer) {
        free(g_yin_buffer);
        g_yin_buffer = NULL;
    }
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", g_final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
