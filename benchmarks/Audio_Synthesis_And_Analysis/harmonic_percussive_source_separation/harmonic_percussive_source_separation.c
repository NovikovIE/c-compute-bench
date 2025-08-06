#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

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

// --- BENCHMARK DATA AND GLOBALS ---
typedef struct {
    int num_freq_bins;
    int num_time_frames;
    float margin;
    int kernel_size;

    float* spectrogram;
    float* harmonic_spectrogram;
    float* percussive_spectrogram;
    float* harmonic_mask;
    float* percussive_mask;
    
    float* temp_window; // For median filter

    double final_result;
} BenchmarkData;

BenchmarkData g_data;

int compare_floats(const void* a, const void* b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

// --- BENCHMARK FUNCTIONS ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_freq_bins num_time_frames margin kernel_size seed\n", argv[0]);
        exit(1);
    }

    g_data.num_freq_bins = atoi(argv[1]);
    g_data.num_time_frames = atoi(argv[2]);
    g_data.margin = atof(argv[3]);
    g_data.kernel_size = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    if (g_data.kernel_size % 2 == 0) {
        fprintf(stderr, "Kernel size must be odd.\n");
        exit(1);
    }

    mt_seed(seed);

    size_t spec_size = (size_t)g_data.num_freq_bins * g_data.num_time_frames;
    g_data.spectrogram = (float*)malloc(spec_size * sizeof(float));
    g_data.harmonic_spectrogram = (float*)malloc(spec_size * sizeof(float));
    g_data.percussive_spectrogram = (float*)malloc(spec_size * sizeof(float));
    g_data.harmonic_mask = (float*)malloc(spec_size * sizeof(float));
    g_data.percussive_mask = (float*)malloc(spec_size * sizeof(float));
    g_data.temp_window = (float*)malloc(g_data.kernel_size * sizeof(float));

    if (!g_data.spectrogram || !g_data.harmonic_spectrogram || !g_data.percussive_spectrogram || 
        !g_data.harmonic_mask || !g_data.percussive_mask || !g_data.temp_window) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (size_t i = 0; i < spec_size; i++) {
        g_data.spectrogram[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    g_data.final_result = 0.0;
}

void run_computation() {
    int n_freq = g_data.num_freq_bins;
    int n_time = g_data.num_time_frames;
    int k_size = g_data.kernel_size;
    int k_half = k_size / 2;

    // 1. Horizontal median filtering for harmonic component
    for (int f = 0; f < n_freq; ++f) {
        for (int t = 0; t < n_time; ++t) {
            for (int i = 0; i < k_size; ++i) {
                int sample_t = t - k_half + i;
                if (sample_t < 0) sample_t = 0;
                if (sample_t >= n_time) sample_t = n_time - 1;
                g_data.temp_window[i] = g_data.spectrogram[f * n_time + sample_t];
            }
            qsort(g_data.temp_window, k_size, sizeof(float), compare_floats);
            g_data.harmonic_spectrogram[f * n_time + t] = g_data.temp_window[k_half];
        }
    }

    // 2. Vertical median filtering for percussive component
    for (int t = 0; t < n_time; ++t) {
        for (int f = 0; f < n_freq; ++f) {
            for (int i = 0; i < k_size; ++i) {
                int sample_f = f - k_half + i;
                if (sample_f < 0) sample_f = 0;
                if (sample_f >= n_freq) sample_f = n_freq - 1;
                g_data.temp_window[i] = g_data.spectrogram[sample_f * n_time + t];
            }
            qsort(g_data.temp_window, k_size, sizeof(float), compare_floats);
            g_data.percussive_spectrogram[f * n_time + t] = g_data.temp_window[k_half];
        }
    }

    // 3. Create masks and accumulate result
    size_t spec_size = (size_t)n_freq * n_time;
    double sum = 0.0;
    for (size_t i = 0; i < spec_size; ++i) {
        float h_val = g_data.harmonic_spectrogram[i];
        float p_val = g_data.percussive_spectrogram[i];
        
        if (h_val >= p_val) {
            g_data.harmonic_mask[i] = 1.0f;
            g_data.percussive_mask[i] = 0.0f;
        } else {
            g_data.harmonic_mask[i] = 0.0f;
            g_data.percussive_mask[i] = 1.0f;
        }

        // A slightly more complex Wiener-style filtering for soft masking
        float h_sq = h_val * h_val;
        float p_sq = p_val * p_val;
        float total_sq = h_sq + p_sq + 1e-9f; // Add epsilon to avoid division by zero

        g_data.harmonic_mask[i] = h_sq / total_sq;
        g_data.percussive_mask[i] = p_sq / total_sq;

        sum += g_data.harmonic_mask[i] + g_data.percussive_mask[i];
    }
    g_data.final_result = sum;
}

void cleanup() {
    free(g_data.spectrogram);
    free(g_data.harmonic_spectrogram);
    free(g_data.percussive_spectrogram);
    free(g_data.harmonic_mask);
    free(g_data.percussive_mask);
    free(g_data.temp_window);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}