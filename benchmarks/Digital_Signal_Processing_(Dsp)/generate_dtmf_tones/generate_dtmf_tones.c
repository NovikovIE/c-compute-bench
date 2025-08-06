#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
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

// Benchmark parameters
static long num_samples;
static int tone_duration_ms;

// Data structures
static float* output_signal;
static int* tone_sequence;
static int num_tones;
static int samples_per_tone;

// Final result
static double final_result;

// DSP Constants
static const int SAMPLE_RATE = 44100;
static const double PI = 3.14159265358979323846;

// DTMF frequencies (4x4 matrix)
static const float dtmf_low_freq[4] = { 697, 770, 852, 941 };
static const float dtmf_high_freq[4] = { 1209, 1336, 1477, 1633 };

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_samples> <tone_duration_ms> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atol(argv[1]);
    tone_duration_ms = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (num_samples <= 0 || tone_duration_ms <= 0) {
        fprintf(stderr, "Parameters must be positive.\n");
        exit(1);
    }

    // Calculate dependent parameters
    samples_per_tone = (int)((tone_duration_ms / 1000.0) * SAMPLE_RATE);
    if (samples_per_tone == 0) {
        fprintf(stderr, "Tone duration too short for sample rate.\n");
        exit(1);
    }
    num_tones = (int)(num_samples / samples_per_tone) + 1;

    // Allocate memory
    output_signal = (float*)malloc(num_samples * sizeof(float));
    if (output_signal == NULL) {
        perror("Failed to allocate output_signal");
        exit(1);
    }

    tone_sequence = (int*)malloc(num_tones * sizeof(int));
    if (tone_sequence == NULL) {
        perror("Failed to allocate tone_sequence");
        free(output_signal);
        exit(1);
    }

    // Generate the random sequence of tones (0-15 for the 4x4 keypad)
    for (int i = 0; i < num_tones; ++i) {
        tone_sequence[i] = mt_rand() % 16;
    }
}

void run_computation() {
    double checksum = 0.0;
    double two_pi_over_sr = 2.0 * PI / SAMPLE_RATE;

    for (long i = 0; i < num_samples; ++i) {
        int tone_idx_in_seq = i / samples_per_tone;
        int current_tone = tone_sequence[tone_idx_in_seq];

        int row_idx = current_tone / 4;
        int col_idx = current_tone % 4;

        float f1 = dtmf_low_freq[row_idx];
        float f2 = dtmf_high_freq[col_idx];

        // Generate two sine waves and sum them
        // The time t is simply i / SAMPLE_RATE
        double t_i = (double)i;
        float sample = 0.5f * sinf(f1 * two_pi_over_sr * t_i) + 
                       0.5f * sinf(f2 * two_pi_over_sr * t_i);
        
        output_signal[i] = sample;
        checksum += sample;
    }

    final_result = checksum;
}

void cleanup() {
    free(output_signal);
    free(tone_sequence);
    output_signal = NULL;
    tone_sequence = NULL;
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
