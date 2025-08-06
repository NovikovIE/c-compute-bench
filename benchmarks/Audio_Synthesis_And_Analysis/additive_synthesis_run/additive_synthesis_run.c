#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
// Do Not Modify
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
// --- End of Mersenne Twister ---

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Global Benchmark Data ---
static int DURATION_SECONDS;
static int SAMPLE_RATE_HZ;
static int NUM_HARMONICS;
static long long NUM_SAMPLES;

static double* harmonics_freqs;
static double* harmonics_amps;
static double* output_signal;

static double final_result;

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <duration_seconds> <sample_rate_hz> <num_harmonics> <seed>\n", argv[0]);
        exit(1);
    }

    DURATION_SECONDS = atoi(argv[1]);
    SAMPLE_RATE_HZ = atoi(argv[2]);
    NUM_HARMONICS = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);
    mt_seed(seed);

    if (DURATION_SECONDS <= 0 || SAMPLE_RATE_HZ <= 0 || NUM_HARMONICS <= 0) {
        fprintf(stderr, "FATAL: Parameters must be positive integers.\n");
        exit(1);
    }

    NUM_SAMPLES = (long long)DURATION_SECONDS * SAMPLE_RATE_HZ;

    harmonics_freqs = (double*)malloc(NUM_HARMONICS * sizeof(double));
    harmonics_amps = (double*)malloc(NUM_HARMONICS * sizeof(double));
    output_signal = (double*)malloc(NUM_SAMPLES * sizeof(double));

    if (!harmonics_freqs || !harmonics_amps || !output_signal) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        free(harmonics_freqs);
        free(harmonics_amps);
        free(output_signal);
        exit(1);
    }

    double fundamental_freq = 440.0; // A4 note
    double total_amp = 0.0;

    for (int i = 0; i < NUM_HARMONICS; i++) {
        harmonics_freqs[i] = fundamental_freq * (i + 1);
        double random_factor = (mt_rand() / (double)UINT32_MAX); // [0.0, 1.0]
        // Amplitudes decay for higher harmonics, typical for natural sounds
        harmonics_amps[i] = random_factor / (double)(i + 1);
        total_amp += harmonics_amps[i];
    }
    
    // Normalize amplitudes to prevent clipping (total amp > 1.0)
    if (total_amp > 1.0) {
        for (int i = 0; i < NUM_HARMONICS; i++) {
            harmonics_amps[i] /= total_amp;
        }
    }
}

void run_computation() {
    double checksum = 0.0;
    const double two_pi = 2.0 * M_PI;
    const double inv_sample_rate = 1.0 / (double)SAMPLE_RATE_HZ;

    for (long long i = 0; i < NUM_SAMPLES; i++) {
        double time = (double)i * inv_sample_rate;
        double sample_value = 0.0;
        for (int j = 0; j < NUM_HARMONICS; j++) {
            sample_value += harmonics_amps[j] * sin(two_pi * harmonics_freqs[j] * time);
        }
        output_signal[i] = sample_value;
        checksum += sample_value;
    }
    final_result = checksum;
}

void cleanup() {
    free(harmonics_freqs);
    free(harmonics_amps);
    free(output_signal);
}

// --- Main Function ---
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
