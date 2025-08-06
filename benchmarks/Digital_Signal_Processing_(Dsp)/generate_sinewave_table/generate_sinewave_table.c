#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
// --- END MERSENNE TWISTER ---

// Benchmark parameters
int num_samples;
int num_harmonics;

// Data structures
float *sine_table;
float *amplitudes;
float *frequencies;

// Result accumulator
double final_result = 0.0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_samples> <num_harmonics> <seed>\n", argv[0]);
        exit(1);
    }

    num_samples = atoi(argv[1]);
    num_harmonics = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    sine_table = (float *)malloc(num_samples * sizeof(float));
    if (!sine_table) {
        perror("Failed to allocate sine_table");
        exit(1);
    }

    amplitudes = (float *)malloc(num_harmonics * sizeof(float));
    if (!amplitudes) {
        perror("Failed to allocate amplitudes");
        exit(1);
    }

    frequencies = (float *)malloc(num_harmonics * sizeof(float));
    if (!frequencies) {
        perror("Failed to allocate frequencies");
        exit(1);
    }

    // Initialize harmonic properties with random values
    for (int i = 0; i < num_harmonics; ++i) {
        amplitudes[i] = (float)mt_rand() / (float)UINT32_MAX;         // Amplitude between 0.0 and 1.0
        frequencies[i] = ((float)mt_rand() / (float)UINT32_MAX) * 10.0f + 1.0f; // Freq multiplier between 1.0 and 11.0
    }
}

void run_computation() {
    // Initialize sine table to zero
    for (int i = 0; i < num_samples; ++i) {
        sine_table[i] = 0.0f;
    }

    // Generate sine wave by summing harmonics
    for (int h = 0; h < num_harmonics; ++h) {
        for (int i = 0; i < num_samples; ++i) {
            float angle = 2.0f * (float)M_PI * frequencies[h] * (float)i / (float)num_samples;
            sine_table[i] += amplitudes[h] * sinf(angle);
        }
    }

    // Calculate a final checksum to prevent dead code elimination
    double sum = 0.0;
    for (int i = 0; i < num_samples; ++i) {
        sum += sine_table[i];
    }
    final_result = sum;
}

void cleanup() {
    free(sine_table);
    free(amplitudes);
    free(frequencies);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
