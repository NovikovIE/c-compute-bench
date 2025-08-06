#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- MERSENNE TWISTER (MT19937) ---
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


// --- BENCHMARK DATA AND FUNCTIONS ---
typedef struct {
    int num_samples;
    int num_transforms_to_run;
    double* input_signals;
    double* output_transforms;
    double final_result;
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_samples> <num_transforms_to_run> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_transforms_to_run = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (g_data.num_samples <= 0 || g_data.num_transforms_to_run <= 0) {
        fprintf(stderr, "FATAL: Parameters must be positive integers.\n");
        exit(1);
    }

    size_t total_elements = (size_t)g_data.num_transforms_to_run * g_data.num_samples;
    g_data.input_signals = (double*)malloc(total_elements * sizeof(double));
    g_data.output_transforms = (double*)malloc(total_elements * sizeof(double));
    if (!g_data.input_signals || !g_data.output_transforms) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (size_t i = 0; i < total_elements; i++) {
        g_data.input_signals[i] = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }
    g_data.final_result = 0.0;
}

void run_computation() {
    int N = g_data.num_samples;
    double pi_over_N = M_PI / N;
    double accumulator = 0.0;

    for (int t = 0; t < g_data.num_transforms_to_run; ++t) {
        double* current_input = g_data.input_signals + (size_t)t * N;
        double* current_output = g_data.output_transforms + (size_t)t * N;

        // Perform DCT-II
        for (int k = 0; k < N; ++k) {
            double sum = 0.0;
            for (int n = 0; n < N; ++n) {
                sum += current_input[n] * cos(pi_over_N * (n + 0.5) * k);
            }
            current_output[k] = sum;
        }

        // Accumulate a value from the output to prevent dead code elimination
        accumulator += current_output[N - 1];
    }

    g_data.final_result = accumulator;
}

void cleanup() {
    free(g_data.input_signals);
    free(g_data.output_transforms);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}