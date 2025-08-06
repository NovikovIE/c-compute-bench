#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Global Variables ---
float* g_input_audio = NULL;
float* g_impulse_response = NULL;
float* g_output_audio = NULL;

int g_input_audio_length;
int g_impulse_response_length;
int g_output_audio_length;

float g_checksum = 0.0f;

// --- Mersenne Twister (Do Not Modify - Include This Verbatim) ---
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

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_audio_length> <impulse_response_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_input_audio_length = atoi(argv[1]);
    g_impulse_response_length = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_input_audio_length <= 0 || g_impulse_response_length <= 0) {
        fprintf(stderr, "Error: Lengths must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    g_output_audio_length = g_input_audio_length + g_impulse_response_length - 1;

    g_input_audio = (float*)malloc((size_t)g_input_audio_length * sizeof(float));
    g_impulse_response = (float*)malloc((size_t)g_impulse_response_length * sizeof(float));
    g_output_audio = (float*)malloc((size_t)g_output_audio_length * sizeof(float));

    if (!g_input_audio || !g_impulse_response || !g_output_audio) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        free(g_input_audio);
        free(g_impulse_response);
        free(g_output_audio);
        exit(1);
    }

    for (int i = 0; i < g_input_audio_length; ++i) {
        g_input_audio[i] = 2.0f * ((float)mt_rand() / (float)UINT32_MAX) - 1.0f;
    }

    for (int i = 0; i < g_impulse_response_length; ++i) {
        g_impulse_response[i] = 2.0f * ((float)mt_rand() / (float)UINT32_MAX) - 1.0f;
    }
}

void run_computation() {
    // Direct convolution: output[n] = sum_{k} (input[n-k] * impulse_response[k])
    for (int n = 0; n < g_output_audio_length; ++n) {
        float sum = 0.0f;
        for (int k = 0; k < g_impulse_response_length; ++k) {
            int input_idx = n - k;
            if (input_idx >= 0 && input_idx < g_input_audio_length) {
                sum += g_input_audio[input_idx] * g_impulse_response[k];
            }
        }
        g_output_audio[n] = sum;
    }

    // Calculate checksum to prevent dead code elimination
    float total_sum = 0.0f;
    for (int i = 0; i < g_output_audio_length; ++i) {
        total_sum += g_output_audio[i];
    }
    g_checksum = total_sum;
}

void cleanup() {
    free(g_input_audio);
    free(g_impulse_response);
    free(g_output_audio);
    g_input_audio = NULL;
    g_impulse_response = NULL;
    g_output_audio = NULL;
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
    printf("%f\n", g_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
