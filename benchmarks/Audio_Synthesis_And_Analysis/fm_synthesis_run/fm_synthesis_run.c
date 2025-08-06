#define _POSIX_C_SOURCE 199309L
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// Fix for M_PI not being defined by some compilers when POSIX standards are enforced
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA STRUCTURES ---
typedef struct {
    float frequency;        // Base Hertz
    float amplitude;        // Output amplitude [0.0, 1.0]
    float modulation_index; // Depth of modulation from previous operator
} Operator;

// Global parameters
long g_num_samples;
int g_num_operators;
int g_sample_rate_hz;

// Global data
Operator *g_operators = NULL;
float *g_output_buffer = NULL;
float g_final_result = 0.0f;

// --- BENCHMARK FUNCTIONS ---

// Helper to generate a random float in a given range
float rand_float(float min, float max) {
    return min + ((float)mt_rand() / (float)UINT32_MAX) * (max - min);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <duration_seconds> <sample_rate_hz> <num_operators> <seed>\n", argv[0]);
        exit(1);
    }

    int duration_seconds = atoi(argv[1]);
    g_sample_rate_hz = atoi(argv[2]);
    g_num_operators = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    g_num_samples = (long)duration_seconds * g_sample_rate_hz;
    if (g_num_samples <= 0 || g_num_operators <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters. Samples and operators must be positive.\n");
        exit(1);
    }

    g_operators = (Operator*)malloc(g_num_operators * sizeof(Operator));
    g_output_buffer = (float*)malloc(g_num_samples * sizeof(float));

    if (!g_operators || !g_output_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize operators with random parameters
    for (int i = 0; i < g_num_operators; ++i) {
        // Operator 0 is the carrier, others are modulators
        if (i == 0) {
            g_operators[i].frequency = rand_float(110.0f, 440.0f); // Carrier frequency (A2 to A4)
            g_operators[i].amplitude = 1.0f; // Carrier amplitude is fixed
            g_operators[i].modulation_index = 0.0f; // No input modulation
        } else {
             // Modulator frequencies are often multiples of the carrier
            g_operators[i].frequency = g_operators[0].frequency * rand_float(0.5f, 8.0f);
            g_operators[i].amplitude = rand_float(0.1f, 1.0f);
            g_operators[i].modulation_index = rand_float(1.0f, 20.0f); // Modulation depth
        }
    }
}

void run_computation() {
    const float TWO_PI = 2.0f * (float)M_PI;
    float accumulated_sum = 0.0f;

    for (long i = 0; i < g_num_samples; ++i) {
        float t = (float)i / g_sample_rate_hz;
        float modulator_output = 0.0f;
        
        // The FM chain is processed backwards: the previous operator's output modulates the current one.
        // We calculate from operator N-1 down to 0.
        for (int j = g_num_operators - 1; j >= 0; --j) {
            Operator* op = &g_operators[j];
            float phase_modulation = (j > 0) ? (op->modulation_index * modulator_output) : 0.0f;
            modulator_output = op->amplitude * sinf(TWO_PI * op->frequency * t + phase_modulation);
        }

        g_output_buffer[i] = modulator_output;
        accumulated_sum += modulator_output;
    }
    g_final_result = accumulated_sum;
}

void cleanup() {
    if (g_operators) {
        free(g_operators);
        g_operators = NULL;
    }
    if (g_output_buffer) {
        free(g_output_buffer);
        g_output_buffer = NULL;
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

    // Print result to stdout
    printf("%f\n", g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
