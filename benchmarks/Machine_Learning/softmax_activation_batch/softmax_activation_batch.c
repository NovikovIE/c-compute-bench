#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) --- Do Not Modify ---
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

// Benchmark parameters and data
int NUM_CLASSES;
int BATCH_SIZE;
float *logits;
float *probabilities;
double final_result;

// Function to generate a random float between -5.0 and 5.0
float random_float() {
    return ((float)mt_rand() / (float)UINT32_MAX) * 10.0f - 5.0f;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_classes> <batch_size> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_CLASSES = atoi(argv[1]);
    BATCH_SIZE = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (NUM_CLASSES <= 0 || BATCH_SIZE <= 0) {
        fprintf(stderr, "Error: num_classes and batch_size must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    size_t total_elements = (size_t)BATCH_SIZE * NUM_CLASSES;
    logits = (float *)malloc(total_elements * sizeof(float));
    probabilities = (float *)malloc(total_elements * sizeof(float));

    if (!logits || !probabilities) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    for (size_t i = 0; i < total_elements; ++i) {
        logits[i] = random_float();
    }
}

void run_computation() {
    double accumulator = 0.0;

    for (int i = 0; i < BATCH_SIZE; ++i) {
        float *current_logits = &logits[(size_t)i * NUM_CLASSES];
        float *current_probs = &probabilities[(size_t)i * NUM_CLASSES];

        // 1. Find max logit for numerical stability
        float max_logit = current_logits[0];
        for (int j = 1; j < NUM_CLASSES; ++j) {
            if (current_logits[j] > max_logit) {
                max_logit = current_logits[j];
            }
        }

        // 2. Compute exponentials and sum
        double sum_exp = 0.0;
        for (int j = 0; j < NUM_CLASSES; ++j) {
            // Use double for intermediate sum to maintain precision
            sum_exp += expf(current_logits[j] - max_logit);
        }

        // 3. Normalize to get probabilities
        for (int j = 0; j < NUM_CLASSES; ++j) {
            current_probs[j] = expf(current_logits[j] - max_logit) / (float)sum_exp;
        }

        // 4. Accumulate a value to prevent dead code elimination
        accumulator += current_probs[0];
    }

    final_result = accumulator;
}

void cleanup() {
    free(logits);
    free(probabilities);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
