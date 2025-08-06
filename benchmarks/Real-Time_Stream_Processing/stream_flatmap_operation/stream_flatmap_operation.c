#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// Benchmark-specific data structures
typedef struct {
    uint32_t payload;
    uint32_t user_id;
} InputEvent;

// Global parameters
static long events_per_second;
static int avg_outputs_per_input;

// Global data structures
static InputEvent* input_stream;

// Global result accumulator
static volatile long long final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <events_per_second> <avg_outputs_per_input> <seed>\n", argv[0]);
        exit(1);
    }

    events_per_second = atol(argv[1]);
    avg_outputs_per_input = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    input_stream = (InputEvent*)malloc(events_per_second * sizeof(InputEvent));
    if (input_stream == NULL) {
        fprintf(stderr, "Memory allocation failed for input_stream\n");
        exit(1);
    }

    for (long i = 0; i < events_per_second; i++) {
        input_stream[i].payload = mt_rand();
        input_stream[i].user_id = mt_rand() % 1000; // Simulate a smaller set of users
    }
}

void run_computation() {
    long long accumulator = 0;
    // The number of outputs per input is random, centered on the average.
    // Using `avg * 2` ensures the average is `avg`.
    const int max_outputs = avg_outputs_per_input > 0 ? (2 * avg_outputs_per_input) : 1;

    for (long i = 0; i < events_per_second; i++) {
        InputEvent current_event = input_stream[i];
        int num_outputs = mt_rand() % max_outputs;

        // This is the flatMap operation: one input event generates multiple output values.
        for (int j = 0; j < num_outputs; j++) {
            // Simulate a stateful aggregation or transformation.
            // The calculation is designed to be non-trivial to prevent over-optimization.
            long long output_value = (current_event.payload ^ (long long)i) + (current_event.user_id * (j + 1));
            accumulator += output_value & 0xFFFFFFFF; // Accumulate lower 32 bits
        }
    }
    final_result = accumulator;
}

void cleanup() {
    if (input_stream) {
        free(input_stream);
        input_stream = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%lld\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
