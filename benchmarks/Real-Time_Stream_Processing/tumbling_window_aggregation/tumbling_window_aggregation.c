#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER (MT19937) --- Do Not Modify ---
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

// --- BENCHMARK DATA STRUCTURES AND GLOBALS ---
typedef struct {
    uint64_t timestamp_ns;
    uint32_t value;
} Event;

// Parameters
static int simulation_duration_s;
static int events_per_second;
static int window_duration_ms;
static int aggregation_complexity;

// Data
static uint64_t total_events;
static Event* event_stream;

// Result
static uint64_t final_result;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s simulation_duration_s events_per_second window_duration_ms aggregation_complexity seed\n", argv[0]);
        exit(1);
    }

    simulation_duration_s = atoi(argv[1]);
    events_per_second = atoi(argv[2]);
    window_duration_ms = atoi(argv[3]);
    aggregation_complexity = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    total_events = (uint64_t)simulation_duration_s * events_per_second;
    event_stream = (Event*)malloc(total_events * sizeof(Event));
    if (!event_stream) {
        fprintf(stderr, "FATAL: Failed to allocate memory for event stream.\n");
        exit(1);
    }

    uint64_t current_timestamp_ns = 0;
    uint64_t avg_interval_ns = 1000000000ULL / events_per_second;

    for (uint64_t i = 0; i < total_events; i++) {
        event_stream[i].timestamp_ns = current_timestamp_ns;
        event_stream[i].value = mt_rand();
        
        // Simulate jitter by adding a small random value to the interval
        uint64_t jitter_ns = (avg_interval_ns > 10) ? (mt_rand() % (avg_interval_ns / 10)) : 0;
        current_timestamp_ns += avg_interval_ns + jitter_ns;
    }
}

void run_computation() {
    uint64_t window_duration_ns = (uint64_t)window_duration_ms * 1000000ULL;
    if (window_duration_ns == 0) {
        final_result = 0;
        return;
    }

    uint64_t accumulator = 0;
    uint64_t current_window_aggregator = 0;
    uint64_t current_window_end_ns = window_duration_ns;

    for (uint64_t i = 0; i < total_events; i++) {
        // If event is outside the current window, finalize the current window(s)
        while (event_stream[i].timestamp_ns >= current_window_end_ns) {
            accumulator += current_window_aggregator; // Add window result to total
            current_window_aggregator = 0;           // Reset for the new window
            current_window_end_ns += window_duration_ns; // Advance to the next window boundary
        }

        // Process the current event within its window
        uint32_t temp_val = event_stream[i].value;
        for (int k = 0; k < aggregation_complexity; k++) {
            // Arbitrary computation to simulate aggregation work
            temp_val = (temp_val * 19937u + 12345u) ^ (temp_val >> 13);
        }
        current_window_aggregator += temp_val;
    }

    // Add the aggregation result from the last (potentially partial) window
    accumulator += current_window_aggregator;

    final_result = accumulator;
}

void cleanup() {
    free(event_stream);
    event_stream = NULL;
}

// --- MAIN FUNCTION ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (double)(end.tv_sec - start.tv_sec) +
                        (double)(end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%llu\n", (unsigned long long)final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
