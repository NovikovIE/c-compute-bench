#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// Benchmark Constants
const double TOTAL_SIMULATION_SECONDS = 5.0;

// Data Structures
typedef struct {
    uint64_t timestamp_us;
    int value;
} Event;

// Global Benchmark Data
Event* event_stream;
long num_events;
uint64_t window_duration_us;
uint64_t slide_duration_us;
uint64_t total_simulation_time_us;
long long final_result = 0; // Use long long to prevent overflow


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <events_per_second> <window_duration_ms> <slide_duration_ms> <seed>\n", argv[0]);
        exit(1);
    }

    long events_per_second = atol(argv[1]);
    long window_duration_ms = atol(argv[2]);
    long slide_duration_ms = atol(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    num_events = (long)(events_per_second * TOTAL_SIMULATION_SECONDS);
    window_duration_us = (uint64_t)window_duration_ms * 1000;
    slide_duration_us = (uint64_t)slide_duration_ms * 1000;
    total_simulation_time_us = (uint64_t)(TOTAL_SIMULATION_SECONDS * 1e6);

    event_stream = (Event*)malloc(num_events * sizeof(Event));
    if (event_stream == NULL) {
        fprintf(stderr, "Failed to allocate memory for event stream.\n");
        exit(1);
    }

    double time_increment_us = 1e6 / (double)events_per_second;
    for (long i = 0; i < num_events; i++) {
        event_stream[i].timestamp_us = (uint64_t)(i * time_increment_us);
        event_stream[i].value = mt_rand() % 1000; // Event values between 0 and 999
    }
}

void run_computation() {
    long long total_aggregated_value = 0;

    for (uint64_t window_start_us = 0;
         window_start_us <= total_simulation_time_us - window_duration_us;
         window_start_us += slide_duration_us) {
        
        uint64_t window_end_us = window_start_us + window_duration_us;
        long long current_window_sum = 0;

        // Naive scan for all events in the current window.
        // This is computationally intensive, simulating a simple but heavy workload.
        for (long i = 0; i < num_events; i++) {
            if (event_stream[i].timestamp_us >= window_start_us &&
                event_stream[i].timestamp_us < window_end_us) {
                current_window_sum += event_stream[i].value;
            }
        }
        total_aggregated_value += current_window_sum;
    }

    final_result = total_aggregated_value;
}


void cleanup() {
    free(event_stream);
    event_stream = NULL;
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
