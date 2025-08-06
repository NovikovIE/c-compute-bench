#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

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
static int events_per_second;
static int pattern_sequence_length;
static long time_window_ms;
static int simulation_duration_s;
static int num_event_types;
static uint64_t time_window_ns;

// Data structures
typedef struct {
    int type;
    int value;
    uint64_t timestamp_ns;
} Event;

typedef struct {
    int active;
    int next_pattern_idx;
    uint64_t start_timestamp_ns;
    long value_accumulator;
} PatternMatcherState;

// Global pointers for data
static Event *event_stream = NULL;
static int *pattern_sequence = NULL;
static PatternMatcherState *active_matchers = NULL;

static long total_events;
static int max_active_matchers;

// Result
static long long final_matched_value = 0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s events_per_second pattern_sequence_length time_window_ms simulation_duration_s num_event_types seed\n", argv[0]);
        exit(1);
    }

    events_per_second = atoi(argv[1]);
    pattern_sequence_length = atoi(argv[2]);
    time_window_ms = atol(argv[3]);
    simulation_duration_s = atoi(argv[4]);
    num_event_types = atoi(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);

    mt_seed(seed);

    total_events = (long)events_per_second * simulation_duration_s;
    time_window_ns = (uint64_t)time_window_ms * 1000000ULL;

    double avg_concurrent_patterns = ((double)events_per_second / (num_event_types > 0 ? num_event_types : 1)) * (time_window_ms / 1000.0);
    max_active_matchers = (int)(avg_concurrent_patterns * 5.0) + 100; // 5x average + buffer for variance

    event_stream = (Event *)malloc(total_events * sizeof(Event));
    pattern_sequence = (int *)malloc(pattern_sequence_length * sizeof(int));
    active_matchers = (PatternMatcherState *)malloc(max_active_matchers * sizeof(PatternMatcherState));

    if (!event_stream || !pattern_sequence || !active_matchers) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate the pattern to match
    for (int i = 0; i < pattern_sequence_length; ++i) {
        pattern_sequence[i] = mt_rand() % num_event_types;
    }

    // Generate the event stream
    uint64_t ns_per_event = 1000000000ULL / events_per_second;
    for (long i = 0; i < total_events; ++i) {
        event_stream[i].type = mt_rand() % num_event_types;
        event_stream[i].value = mt_rand() % 1000;
        event_stream[i].timestamp_ns = i * ns_per_event;
    }
}

void run_computation() {
    memset(active_matchers, 0, max_active_matchers * sizeof(PatternMatcherState));
    final_matched_value = 0;

    for (long i = 0; i < total_events; ++i) {
        Event *current_event = &event_stream[i];

        // Process existing matchers
        for (int j = 0; j < max_active_matchers; ++j) {
            if (active_matchers[j].active) {
                PatternMatcherState *matcher = &active_matchers[j];

                // Check for timeout
                if (current_event->timestamp_ns > matcher->start_timestamp_ns + time_window_ns) {
                    matcher->active = 0; // Deactivate due to timeout
                    continue;
                }

                // Check if the event type matches the next in the sequence
                if (current_event->type == pattern_sequence[matcher->next_pattern_idx]) {
                    matcher->value_accumulator += current_event->value;
                    matcher->next_pattern_idx++;

                    if (matcher->next_pattern_idx == pattern_sequence_length) {
                        final_matched_value += matcher->value_accumulator;
                        matcher->active = 0; // Pattern complete, deactivate
                    }
                } else { // Does not match sequence
                   matcher->active = 0;
                }
            }
        }

        // Check if this event starts a new pattern
        if (current_event->type == pattern_sequence[0]) {
            for (int j = 0; j < max_active_matchers; ++j) {
                if (!active_matchers[j].active) {
                    PatternMatcherState *new_matcher = &active_matchers[j];
                    new_matcher->active = 1;
                    new_matcher->next_pattern_idx = 1;
                    new_matcher->start_timestamp_ns = current_event->timestamp_ns;
                    new_matcher->value_accumulator = current_event->value;

                    // Check if this one-event pattern is the entire pattern
                    if (new_matcher->next_pattern_idx == pattern_sequence_length) {
                        final_matched_value += new_matcher->value_accumulator;
                        new_matcher->active = 0;
                    }
                    break; // Found a slot, stop searching
                }
            }
            // If no free slot, the potential pattern is dropped (as happens in real systems)
        }
    }
}

void cleanup() {
    free(event_stream);
    free(pattern_sequence);
    free(active_matchers);
    event_stream = NULL;
    pattern_sequence = NULL;
    active_matchers = NULL;
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
    printf("%lld\n", final_matched_value);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
