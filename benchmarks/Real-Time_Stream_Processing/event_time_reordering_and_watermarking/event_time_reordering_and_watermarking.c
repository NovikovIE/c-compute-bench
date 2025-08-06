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

// Data structure for a single event
typedef struct {
    long long event_time_ms;      // Time the event occurred
    long long processing_time_ms; // Time the event arrived for processing
    int data;
} Event;

// Global benchmark data
static long long total_events;
static int events_per_second;
static int max_out_of_orderness_ms;
static int duration_seconds;
static Event* events = NULL;

// Computation-specific data
static Event* pending_events_buffer = NULL;
static long long pending_events_count = 0;

// Final result
static long long final_result;

// Comparison function for qsort to sort events by their processing time
int compare_events(const void* a, const void* b) {
    Event* eventA = (Event*)a;
    Event* eventB = (Event*)b;
    if (eventA->processing_time_ms < eventB->processing_time_ms) return -1;
    if (eventA->processing_time_ms > eventB->processing_time_ms) return 1;
    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <events_per_second> <max_out_of_orderness_ms> <duration_seconds> <seed>\n", argv[0]);
        exit(1);
    }

    events_per_second = atoi(argv[1]);
    max_out_of_orderness_ms = atoi(argv[2]);
    duration_seconds = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    total_events = (long long)events_per_second * duration_seconds;
    if (total_events <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters. Total events must be positive.\n");
        exit(1);
    }

    events = (Event*)malloc(total_events * sizeof(Event));
    pending_events_buffer = (Event*)malloc(total_events * sizeof(Event));
    if (!events || !pending_events_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        free(events);
        free(pending_events_buffer);
        exit(1);
    }

    double event_time_increment = 1000.0 / events_per_second;

    for (long long i = 0; i < total_events; i++) {
        events[i].event_time_ms = (long long)(i * event_time_increment);
        long long delay_ms = mt_rand() % (max_out_of_orderness_ms + 1);
        events[i].processing_time_ms = events[i].event_time_ms + delay_ms;
        events[i].data = mt_rand() % 1000;
    }

    // Sort events by processing_time_ms to simulate out-of-order arrival
    qsort(events, total_events, sizeof(Event), compare_events);
}

void run_computation() {
    const int WINDOW_SIZE_MS = 1000;
    long long max_event_time = 0;
    long long watermark = 0;
    long long next_window_end_time = WINDOW_SIZE_MS;
    long long aggregated_value = 0;

    pending_events_count = 0;

    for (long long i = 0; i < total_events; ++i) {
        Event current_event = events[i];

        pending_events_buffer[pending_events_count++] = current_event;

        if (current_event.event_time_ms > max_event_time) {
            max_event_time = current_event.event_time_ms;
        }

        watermark = max_event_time - max_out_of_orderness_ms;

        while (watermark >= next_window_end_time) {
            long long window_sum = 0;
            long long new_pending_count = 0;

            for (long long j = 0; j < pending_events_count; ++j) {
                if (pending_events_buffer[j].event_time_ms < next_window_end_time) {
                    window_sum += pending_events_buffer[j].data;
                } else {
                    pending_events_buffer[new_pending_count++] = pending_events_buffer[j];
                }
            }
            aggregated_value += window_sum;
            pending_events_count = new_pending_count;

            next_window_end_time += WINDOW_SIZE_MS;
        }
    }

    // Final flush for any windows the final watermark can close
    while (watermark >= next_window_end_time && next_window_end_time <= (duration_seconds * 1000LL)) {
        long long window_sum = 0;
        long long new_pending_count = 0;
        for (long long j = 0; j < pending_events_count; ++j) {
            if (pending_events_buffer[j].event_time_ms < next_window_end_time) {
                window_sum += pending_events_buffer[j].data;
            } else {
                pending_events_buffer[new_pending_count++] = pending_events_buffer[j];
            }
        }
        aggregated_value += window_sum;
        pending_events_count = new_pending_count;
        next_window_end_time += WINDOW_SIZE_MS;
    }

    final_result = aggregated_value;
}

void cleanup() {
    free(events);
    free(pending_events_buffer);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
