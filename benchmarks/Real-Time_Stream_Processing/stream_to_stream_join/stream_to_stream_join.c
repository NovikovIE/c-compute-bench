#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// A constant defining the number of unique keys for events.
// This affects the size of the state and the probability of joins.
#define NUM_KEYS 4096

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

// Represents a single event in a data stream.
typedef struct {
    uint32_t timestamp_ms;
    uint32_t key;
    uint32_t value;
} Event;

// Represents an event from Stream 2 held in the state window.
typedef struct StateNode {
    uint32_t timestamp_ms;
    uint32_t value;
    struct StateNode *next;
} StateNode;

// A global struct to hold all benchmark data and parameters.
struct {
    // Parameters
    long long events_per_second_stream1;
    long long events_per_second_stream2;
    int window_duration_ms;
    int simulation_duration_s;

    // Data
    long num_events1;
    long num_events2;
    Event *stream1;
    Event *stream2;
    StateNode **key_state; // Hash map for storing stream2 events, keyed by `event.key`

    // Result
    unsigned long long join_result_accumulator;
} g_data;

// Comparison function for qsort to sort events by timestamp.
int compare_events(const void *a, const void *b) {
    Event *eventA = (Event *)a;
    Event *eventB = (Event *)b;
    if (eventA->timestamp_ms < eventB->timestamp_ms) return -1;
    if (eventA->timestamp_ms > eventB->timestamp_ms) return 1;
    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <events_per_second_stream1> <events_per_second_stream2> <window_duration_ms> <simulation_duration_s> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.events_per_second_stream1 = atoll(argv[1]);
    g_data.events_per_second_stream2 = atoll(argv[2]);
    g_data.window_duration_ms = atoi(argv[3]);
    g_data.simulation_duration_s = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    g_data.num_events1 = g_data.events_per_second_stream1 * g_data.simulation_duration_s;
    g_data.num_events2 = g_data.events_per_second_stream2 * g_data.simulation_duration_s;

    g_data.stream1 = (Event *)malloc(g_data.num_events1 * sizeof(Event));
    g_data.stream2 = (Event *)malloc(g_data.num_events2 * sizeof(Event));
    g_data.key_state = (StateNode **)calloc(NUM_KEYS, sizeof(StateNode *));

    if (!g_data.stream1 || !g_data.stream2 || !g_data.key_state) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    uint32_t max_timestamp_ms = g_data.simulation_duration_s * 1000;

    for (long i = 0; i < g_data.num_events1; ++i) {
        g_data.stream1[i] = (Event){
            .timestamp_ms = mt_rand() % max_timestamp_ms,
            .key = mt_rand() % NUM_KEYS,
            .value = mt_rand()
        };
    }

    for (long i = 0; i < g_data.num_events2; ++i) {
        g_data.stream2[i] = (Event){
            .timestamp_ms = mt_rand() % max_timestamp_ms,
            .key = mt_rand() % NUM_KEYS,
            .value = mt_rand()
        };
    }

    qsort(g_data.stream1, g_data.num_events1, sizeof(Event), compare_events);
    qsort(g_data.stream2, g_data.num_events2, sizeof(Event), compare_events);

    g_data.join_result_accumulator = 0;
}

void run_computation() {
    long idx1 = 0, idx2 = 0;

    while (idx1 < g_data.num_events1) {
        Event *e1 = &g_data.stream1[idx1];

        // Ingest events from stream2 that occurred before or at the same time as the current stream1 event.
        while (idx2 < g_data.num_events2 && g_data.stream2[idx2].timestamp_ms <= e1->timestamp_ms) {
            Event *e2 = &g_data.stream2[idx2];
            StateNode *new_node = (StateNode *)malloc(sizeof(StateNode));
            if (!new_node) { exit(1); } // Should not happen in benchmark context
            new_node->timestamp_ms = e2->timestamp_ms;
            new_node->value = e2->value;
            new_node->next = g_data.key_state[e2->key];
            g_data.key_state[e2->key] = new_node;
            idx2++;
        }

        // Find joins for e1 and purge expired state for its key.
        StateNode **head_ptr = &g_data.key_state[e1->key];
        StateNode *current = *head_ptr;
        StateNode *prev = NULL;

        while (current != NULL) {
            // Since we prepend to the list, events are ordered by descending timestamp.
            // If we find an expired event, all subsequent events in this list are also expired.
            if (e1->timestamp_ms - current->timestamp_ms > g_data.window_duration_ms) {
                if (prev) {
                    prev->next = NULL;
                } else {
                    *head_ptr = NULL;
                }
                // Free the rest of the list that has expired.
                while (current) {
                    StateNode *next_to_free = current->next;
                    free(current);
                    current = next_to_free;
                }
                break; // Stop processing this list
            }
            
            // Join condition met.
            g_data.join_result_accumulator += e1->value ^ current->value;
            
            prev = current;
            current = current->next;
        }

        idx1++;
    }
}

void cleanup() {
    free(g_data.stream1);
    free(g_data.stream2);

    for (int i = 0; i < NUM_KEYS; ++i) {
        StateNode *current = g_data.key_state[i];
        while (current != NULL) {
            StateNode *next = current->next;
            free(current);
            current = next;
        }
    }
    free(g_data.key_state);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination.
    printf("%llu\n", g_data.join_result_accumulator);

    // Print timing information to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
