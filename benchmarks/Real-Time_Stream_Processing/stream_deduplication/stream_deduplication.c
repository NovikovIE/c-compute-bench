#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of MT19937 ---

// Event struct representing a data point in the stream
typedef struct {
    uint32_t key;          // The key to deduplicate on
    uint64_t timestamp_ns; // Timestamp of the event in nanoseconds
    int value;             // A value to aggregate
} Event;

// Hash table entry for tracking last seen keys
typedef struct {
    uint32_t key;
    uint64_t last_seen_ns;
} HashTableEntry;

// --- Global Benchmark State ---
struct {
    long events_per_second;
    long num_unique_keys;
    uint64_t time_window_ns;

    long long total_events;
    Event* events;

    HashTableEntry* seen_keys_table;
    size_t table_size; // Must be a power of 2
    size_t table_mask; // table_size - 1

    long long final_result; // Accumulated result to prevent dead code elimination
} g_state;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
size_t next_power_of_2(size_t n);


// Finds the next power of 2 >= n. Hleper for hash table sizing.
size_t next_power_of_2(size_t n) {
    if (n == 0) return 1;
    size_t p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <events_per_second> <time_window_ms> <num_unique_keys> <seed>\n", argv[0]);
        exit(1);
    }

    g_state.events_per_second = atol(argv[1]);
    long time_window_ms = atol(argv[2]);
    g_state.num_unique_keys = atol(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    g_state.time_window_ns = (uint64_t)time_window_ms * 1000000;

    // Simulate a stream of 5 seconds for a substantial workload
    int simulation_duration_s = 5;
    g_state.total_events = g_state.events_per_second * simulation_duration_s;

    // Allocate memory for the event stream
    g_state.events = (Event*)malloc(g_state.total_events * sizeof(Event));
    if (!g_state.events) {
        fprintf(stderr, "Failed to allocate memory for events\n");
        exit(1);
    }

    // Setup the hash table for deduplication
    g_state.table_size = next_power_of_2(g_state.num_unique_keys);
    g_state.table_mask = g_state.table_size - 1;
    g_state.seen_keys_table = (HashTableEntry*)calloc(g_state.table_size, sizeof(HashTableEntry));
    if (!g_state.seen_keys_table) {
        fprintf(stderr, "Failed to allocate memory for hash table\n");
        free(g_state.events);
        exit(1);
    }

    // Generate event data
    uint64_t ns_per_event = 1000000000 / g_state.events_per_second;
    for (long long i = 0; i < g_state.total_events; ++i) {
        g_state.events[i].key = (mt_rand() % g_state.num_unique_keys) + 1; // Keys are > 0
        g_state.events[i].timestamp_ns = i * ns_per_event;
        g_state.events[i].value = mt_rand() % 100;
    }
}

void run_computation() {
    long long processed_value_sum = 0;

    for (long long i = 0; i < g_state.total_events; ++i) {
        Event* current_event = &g_state.events[i];
        uint32_t key = current_event->key;
        
        size_t hash_index = key & g_state.table_mask;
        int is_duplicate = 0;

        // Linear probing to find the key or an empty slot
        while (1) {
            HashTableEntry* entry = &g_state.seen_keys_table[hash_index];

            if (entry->key == 0) { // Empty slot, key is new
                entry->key = key;
                entry->last_seen_ns = current_event->timestamp_ns;
                break;
            }
            
            if (entry->key == key) { // Found the key
                if (current_event->timestamp_ns - entry->last_seen_ns <= g_state.time_window_ns) {
                    is_duplicate = 1;
                } else {
                    // Stale entry, update timestamp
                    entry->last_seen_ns = current_event->timestamp_ns;
                }
                break;
            }

            // Collision, move to the next slot
            hash_index = (hash_index + 1) & g_state.table_mask;
        }

        if (!is_duplicate) {
            processed_value_sum += current_event->value;
        }
    }

    g_state.final_result = processed_value_sum;
}

void cleanup() {
    free(g_state.events);
    free(g_state.seen_keys_table);
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
    printf("%lld\n", g_state.final_result);

    // Print the timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
