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

// --- Benchmark Data Structures ---
typedef struct {
    uint32_t key;
    uint16_t payload;
} Event;

typedef struct {
    uint32_t value;
} TableRecord;

// Global pointers for benchmark data
Event *event_stream;
TableRecord *lookup_table;
long num_events;
long table_size;

// Global variable to store the final result
volatile int final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <events_per_second> <table_size_records> <seed>\n", argv[0]);
        exit(1);
    }

    long events_per_second = atol(argv[1]);
    table_size = atol(argv[2]);
    uint32_t seed = atoi(argv[3]);

    // For this benchmark, we simulate 1 second's worth of events
    num_events = events_per_second;

    mt_seed(seed);

    // Allocate lookup table
    lookup_table = (TableRecord*)malloc(table_size * sizeof(TableRecord));
    if (lookup_table == NULL) {
        fprintf(stderr, "Failed to allocate memory for lookup_table.\n");
        exit(1);
    }

    // Allocate event stream
    event_stream = (Event*)malloc(num_events * sizeof(Event));
    if (event_stream == NULL) {
        fprintf(stderr, "Failed to allocate memory for event_stream.\n");
        free(lookup_table);
        exit(1);
    }

    // Populate lookup table with random data
    for (long i = 0; i < table_size; i++) {
        lookup_table[i].value = mt_rand();
    }

    // Populate event stream with random keys that point into the lookup table
    for (long i = 0; i < num_events; i++) {
        event_stream[i].key = mt_rand() % table_size;
        event_stream[i].payload = (uint16_t)(mt_rand() % 1000);
    }
}

void run_computation() {
    unsigned long long accumulator = 0;
    for (long i = 0; i < num_events; i++) {
        // Get the event's key
        uint32_t key = event_stream[i].key;

        // Perform the "join": look up the value in the table
        uint32_t table_value = lookup_table[key].value;

        // Perform some stateful aggregation
        // Add the table value and the event's payload to the accumulator
        accumulator += table_value + event_stream[i].payload;
    }
    // Store the result in a volatile global variable to prevent dead code elimination
    final_result = (int)(accumulator & 0xFFFFFFFF);
}

void cleanup() {
    free(event_stream);
    free(lookup_table);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print the final result to stdout to ensure computation is not optimized away
    printf("%d\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
