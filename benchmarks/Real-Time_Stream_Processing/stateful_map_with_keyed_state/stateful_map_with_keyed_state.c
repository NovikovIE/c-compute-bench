#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator --- DO NOT MODIFY
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

// --- Benchmark Globals ---
typedef struct
{
    uint32_t key;
    uint32_t value;
} Event;

struct
{
    // Parameters
    long num_events;
    int num_unique_keys;
    int state_size_per_key_bytes;

    // Data structures
    Event *event_stream;
    char **state_table;

    // Result
    unsigned long long final_checksum;
} benchmark_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[])
{
    if (argc != 5)
    {
        fprintf(stderr, "Usage: %s <events_per_second> <num_unique_keys> <state_size_per_key_bytes> <seed>\n", argv[0]);
        exit(1);
    }

    benchmark_data.num_events = atol(argv[1]);
    benchmark_data.num_unique_keys = atoi(argv[2]);
    benchmark_data.state_size_per_key_bytes = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (benchmark_data.state_size_per_key_bytes <= 0 || benchmark_data.state_size_per_key_bytes % sizeof(uint32_t) != 0)
    {
        fprintf(stderr, "FATAL: state_size_per_key_bytes must be a positive multiple of %zu.\n", sizeof(uint32_t));
        exit(1);
    }

    mt_seed(seed);

    // Allocate event stream
    benchmark_data.event_stream = (Event *)malloc(benchmark_data.num_events * sizeof(Event));
    if (!benchmark_data.event_stream)
    {
        perror("FATAL: Memory allocation failed for event stream");
        exit(1);
    }

    // Allocate and initialize state table
    benchmark_data.state_table = (char **)malloc(benchmark_data.num_unique_keys * sizeof(char *));
    if (!benchmark_data.state_table)
    {
        perror("FATAL: Memory allocation failed for state table");
        exit(1);
    }

    for (int i = 0; i < benchmark_data.num_unique_keys; ++i)
    {
        benchmark_data.state_table[i] = (char *)malloc(benchmark_data.state_size_per_key_bytes);
        if (!benchmark_data.state_table[i])
        {
            perror("FATAL: Memory allocation failed for state entry");
            exit(1);
        }
        memset(benchmark_data.state_table[i], 0, benchmark_data.state_size_per_key_bytes);
    }

    // Generate random events
    for (long i = 0; i < benchmark_data.num_events; ++i)
    {
        benchmark_data.event_stream[i].key = mt_rand() % benchmark_data.num_unique_keys;
        benchmark_data.event_stream[i].value = mt_rand();
    }
}

void run_computation()
{
    unsigned long long local_checksum = 0;
    int num_ints_in_state = benchmark_data.state_size_per_key_bytes / sizeof(uint32_t);

    for (long i = 0; i < benchmark_data.num_events; ++i)
    {
        Event current_event = benchmark_data.event_stream[i];
        uint32_t key = current_event.key;
        char *key_state = benchmark_data.state_table[key];
        uint32_t *state_as_ints = (uint32_t *)key_state;

        for (int j = 0; j < num_ints_in_state; ++j)
        {
            // Simple state update function: XOR with event value and position
            state_as_ints[j] ^= (current_event.value + j);
            local_checksum += state_as_ints[j];
        }
    }

    benchmark_data.final_checksum = local_checksum;
}

void cleanup()
{
    for (int i = 0; i < benchmark_data.num_unique_keys; ++i)
    {
        free(benchmark_data.state_table[i]);
    }
    free(benchmark_data.state_table);
    free(benchmark_data.event_stream);
}

// --- Main Function ---
int main(int argc, char *argv[])
{
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%llu\n", benchmark_data.final_checksum);
    
    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
