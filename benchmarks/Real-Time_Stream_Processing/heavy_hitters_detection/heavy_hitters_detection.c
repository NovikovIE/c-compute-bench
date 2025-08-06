#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// --- Benchmark-specific data structures and globals ---

// Hash entry for counting frequencies
typedef struct {
    uint32_t key;
    uint32_t count;
} HashEntry;

// Structure for holding top K results for sorting
typedef struct {
    uint32_t id;
    uint32_t count;
} ResultEntry;

// Global struct to hold all benchmark data
typedef struct {
    long total_events;
    int top_k;
    uint32_t seed;

    uint32_t *event_stream;
    HashEntry *count_map;
    size_t count_map_size;

    long long final_result; // For stdout, prevents dead code elimination
} BenchmarkData;

static BenchmarkData* g_data;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
int compare_results(const void *a, const void *b);

// --- Benchmark Implementation ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <total_events> <top_k_elements_to_find> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("Failed to allocate memory for BenchmarkData");
        exit(1);
    }

    g_data->total_events = atol(argv[1]);
    g_data->top_k = atoi(argv[2]);
    g_data->seed = atoi(argv[3]);

    mt_seed(g_data->seed);

    // Allocate the event stream
    g_data->event_stream = (uint32_t*)malloc(g_data->total_events * sizeof(uint32_t));
    if (!g_data->event_stream) {
        perror("Failed to allocate event stream");
        exit(1);
    }

    // Use a fixed vocabulary size for events and size the hash map accordingly
    // Using a prime number for size helps reduce collisions
    const size_t VOCABULARY_SIZE = 250000;
    g_data->count_map_size = 500009;
    g_data->count_map = (HashEntry*)calloc(g_data->count_map_size, sizeof(HashEntry));
    if (!g_data->count_map) {
        perror("Failed to allocate count map");
        exit(1);
    }
    
    // Generate a skewed data stream to create 'heavy hitters'
    // 80% of events will be from the top 1% of the vocabulary
    const uint32_t popular_item_threshold = VOCABULARY_SIZE * 0.01;
    for (long i = 0; i < g_data->total_events; ++i) {
        uint32_t item_id;
        if ((mt_rand() % 100) < 80) {
            // Generate a 'popular' item
            item_id = (mt_rand() % popular_item_threshold) + 1;
        } else {
            // Generate a 'rare' item
            item_id = (mt_rand() % (VOCABULARY_SIZE - popular_item_threshold)) + popular_item_threshold + 1;
        }
        g_data->event_stream[i] = item_id;
    }
    g_data->final_result = 0;
}

void run_computation() {
    // Part 1: Process stream and count frequencies using a hash map
    for (long i = 0; i < g_data->total_events; ++i) {
        uint32_t key = g_data->event_stream[i];
        size_t index = key % g_data->count_map_size;
        
        // Linear probing to find spot or existing entry
        while (g_data->count_map[index].key != 0 && g_data->count_map[index].key != key) {
            index = (index + 1) % g_data->count_map_size;
        }

        if (g_data->count_map[index].key == 0) {
            g_data->count_map[index].key = key;
            g_data->count_map[index].count = 1;
        } else {
            g_data->count_map[index].count++;
        }
    }

    // Part 2: Extract unique items and their counts
    size_t unique_item_count = 0;
    for (size_t i = 0; i < g_data->count_map_size; ++i) {
        if (g_data->count_map[i].key != 0) {
            unique_item_count++;
        }
    }

    ResultEntry *results = (ResultEntry*)malloc(unique_item_count * sizeof(ResultEntry));
    if (!results) {
        perror("Failed to allocate results array");
        exit(1);
    }

    size_t result_idx = 0;
    for (size_t i = 0; i < g_data->count_map_size; ++i) {
        if (g_data->count_map[i].key != 0) {
            results[result_idx].id = g_data->count_map[i].key;
            results[result_idx].count = g_data->count_map[i].count;
            result_idx++;
        }
    }

    // Part 3: Sort the results to find the top K heavy hitters
    qsort(results, unique_item_count, sizeof(ResultEntry), compare_results);

    // Part 4: Calculate a final result checksum from the top K elements
    long long checksum = 0;
    int limit = g_data->top_k < unique_item_count ? g_data->top_k : unique_item_count;
    for (int i = 0; i < limit; ++i) {
        checksum += (long long)results[i].id * results[i].count;
    }
    g_data->final_result = checksum;

    // Clean up temporary allocations within this function
    free(results);
}

void cleanup() {
    free(g_data->event_stream);
    free(g_data->count_map);
    free(g_data);
}

// Comparison function for qsort (descending order by count)
int compare_results(const void *a, const void *b) {
    ResultEntry *entryA = (ResultEntry *)a;
    ResultEntry *entryB = (ResultEntry *)b;
    if (entryB->count > entryA->count) return 1;
    if (entryB->count < entryA->count) return -1;
    return 0; // Or sort by ID for tie-breaking: return entryA->id - entryB->id;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    long long result = g_data->final_result;

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
