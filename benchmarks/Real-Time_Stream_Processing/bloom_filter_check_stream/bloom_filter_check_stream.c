#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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
// --- End of Mersenne Twister ---

// --- Global Benchmark State ---
typedef struct {
    // Parameters
    long events_per_second;
    long num_items_in_filter;
    int hash_functions;

    // Data structures
    unsigned char* bloom_filter; // Bit array
    size_t filter_size_bits;
    size_t filter_size_bytes;
    uint32_t* event_stream; // Stream of items to check

    // Result
    int members_found;
} BenchmarkData;

static BenchmarkData g_data;

// --- Helper Functions - Hashing and Bit Manipulation ---
// Two simple hash functions used for double hashing
uint32_t hash1(uint32_t key) {
    key ^= key >> 16;
    key *= 0x85ebca6b;
    key ^= key >> 13;
    key *= 0xc2b2ae35;
    key ^= key >> 16;
    return key;
}

uint32_t hash2(uint32_t key) {
    key ^= key >> 15;
    key *= 0x1b873593;
    key ^= key >> 12;
    return key;
}

// Set a bit in the bloom filter
void set_bit(size_t i) {
    g_data.bloom_filter[i / 8] |= (1 << (i % 8));
}

// Get a bit from the bloom filter
int get_bit(size_t i) {
    return (g_data.bloom_filter[i / 8] & (1 << (i % 8))) != 0;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s events_per_second num_items_in_filter hash_functions seed\n", argv[0]);
        exit(1);
    }

    g_data.events_per_second = atol(argv[1]);
    g_data.num_items_in_filter = atol(argv[2]);
    g_data.hash_functions = atoi(argv[3]);
    uint32_t seed = (uint32_t)atol(argv[4]);
    mt_seed(seed);

    // 1. Initialize Bloom Filter
    // Rule of thumb: 8 bits per item for a false positive rate of ~2%
    g_data.filter_size_bits = g_data.num_items_in_filter * 8;
    g_data.filter_size_bytes = (g_data.filter_size_bits + 7) / 8;
    g_data.bloom_filter = (unsigned char*)calloc(g_data.filter_size_bytes, sizeof(unsigned char));
    if (!g_data.bloom_filter) {
        fprintf(stderr, "Error: Failed to allocate memory for Bloom filter.\n");
        exit(1);
    }

    // 2. Populate the filter with random items
    uint32_t* items_to_insert = (uint32_t*)malloc(g_data.num_items_in_filter * sizeof(uint32_t));
    if (!items_to_insert) {
        fprintf(stderr, "Error: Failed to allocate memory for items to insert.\n");
        exit(1);
    }

    for (long i = 0; i < g_data.num_items_in_filter; ++i) {
        items_to_insert[i] = mt_rand();
        uint32_t h1 = hash1(items_to_insert[i]);
        uint32_t h2 = hash2(items_to_insert[i]);
        for (int k = 0; k < g_data.hash_functions; ++k) {
            set_bit((h1 + k * h2) % g_data.filter_size_bits);
        }
    }

    // 3. Create the event stream to check against the filter
    g_data.event_stream = (uint32_t*)malloc(g_data.events_per_second * sizeof(uint32_t));
    if (!g_data.event_stream) {
        fprintf(stderr, "Error: Failed to allocate memory for event stream.\n");
        free(items_to_insert);
        exit(1);
    }

    // First half of stream: items known to be in the filter
    long stream_half = g_data.events_per_second / 2;
    for (long i = 0; i < stream_half; ++i) {
        g_data.event_stream[i] = items_to_insert[mt_rand() % g_data.num_items_in_filter];
    }
    
    // Second half of stream: new random items (likely not in the filter)
    for (long i = stream_half; i < g_data.events_per_second; ++i) {
        g_data.event_stream[i] = mt_rand();
    }

    free(items_to_insert);
}

void run_computation() {
    g_data.members_found = 0;
    for (long i = 0; i < g_data.events_per_second; ++i) {
        uint32_t item = g_data.event_stream[i];
        int is_member = 1;
        
        uint32_t h1 = hash1(item);
        uint32_t h2 = hash2(item);

        for (int k = 0; k < g_data.hash_functions; ++k) {
            size_t hash_index = (h1 + k * h2) % g_data.filter_size_bits;
            if (!get_bit(hash_index)) {
                is_member = 0;
                break;
            }
        }
        if (is_member) {
            g_data.members_found++;
        }
    }
}

void cleanup() {
    free(g_data.bloom_filter);
    free(g_data.event_stream);
    g_data.bloom_filter = NULL;
    g_data.event_stream = NULL;
}


// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", g_data.members_found);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
