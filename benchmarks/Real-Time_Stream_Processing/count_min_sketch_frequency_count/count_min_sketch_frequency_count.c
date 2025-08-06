#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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
static uint32_t total_events;
static uint32_t sketch_width;
static uint32_t sketch_depth;

// Data structures
static uint32_t** sketch;      // The Count-Min sketch table
static uint32_t* events;        // The input stream of events
static uint32_t* hash_seeds;    // Seeds for the d hash functions

// Result
static long long final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s total_events sketch_width sketch_depth seed\n", argv[0]);
        exit(1);
    }

    total_events = strtoul(argv[1], NULL, 10);
    sketch_width = strtoul(argv[2], NULL, 10);
    sketch_depth = strtoul(argv[3], NULL, 10);
    uint32_t seed = strtoul(argv[4], NULL, 10);

    mt_seed(seed);

    // Allocate event stream
    events = (uint32_t*)malloc(total_events * sizeof(uint32_t));
    if (!events) { perror("malloc events"); exit(1); }

    // Allocate hash function seeds
    hash_seeds = (uint32_t*)malloc(sketch_depth * sizeof(uint32_t));
    if (!hash_seeds) { perror("malloc hash_seeds"); exit(1); }

    // Allocate the Count-Min sketch (a 2D array)
    sketch = (uint32_t**)malloc(sketch_depth * sizeof(uint32_t*));
    if (!sketch) { perror("malloc sketch"); exit(1); }
    for (uint32_t i = 0; i < sketch_depth; ++i) {
        sketch[i] = (uint32_t*)calloc(sketch_width, sizeof(uint32_t));
        if (!sketch[i]) { perror("calloc sketch[i]"); exit(1); }
    }

    // Generate random events (items in the stream)
    for (uint32_t i = 0; i < total_events; ++i) {
        events[i] = mt_rand();
    }

    // Generate random seeds for the hash functions
    for (uint32_t i = 0; i < sketch_depth; ++i) {
        // Ensure seeds are non-zero and odd for better multiplicative hashing
        hash_seeds[i] = mt_rand() | 1;
    }
}

void run_computation() {
    // Phase 1: Update the sketch with all events
    for (uint32_t i = 0; i < total_events; ++i) {
        uint32_t item = events[i];
        for (uint32_t d = 0; d < sketch_depth; ++d) {
            uint32_t hash = item;
            hash *= hash_seeds[d]; // Simple multiplicative hash
            hash ^= hash >> 16;      // Mix bits
            uint32_t index = hash % sketch_width;
            sketch[d][index]++;
        }
    }

    // Phase 2: Query a subset of items to produce a result and prevent dead code elimination
    final_result = 0;
    uint32_t num_queries = total_events < 1000 ? total_events : 1000;
    for (uint32_t i = 0; i < num_queries; ++i) {
        uint32_t item = events[i];
        uint32_t min_count = -1; // Equivalent to UINT32_MAX

        for (uint32_t d = 0; d < sketch_depth; ++d) {
            uint32_t hash = item;
            hash *= hash_seeds[d];
            hash ^= hash >> 16;
            uint32_t index = hash % sketch_width;
            uint32_t count = sketch[d][index];
            if (count < min_count) {
                min_count = count;
            }
        }
        final_result += min_count;
    }
}

void cleanup() {
    for (uint32_t i = 0; i < sketch_depth; ++i) {
        free(sketch[i]);
    }
    free(sketch);
    free(events);
    free(hash_seeds);
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
