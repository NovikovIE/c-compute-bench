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
// This benchmark simulates probing a large endgame database for a game like checkers.
// The 'database' is a large array of pre-calculated results for board states.
// The 'probes' are a sequence of lookups into this database, simulating an AI
// analyzing game positions.

// Parameters
static size_t g_num_probes;
static size_t g_num_database_entries;

// Data structures
static uint32_t* g_database = NULL;
static size_t* g_probe_indices = NULL;

// Result
static uint64_t g_final_result = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    // 1. Argument parsing
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <database_size_mb> <num_probes> <seed>\n", argv[0]);
        exit(1);
    }
    long long database_size_mb = atoll(argv[1]);
    g_num_probes = (size_t)atoll(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    
    // 2. Seed the random number generator
    mt_seed(seed);

    // 3. Calculate database size in entries and allocate memory
    // The database simulates pre-computed endgame states. Each entry is a 32-bit value.
    g_num_database_entries = (size_t)database_size_mb * 1024 * 1024 / sizeof(uint32_t);
    if (g_num_database_entries == 0) {
        fprintf(stderr, "Error: database_size_mb is too small to create any entries.\n");
        exit(1);
    }
    g_database = (uint32_t*)malloc(g_num_database_entries * sizeof(uint32_t));
    if (g_database == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for the database.\n");
        exit(1);
    }

    // 4. Allocate memory for probe indices.
    // These indices represent the sequence of board states to look up.
    g_probe_indices = (size_t*)malloc(g_num_probes * sizeof(size_t));
    if (g_probe_indices == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for probe indices.\n");
        free(g_database);
        exit(1);
    }

    // 5. Populate the database with random values to simulate game state outcomes
    for (size_t i = 0; i < g_num_database_entries; i++) {
        g_database[i] = mt_rand();
    }

    // 6. Generate the random probe sequence
    for (size_t i = 0; i < g_num_probes; i++) {
        // Generate a random index within the bounds of the database
        g_probe_indices[i] = (size_t)mt_rand() % g_num_database_entries;
    }
}

void run_computation() {
    uint64_t accumulator = 0;
    // Iterate through the pre-generated probe indices
    for (size_t i = 0; i < g_num_probes; ++i) {
        // Get the index for the current probe
        size_t probe_idx = g_probe_indices[i];
        
        // Access the database at that index and add to the accumulator
        // This simulates looking up a board state in the endgame database.
        // The addition prevents the compiler from optimizing out the memory access.
        accumulator += g_database[probe_idx];
    }
    g_final_result = accumulator;
}

void cleanup() {
    free(g_database);
    free(g_probe_indices);
    g_database = NULL;
    g_probe_indices = NULL;
}

// --- Main function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%llu\n", (unsigned long long)g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
