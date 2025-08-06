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
// --- END MERSENNE TWISTER ---

// Using open addressing (linear probing) for the dictionary to avoid mallocs during computation.
// A key is the combination of a prefix code and the next character.
typedef struct {
    int prefix_code;
    unsigned char new_char;
} DictKey;

// Hash table entry.
typedef struct {
    DictKey key;
    int value_code; // The new code for the string represented by the key.
} DictEntry;

// Global structure to hold all benchmark data
typedef struct {
    unsigned char *input_data;
    size_t input_size;
    int *compressed_data;
    int *compressed_size_ptr; // Pointer to the final compressed size
    size_t max_dictionary_size;
    DictEntry *dictionary_table; 
    size_t dictionary_capacity;  // Actual allocated size of the hash table
    long long final_checksum;
} BenchmarkData;

static BenchmarkData g_data;

// --- LZW HELPER FUNCTIONS ---

// FNV-1a hash variant for dictionary lookups.
static inline size_t hash_key(int prefix, unsigned char c) {
    size_t hash = 0xcbf29ce484222325;
    hash ^= (size_t)prefix;
    hash *= 0x100000001b3;
    hash ^= (size_t)c;
    hash *= 0x100000001b3;
    return hash;
}

// Finds the code for a given string (prefix_code + c) in the dictionary.
// Returns the code if found, -1 otherwise.
static int dictionary_find(int prefix, unsigned char c) {
    size_t index = hash_key(prefix, c) & (g_data.dictionary_capacity - 1);
    while (g_data.dictionary_table[index].value_code != -1) { // -1 marks an empty slot
        if (g_data.dictionary_table[index].key.prefix_code == prefix && g_data.dictionary_table[index].key.new_char == c) {
            return g_data.dictionary_table[index].value_code;
        }
        index = (index + 1) & (g_data.dictionary_capacity - 1); // Linear probe
    }
    return -1;
}

// Adds a new entry (prefix_code + c -> code) to the dictionary.
static void dictionary_insert(int prefix, unsigned char c, int code) {
    size_t index = hash_key(prefix, c) & (g_data.dictionary_capacity - 1);
    while (g_data.dictionary_table[index].value_code != -1) { // Find an empty slot
        index = (index + 1) & (g_data.dictionary_capacity - 1); // Linear probe
    }
    g_data.dictionary_table[index].key.prefix_code = prefix;
    g_data.dictionary_table[index].key.new_char = c;
    g_data.dictionary_table[index].value_code = code;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <max_dictionary_size> <seed>\n", argv[0]);
        exit(1);
    }

    size_t input_size_mb = atol(argv[1]);
    g_data.max_dictionary_size = atol(argv[2]);
    uint32_t seed = atoi(argv[3]);

    g_data.input_size = input_size_mb * 1024 * 1024;

    // Allocate memory for input and compressed output.
    // A safe upper bound for compressed size is the input size (if no compression occurs).
    g_data.input_data = (unsigned char *)malloc(g_data.input_size * sizeof(unsigned char));
    g_data.compressed_data = (int *)malloc(g_data.input_size * sizeof(int));
    g_data.compressed_size_ptr = (int *)malloc(sizeof(int));

    if (!g_data.input_data || !g_data.compressed_data || !g_data.compressed_size_ptr) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Calculate dictionary hash table capacity (next power of 2 * 2 for low load factor)
    g_data.dictionary_capacity = 1;
    while (g_data.dictionary_capacity < g_data.max_dictionary_size) {
        g_data.dictionary_capacity <<= 1;
    }
    g_data.dictionary_capacity *= 2;

    g_data.dictionary_table = (DictEntry*)malloc(g_data.dictionary_capacity * sizeof(DictEntry));
    if (!g_data.dictionary_table) {
        fprintf(stderr, "FATAL: Dictionary memory allocation failed.\n");
        exit(1);
    }

    // Generate semi-compressible random data.
    mt_seed(seed);
    for (size_t i = 0; i < g_data.input_size; ++i) {
        // Use a smaller character set to ensure data is compressible
        g_data.input_data[i] = mt_rand() % 128;
    }
}

void run_computation() {
    if (g_data.input_size == 0) {
        *g_data.compressed_size_ptr = 0;
        g_data.final_checksum = 0;
        return;
    }

    // Reset dictionary table for this run.
    for (size_t i = 0; i < g_data.dictionary_capacity; ++i) {
        g_data.dictionary_table[i].value_code = -1; // Mark all slots as empty.
    }

    int next_code = 256; // Initial dictionary has 256 entries for single bytes.
    int compressed_idx = 0;
    int prefix_code = g_data.input_data[0];

    for (size_t i = 1; i < g_data.input_size; ++i) {
        unsigned char current_char = g_data.input_data[i];
        int found_code = dictionary_find(prefix_code, current_char);

        if (found_code != -1) {
            // String is in the dictionary, extend the current prefix.
            prefix_code = found_code;
        } else {
            // String not found. Output the code for the current prefix.
            g_data.compressed_data[compressed_idx++] = prefix_code;

            // Add the new string (prefix + char) to the dictionary if there's space.
            if ((size_t)next_code < g_data.max_dictionary_size) {
                dictionary_insert(prefix_code, current_char, next_code);
                next_code++;
            }
            
            // Start a new prefix.
            prefix_code = current_char;
        }
    }

    // Output the code for the last remaining prefix.
    g_data.compressed_data[compressed_idx++] = prefix_code;
    *g_data.compressed_size_ptr = compressed_idx;

    // Calculate a checksum of the compressed data to prevent dead code elimination.
    long long checksum = 0;
    for (int i = 0; i < *g_data.compressed_size_ptr; ++i) {
        checksum += g_data.compressed_data[i];
    }
    g_data.final_checksum = checksum;
}

void cleanup() {
    free(g_data.input_data);
    free(g_data.compressed_data);
    free(g_data.compressed_size_ptr);
    free(g_data.dictionary_table);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final checksum to stdout
    printf("%lld\n", g_data.final_checksum);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
