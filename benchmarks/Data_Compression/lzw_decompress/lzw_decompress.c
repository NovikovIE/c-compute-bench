#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// DO NOT MODIFY THIS CODE
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
// END OF MT19937

// --- BENCHMARK SPECIFIC CODE ---

// Decompression Dictionary Entry
typedef struct {
    uint8_t* str;
    size_t len;
} DictEntry;

// Global variables to share data between setup, computation, and cleanup
static uint16_t* g_compressed_data = NULL;
static size_t g_compressed_size = 0;
static uint8_t* g_decompressed_data = NULL;
static size_t g_decompressed_size_max = 0;
static int g_max_dictionary_size = 0;

static DictEntry* g_dictionary = NULL;
static uint8_t* g_string_pool = NULL;

static long long g_verification_sum = 0;


// --- LZW COMPRESSION SIMULATION (for setup only) ---

typedef struct CompTrieNode {
    uint16_t code;
    struct CompTrieNode* children[256];
} CompTrieNode;

CompTrieNode* create_comp_node() {
    CompTrieNode* node = (CompTrieNode*)malloc(sizeof(CompTrieNode));
    if (!node) { perror("malloc CompTrieNode failed"); exit(1); }
    memset(node->children, 0, sizeof(node->children));
    node->code = 0;
    return node;
}

void free_comp_trie(CompTrieNode* node) {
    if (!node) return;
    for (int i = 0; i < 256; i++) {
        if (node->children[i]) {
            free_comp_trie(node->children[i]);
        }
    }
    free(node);
}


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <max_dictionary_size> <seed>\n", argv[0]);
        exit(1);
    }

    size_t input_size_mb = atoi(argv[1]);
    g_max_dictionary_size = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    if (g_max_dictionary_size <= 256) {
        fprintf(stderr, "FATAL: max_dictionary_size must be greater than 256.\n");
        exit(1);
    }

    g_decompressed_size_max = input_size_mb * 1024 * 1024;

    // 1. Generate compressible original data
    uint8_t* original_data = (uint8_t*)malloc(g_decompressed_size_max);
    if (!original_data) { perror("malloc original_data failed"); exit(1); }

    const int num_words = 32;
    const int word_len = 8;
    uint8_t words[num_words][word_len];
    for (int i = 0; i < num_words; i++) {
        for (int j = 0; j < word_len; j++) {
            words[i][j] = mt_rand() % 256;
        }
    }
    for (size_t i = 0; i < g_decompressed_size_max; i += word_len) {
        int word_idx = mt_rand() % num_words;
        size_t len_to_copy = (i + word_len > g_decompressed_size_max) ? (g_decompressed_size_max - i) : word_len;
        memcpy(original_data + i, words[word_idx], len_to_copy);
    }

    // 2. Simulate LZW compression to generate a valid compressed stream
    uint16_t* temp_compressed_data = (uint16_t*)malloc(g_decompressed_size_max * sizeof(uint16_t));
    if (!temp_compressed_data) { perror("malloc temp_compressed_data failed"); exit(1); }

    CompTrieNode* root = create_comp_node();
    uint16_t next_code = 256;
    for (int i = 0; i < 256; i++) {
        root->children[i] = create_comp_node();
        root->children[i]->code = i;
    }

    size_t compressed_idx = 0;
    size_t i = 0;
    while (i < g_decompressed_size_max) {
        CompTrieNode* current_node = root;
        size_t match_len = 0;
        uint16_t last_code = original_data[i];

        while (i + match_len < g_decompressed_size_max) {
            uint8_t byte = original_data[i + match_len];
            if (current_node->children[byte]) {
                current_node = current_node->children[byte];
                last_code = current_node->code;
                match_len++;
            } else {
                if (next_code < g_max_dictionary_size) {
                    current_node->children[byte] = create_comp_node();
                    current_node->children[byte]->code = next_code++;
                }
                break;
            }
        }
        temp_compressed_data[compressed_idx++] = last_code;
        i += match_len;
    }

    g_compressed_size = compressed_idx;
    g_compressed_data = (uint16_t*)malloc(g_compressed_size * sizeof(uint16_t));
    if (!g_compressed_data) { perror("malloc g_compressed_data failed"); exit(1); }
    memcpy(g_compressed_data, temp_compressed_data, g_compressed_size * sizeof(uint16_t));

    free(original_data);
    free(temp_compressed_data);
    free_comp_trie(root);

    // 3. Allocate buffers for the actual computation
    g_decompressed_data = (uint8_t*)malloc(g_decompressed_size_max);
    if (!g_decompressed_data) { perror("malloc g_decompressed_data failed"); exit(1); }

    g_dictionary = (DictEntry*)malloc(g_max_dictionary_size * sizeof(DictEntry));
    if (!g_dictionary) { perror("malloc g_dictionary failed"); exit(1); }

    size_t string_pool_size = g_decompressed_size_max * 1.5; // Heuristic safety margin
    g_string_pool = (uint8_t*)malloc(string_pool_size);
    if (!g_string_pool) { perror("malloc g_string_pool failed"); exit(1); }
}

void run_computation() {
    if (g_compressed_size == 0) {
        g_verification_sum = 0;
        return;
    }

    uint8_t* current_string_pool_ptr = g_string_pool;

    // 1. Initialize dictionary with single-character strings
    for (int i = 0; i < 256; i++) {
        g_dictionary[i].str = current_string_pool_ptr;
        g_dictionary[i].len = 1;
        *current_string_pool_ptr++ = (uint8_t)i;
    }
    int next_code = 256;

    uint8_t* decompressed_ptr = g_decompressed_data;
    
    // 2. Main decompression loop
    uint16_t old_code = g_compressed_data[0];
    DictEntry old_entry = g_dictionary[old_code];
    memcpy(decompressed_ptr, old_entry.str, old_entry.len);
    decompressed_ptr += old_entry.len;

    for (size_t i = 1; i < g_compressed_size; i++) {
        uint16_t new_code = g_compressed_data[i];
        DictEntry new_entry;
        uint8_t C;

        if (new_code < next_code) { // Code is in the dictionary
            new_entry = g_dictionary[new_code];
            memcpy(decompressed_ptr, new_entry.str, new_entry.len);
            decompressed_ptr += new_entry.len;
            C = new_entry.str[0];
        } else { // Special case: KwKwK
            new_entry = old_entry;
            C = old_entry.str[0];
            memcpy(decompressed_ptr, new_entry.str, new_entry.len);
            decompressed_ptr += new_entry.len;
            *decompressed_ptr++ = C;
        }

        if (next_code < g_max_dictionary_size) {
            g_dictionary[next_code].len = old_entry.len + 1;
            g_dictionary[next_code].str = current_string_pool_ptr;
            memcpy(current_string_pool_ptr, old_entry.str, old_entry.len);
            current_string_pool_ptr[old_entry.len] = C;
            current_string_pool_ptr += g_dictionary[next_code].len;
            next_code++;
        }

        old_code = new_code;
        old_entry = g_dictionary[old_code];
    }

    // 3. Compute verification sum
    size_t actual_decompressed_size = decompressed_ptr - g_decompressed_data;
    g_verification_sum = 0;
    for (size_t i = 0; i < actual_decompressed_size; i++) {
        g_verification_sum += g_decompressed_data[i];
    }
}

void cleanup() {
    free(g_compressed_data);
    free(g_decompressed_data);
    free(g_dictionary);
    free(g_string_pool);
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
    printf("%lld\n", g_verification_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
