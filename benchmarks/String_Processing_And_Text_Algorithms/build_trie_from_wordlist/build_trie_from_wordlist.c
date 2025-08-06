#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

#define ALPHABET_SIZE 26

typedef struct TrieNode {
    struct TrieNode* children[ALPHABET_SIZE];
    int is_end_of_word;
} TrieNode;

typedef struct {
    int num_words;
    int avg_word_len;
    uint32_t seed;
    char** word_list;
    TrieNode* trie_root;
    long long total_nodes;
} BenchmarkData;

BenchmarkData g_data;

// Forward declaration for cleanup
void free_trie_recursive(TrieNode* node);

// Helper for run_computation. Its work (malloc, assignments) is part of the computation.
TrieNode* create_new_node() {
    TrieNode* node = (TrieNode*)malloc(sizeof(TrieNode));
    if (!node) {
        perror("Failed to allocate memory for TrieNode");
        exit(EXIT_FAILURE);
    }
    node->is_end_of_word = 0;
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        node->children[i] = NULL;
    }
    g_data.total_nodes++;
    return node;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_words> <average_word_length> <seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    g_data.num_words = atoi(argv[1]);
    g_data.avg_word_len = atoi(argv[2]);
    g_data.seed = (uint32_t)strtoul(argv[3], NULL, 10);
    
    if (g_data.num_words <= 0 || g_data.avg_word_len <= 0) {
        fprintf(stderr, "ERROR: num_words and average_word_length must be positive.\n");
        exit(EXIT_FAILURE);
    }

    g_data.total_nodes = 0;
    g_data.trie_root = NULL;

    mt_seed(g_data.seed);

    g_data.word_list = (char**)malloc(g_data.num_words * sizeof(char*));
    if (!g_data.word_list) {
        perror("Failed to allocate memory for word list");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < g_data.num_words; i++) {
        int length = 1 + (g_data.avg_word_len > 1 ? (mt_rand() % (g_data.avg_word_len * 2 - 1)) : 0);

        g_data.word_list[i] = (char*)malloc((length + 1) * sizeof(char));
        if (!g_data.word_list[i]) {
            perror("Failed to allocate memory for a word");
            for (int j = 0; j < i; j++) {
                free(g_data.word_list[j]);
            }
            free(g_data.word_list);
            exit(EXIT_FAILURE);
        }

        for (int j = 0; j < length; j++) {
            g_data.word_list[i][j] = 'a' + (mt_rand() % ALPHABET_SIZE);
        }
        g_data.word_list[i][length] = '\0';
    }
}

void run_computation() {
    g_data.trie_root = create_new_node(); // Create the root node

    for (int i = 0; i < g_data.num_words; i++) {
        TrieNode* current = g_data.trie_root;
        char* word = g_data.word_list[i];
        for (int j = 0; word[j] != '\0'; j++) {
            int index = word[j] - 'a';
            if (index < 0 || index >= ALPHABET_SIZE) continue; // Safety check
            if (!current->children[index]) {
                current->children[index] = create_new_node();
            }
            current = current->children[index];
        }
        current->is_end_of_word = 1;
    }
}

void cleanup() {
    if (g_data.word_list) {
        for (int i = 0; i < g_data.num_words; i++) {
            free(g_data.word_list[i]);
        }
        free(g_data.word_list);
    }
    
    free_trie_recursive(g_data.trie_root);
}

void free_trie_recursive(TrieNode* node) {
    if (!node) {
        return;
    }
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        free_trie_recursive(node->children[i]);
    }
    free(node);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data.total_nodes);

    fprintf(stderr, "%.6f", time_taken);
    
    cleanup();

    return 0;
}
