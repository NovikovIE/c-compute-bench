#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Specific Code ---

#define ALPHABET_SIZE 26
#define ROOT_NODE 0

typedef struct {
    int next_states[ALPHABET_SIZE];
    int patterns_ending_here;
    int cumulative_matches;
    int failure_link;
} AhoNode;

typedef struct {
    // Parameters
    long text_length;
    int num_patterns;
    int average_pattern_length;

    // Data
    char *text;
    char **patterns;
    AhoNode *trie;
    int trie_nodes_count;
    int *bfs_queue;

    // Result
    long long total_matches;
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s text_length num_patterns average_pattern_length seed\n", argv[0]);
        exit(1);
    }

    g_data.text_length = atol(argv[1]);
    g_data.num_patterns = atoi(argv[2]);
    g_data.average_pattern_length = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    mt_seed(seed);

    // Allocate memory
    g_data.text = (char*) malloc(g_data.text_length + 1);
    g_data.patterns = (char**) malloc(g_data.num_patterns * sizeof(char*));
    if (!g_data.text || !g_data.patterns) {
        fprintf(stderr, "Memory allocation failed for text/patterns.\n");
        exit(1);
    }

    // Generate random text
    for (long i = 0; i < g_data.text_length; ++i) {
        g_data.text[i] = 'a' + (mt_rand() % ALPHABET_SIZE);
    }
    g_data.text[g_data.text_length] = '\0';

    // Generate random patterns
    size_t total_pattern_chars = 0;
    for (int i = 0; i < g_data.num_patterns; ++i) {
        int length = g_data.average_pattern_length / 2 + (mt_rand() % g_data.average_pattern_length);
        if (length == 0) length = 1;
        g_data.patterns[i] = (char*) malloc(length + 1);
        if (!g_data.patterns[i]) {
            fprintf(stderr, "Memory allocation failed for a pattern.\n");
            exit(1);
        }
        for (int j = 0; j < length; ++j) {
            g_data.patterns[i][j] = 'a' + (mt_rand() % ALPHABET_SIZE);
        }
        g_data.patterns[i][length] = '\0';
        total_pattern_chars += length;
    }

    // --- Build Aho-Corasick State Machine ---
    size_t max_trie_nodes = total_pattern_chars + 1;
    g_data.trie = (AhoNode*) malloc(max_trie_nodes * sizeof(AhoNode));
    g_data.bfs_queue = (int*) malloc(max_trie_nodes * sizeof(int));
    if (!g_data.trie || !g_data.bfs_queue) {
        fprintf(stderr, "Memory allocation failed for trie/queue.\n");
        exit(1);
    }

    // Initialize root node
    AhoNode* root = &g_data.trie[ROOT_NODE];
    memset(root->next_states, -1, sizeof(root->next_states));
    root->patterns_ending_here = 0;
    root->cumulative_matches = 0;
    root->failure_link = ROOT_NODE;
    g_data.trie_nodes_count = 1;

    // 1. Build the basic trie
    for (int i = 0; i < g_data.num_patterns; ++i) {
        int current_node_idx = ROOT_NODE;
        for (char *p = g_data.patterns[i]; *p; ++p) {
            int char_idx = *p - 'a';
            if (g_data.trie[current_node_idx].next_states[char_idx] == -1) {
                int new_node_idx = g_data.trie_nodes_count++;
                g_data.trie[current_node_idx].next_states[char_idx] = new_node_idx;
                AhoNode* new_node = &g_data.trie[new_node_idx];
                memset(new_node->next_states, -1, sizeof(new_node->next_states));
                new_node->patterns_ending_here = 0;
                new_node->cumulative_matches = 0;
                new_node->failure_link = ROOT_NODE; 
            }
            current_node_idx = g_data.trie[current_node_idx].next_states[char_idx];
        }
        g_data.trie[current_node_idx].patterns_ending_here++;
    }

    // 2. Build failure links using BFS and precompute cumulative matches
    for (int i=0; i < g_data.trie_nodes_count; ++i) {
       g_data.trie[i].cumulative_matches = g_data.trie[i].patterns_ending_here;
    }

    int q_head = 0, q_tail = 0;
    for (int i = 0; i < ALPHABET_SIZE; ++i) {
        if (root->next_states[i] != -1) {
            g_data.trie[root->next_states[i]].failure_link = ROOT_NODE;
            g_data.bfs_queue[q_tail++] = root->next_states[i];
        }
    }

    while (q_head < q_tail) {
        int u = g_data.bfs_queue[q_head++];
        
        for (int i = 0; i < ALPHABET_SIZE; ++i) {
            int v = g_data.trie[u].next_states[i];
            if (v == -1) continue;

            int f = g_data.trie[u].failure_link;
            while (f != ROOT_NODE && g_data.trie[f].next_states[i] == -1) {
                f = g_data.trie[f].failure_link;
            }
            int failure_dest = (g_data.trie[f].next_states[i] != -1) ? g_data.trie[f].next_states[i] : ROOT_NODE;
            g_data.trie[v].failure_link = failure_dest;
            g_data.trie[v].cumulative_matches += g_data.trie[failure_dest].cumulative_matches;

            g_data.bfs_queue[q_tail++] = v;
        }
    }
}

void run_computation() {
    g_data.total_matches = 0;
    int current_state = ROOT_NODE;

    for (long i = 0; i < g_data.text_length; ++i) {
        int char_idx = g_data.text[i] - 'a';
        
        while (current_state != ROOT_NODE && g_data.trie[current_state].next_states[char_idx] == -1) {
            current_state = g_data.trie[current_state].failure_link;
        }

        if (g_data.trie[current_state].next_states[char_idx] != -1) {
            current_state = g_data.trie[current_state].next_states[char_idx];
        } 

        g_data.total_matches += g_data.trie[current_state].cumulative_matches;
    }
}

void cleanup() {
    for (int i = 0; i < g_data.num_patterns; ++i) {
        free(g_data.patterns[i]);
    }
    free(g_data.patterns);
    free(g_data.text);
    free(g_data.trie);
    free(g_data.bfs_queue);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data.total_matches);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
