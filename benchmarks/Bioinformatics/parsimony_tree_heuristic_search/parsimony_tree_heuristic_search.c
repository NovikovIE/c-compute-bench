#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) ---
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
// --- END Mersenne Twister ---

// --- Benchmark Globals ---
int num_taxa;
int sequence_length;
int num_iterations;
int num_nodes;

char **sequences; // DNA sequences for the leaf nodes (taxa)

// Tree structure represented by arrays
// Nodes 0 to num_taxa-1 are leaves
// Nodes num_taxa to num_nodes-1 are internal
int *tree_parent;
int *tree_left_child;
int *tree_right_child;
int root_node_index;

int best_parsimony_score;

// --- Helper Functions ---

// Recursive helper for Fitch's algorithm
int post_order_fitch(int node_idx, int site_idx, uint8_t* sets) {
    // Leaf node: set is determined by its sequence
    if (node_idx < num_taxa) {
        char base = sequences[node_idx][site_idx];
        switch(base) {
            case 'A': sets[node_idx] = 1; break; // 0001
            case 'C': sets[node_idx] = 2; break; // 0010
            case 'G': sets[node_idx] = 4; break; // 0100
            case 'T': sets[node_idx] = 8; break; // 1000
        }
        return 0; // No parsimony cost at leaves
    }

    // Internal node: recurse on children
    int score_increment = 0;
    score_increment += post_order_fitch(tree_left_child[node_idx], site_idx, sets);
    score_increment += post_order_fitch(tree_right_child[node_idx], site_idx, sets);

    uint8_t left_set = sets[tree_left_child[node_idx]];
    uint8_t right_set = sets[tree_right_child[node_idx]];

    uint8_t intersection = left_set & right_set;
    if (intersection != 0) {
        sets[node_idx] = intersection;
    } else {
        sets[node_idx] = left_set | right_set;
        score_increment++;
    }
    return score_increment;
}

// Calculates the parsimony score for the entire tree
int calculate_parsimony_score() {
    int total_score = 0;
    uint8_t* site_sets = (uint8_t*)malloc(num_nodes * sizeof(uint8_t));
    if (!site_sets) {
        perror("Failed to allocate memory for site sets");
        exit(1);
    }
    for (int i = 0; i < sequence_length; ++i) {
        total_score += post_order_fitch(root_node_index, i, site_sets);
    }
    free(site_sets);
    return total_score;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_taxa sequence_length num_iterations seed\n", argv[0]);
        exit(1);
    }
    num_taxa = atoi(argv[1]);
    sequence_length = atoi(argv[2]);
    num_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    num_nodes = 2 * num_taxa - 1;
    root_node_index = num_nodes - 1;

    sequences = (char**)malloc(num_taxa * sizeof(char*));
    for (int i = 0; i < num_taxa; ++i) {
        sequences[i] = (char*)malloc((sequence_length + 1) * sizeof(char));
    }

    const char bases[] = "ACGT";
    for (int i = 0; i < num_taxa; ++i) {
        for (int j = 0; j < sequence_length; ++j) {
            sequences[i][j] = bases[mt_rand() % 4];
        }
        sequences[i][sequence_length] = '\0';
    }

    tree_parent = (int*)malloc(num_nodes * sizeof(int));
    tree_left_child = (int*)malloc(num_nodes * sizeof(int));
    tree_right_child = (int*)malloc(num_nodes * sizeof(int));

    for (int i = 0; i < num_nodes; ++i) {
        tree_parent[i] = tree_left_child[i] = tree_right_child[i] = -1;
    }

    int current_internal_node = num_taxa;
    tree_left_child[current_internal_node] = 0;
    tree_right_child[current_internal_node] = 1;
    tree_parent[0] = current_internal_node;
    tree_parent[1] = current_internal_node;
    
    for (int i = 2; i < num_taxa; ++i) {
        int prev_internal_node = current_internal_node;
        current_internal_node++;
        tree_left_child[current_internal_node] = prev_internal_node;
        tree_right_child[current_internal_node] = i;
        tree_parent[prev_internal_node] = current_internal_node;
        tree_parent[i] = current_internal_node;
    }
    tree_parent[root_node_index] = -1;
}

void run_computation() {
    best_parsimony_score = calculate_parsimony_score();

    for (int i = 0; i < num_iterations; ++i) {
        int leaf1_idx = mt_rand() % num_taxa;
        int leaf2_idx;
        do {
            leaf2_idx = mt_rand() % num_taxa;
        } while (leaf1_idx == leaf2_idx);

        int p1 = tree_parent[leaf1_idx];
        int p2 = tree_parent[leaf2_idx];
        int p1_is_left = (tree_left_child[p1] == leaf1_idx);
        int p2_is_left = (tree_left_child[p2] == leaf2_idx);

        if (p1 == p2) continue;

        tree_parent[leaf1_idx] = p2;
        tree_parent[leaf2_idx] = p1;
        if (p1_is_left) tree_left_child[p1] = leaf2_idx;
        else tree_right_child[p1] = leaf2_idx;
        if (p2_is_left) tree_left_child[p2] = leaf1_idx;
        else tree_right_child[p2] = leaf1_idx;

        int current_score = calculate_parsimony_score();
        
        if (current_score < best_parsimony_score) {
            best_parsimony_score = current_score;
        } else {
            tree_parent[leaf1_idx] = p1;
            tree_parent[leaf2_idx] = p2;
            if (p1_is_left) tree_left_child[p1] = leaf1_idx;
            else tree_right_child[p1] = leaf1_idx;
            if (p2_is_left) tree_left_child[p2] = leaf2_idx;
            else tree_right_child[p2] = leaf2_idx;
        }
    }
}

void cleanup() {
    for (int i = 0; i < num_taxa; i++) {
        free(sequences[i]);
    }
    free(sequences);
    free(tree_parent);
    free(tree_left_child);
    free(tree_right_child);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", best_parsimony_score);

    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
