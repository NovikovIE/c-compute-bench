#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA AND PARAMETERS ---
#define MAX_ACTIONS 5
#define MIN_ACTIONS 2

typedef struct {
    int num_actions;
    int children[MAX_ACTIONS]; // Indices of children nodes
    double regret_sum[MAX_ACTIONS];
    double strategy_sum[MAX_ACTIONS];
    double payoff; // For terminal nodes
} Node;

static int NUM_GAME_TREE_NODES;
static int NUM_ITERATIONS;
static Node *game_tree;
static double final_result = 0.0;

// --- BENCHMARK FUNCTIONS ---

// Calculate current strategy based on positive regrets
static void get_strategy(Node* node, double* strategy) {
    double normalizing_sum = 0.0;
    for (int a = 0; a < node->num_actions; ++a) {
        strategy[a] = node->regret_sum[a] > 0 ? node->regret_sum[a] : 0.0;
        normalizing_sum += strategy[a];
    }

    if (normalizing_sum > 0) {
        for (int a = 0; a < node->num_actions; ++a) {
            strategy[a] /= normalizing_sum;
        }
    } else {
        // Default to uniform random strategy if no positive regret
        for (int a = 0; a < node->num_actions; ++a) {
            strategy[a] = 1.0 / node->num_actions;
        }
    }
}

// The core CFR recursive function.
// This is a simplified version for a benchmark, not a fully-featured
// game solver. It captures the computational pattern of tree traversal,
// strategy calculation, and regret updates.
double cfr_recursive(int node_idx) {
    Node* node = &game_tree[node_idx];

    // If terminal node, return its payoff
    if (node->num_actions == 0) {
        return node->payoff;
    }

    double strategy[MAX_ACTIONS];
    get_strategy(node, strategy);

    double child_utils[MAX_ACTIONS];
    double node_util = 0.0;
    for (int a = 0; a < node->num_actions; ++a) {
        // Payoff is from the perspective of the next player, so we negate it.
        child_utils[a] = -cfr_recursive(node->children[a]);
        node_util += strategy[a] * child_utils[a];
    }

    // Update regrets for the current player
    for (int a = 0; a < node->num_actions; ++a) {
        double regret = child_utils[a] - node_util;
        node->regret_sum[a] += regret;
        node->strategy_sum[a] += strategy[a]; // Accumulate for average strategy
    }

    return node_util;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_game_tree_nodes num_iterations seed\n", argv[0]);
        exit(1);
    }

    NUM_GAME_TREE_NODES = atoi(argv[1]);
    NUM_ITERATIONS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    game_tree = (Node *)malloc(NUM_GAME_TREE_NODES * sizeof(Node));
    if (!game_tree) {
        fprintf(stderr, "FATAL: Memory allocation failed for game tree.\n");
        exit(1);
    }

    // Procedurally generate a tree structure
    int next_node_idx = 1; // Node 0 is the root
    for (int i = 0; i < NUM_GAME_TREE_NODES; i++) {
        Node *node = &game_tree[i];
        int is_terminal = (next_node_idx >= NUM_GAME_TREE_NODES);

        if (is_terminal) {
            node->num_actions = 0;
        } else {
            node->num_actions = MIN_ACTIONS + (mt_rand() % (MAX_ACTIONS - MIN_ACTIONS + 1));
        }

        for (int a = 0; a < node->num_actions; a++) {
            if (next_node_idx < NUM_GAME_TREE_NODES) {
                node->children[a] = next_node_idx++;
            } else {
                // Ran out of nodes, so this action leads to a terminal state.
                // We must reduce the action count for this node.
                node->num_actions = a;
                break;
            }
        }
        
        // If node ended up being terminal, assign a random payoff
        if (node->num_actions == 0) {
            node->payoff = (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0; // Payoff in [-1, 1]
        } else {
            node->payoff = 0.0;
        }

        // Initialize regret and strategy sums
        for (int a = 0; a < MAX_ACTIONS; a++) {
            node->regret_sum[a] = 0.0;
            node->strategy_sum[a] = 0.0;
        }
    }
}

void run_computation() {
    double utility_sum = 0.0;
    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        utility_sum += cfr_recursive(0); // Start traversal from the root node (index 0)
    }
    final_result = utility_sum;
}

void cleanup() {
    free(game_tree);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
