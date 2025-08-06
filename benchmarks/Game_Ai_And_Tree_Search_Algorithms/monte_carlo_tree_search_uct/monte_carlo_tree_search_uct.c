#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Mersenne Twister (MT19937) PRNG (Do Not Modify)
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

// --- Benchmark-specific code ---

// UCT exploration constant
const double UCT_C = 1.41421356237; // sqrt(2)

// Represents a node in the Monte Carlo Search Tree
typedef struct Node {
    int wins;
    int visits;
    int player; // Player who made the move to get to this state. 1 or -1.
    int num_children;
    struct Node* parent;
    struct Node** children;
} Node;

// Global structure to hold benchmark data and parameters
typedef struct {
    int num_simulations;
    int max_playout_depth;
    int branching_factor;
    Node* root;
    int final_result;
} BenchmarkData;

static BenchmarkData g_data;

// Helper to create a new tree node
Node* create_node(Node* parent, int player) {
    Node* node = (Node*)malloc(sizeof(Node));
    if (!node) {
        fprintf(stderr, "FATAL: Memory allocation failed for node.\n");
        exit(1);
    }
    node->wins = 0;
    node->visits = 0;
    node->player = player;
    node->num_children = 0;
    node->parent = parent;
    node->children = (Node**)calloc(g_data.branching_factor, sizeof(Node*));
    if (!node->children) {
        fprintf(stderr, "FATAL: Memory allocation failed for children.\n");
        free(node);
        exit(1);
    }
    return node;
}

// Recursively free tree nodes
void free_tree(Node* node) {
    if (!node) return;
    for (int i = 0; i < node->num_children; i++) {
        free_tree(node->children[i]);
    }
    free(node->children);
    free(node);
}

// Select the best child based on the UCT formula
Node* select_best_child(Node* node) {
    Node* best_child = NULL;
    double max_score = -1.0;

    for (int i = 0; i < node->num_children; i++) {
        Node* child = node->children[i];
        if (child->visits == 0) {
            // Prioritize unvisited children
            return child;
        }
        double uct_score = ((double)child->wins / child->visits) + 
                             UCT_C * sqrt(log((double)node->visits) / child->visits);
        
        if (uct_score > max_score) {
            max_score = uct_score;
            best_child = child;
        }
    }
    return best_child;
}


// --- Core MCTS Functions ---

// 1. Selection & Expansion: Traverse the tree until a leaf or unexpanded node is found
Node* select_and_expand(Node* node) {
    while (1) {
        // If node has potential children that haven't been created yet
        if (node->num_children < g_data.branching_factor) {
            // Expand: create a new child node
            Node* new_child = create_node(node, -node->player);
            node->children[node->num_children++] = new_child;
            return new_child;
        }

        // If node is terminal (all children explored), stop
        if (node->num_children == 0) {
             return node;
        }
        
        // Selection: move to the best child
        node = select_best_child(node);
    }
}

// 2. Simulation (Playout): From a node, simulate a random game to a terminal state
double simulate_playout(Node* node) {
    int current_player = node->player;
    for (int i = 0; i < g_data.max_playout_depth; i++) {
        // In this simple simulation, the outcome is random
        if ((mt_rand() % 2) == 0) {
            // Current player wins
            return (current_player == 1) ? 1.0 : -1.0;
        }
        current_player = -current_player;
    }
    return 0.0; // Draw
}

// 3. Backpropagation: Update statistics from the result of the playout back to the root
void backpropagate(Node* node, double result) {
    while (node != NULL) {
        node->visits++;
        // The result is from player 1's perspective. 
        // If the current node's move was by player -1, the result is inverted for them.
        if (node->player * result > 0) {
            node->wins++;
        }
        node = node->parent;
    }
}

// --- Benchmark Setup, Computation, and Cleanup ---

void setup_benchmark(int argc, char** argv) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_simulations> <max_playout_depth> <branching_factor> <seed>\n", argv[0]);
        exit(1);
    }
    g_data.num_simulations = atoi(argv[1]);
    g_data.max_playout_depth = atoi(argv[2]);
    g_data.branching_factor = atoi(argv[3]);
    unsigned int seed = atoi(argv[4]);

    mt_seed(seed);

    g_data.root = create_node(NULL, -1); // Root represents state before player 1's move
    g_data.root->visits = 1; // Visit root once to avoid log(0)
    g_data.final_result = 0;
}

void run_computation() {
    for (int i = 0; i < g_data.num_simulations; i++) {
        Node* node = select_and_expand(g_data.root);
        double result = simulate_playout(node);
        backpropagate(node, result);
    }

    // After all simulations, determine the best move from the root
    int best_move_index = -1;
    int max_visits = -1;
    for (int i = 0; i < g_data.root->num_children; i++) {
        if (g_data.root->children[i]->visits > max_visits) {
            max_visits = g_data.root->children[i]->visits;
            best_move_index = i;
        }
    }
    g_data.final_result = best_move_index;
}

void cleanup() {
    free_tree(g_data.root);
}

// --- Main Function ---

int main(int argc, char** argv) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
