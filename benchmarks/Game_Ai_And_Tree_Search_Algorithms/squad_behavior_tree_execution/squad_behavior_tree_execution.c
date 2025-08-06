/*
 * squad_behavior_tree_execution: Game AI and Tree Search Algorithms Benchmark
 *
 * Description:
 * This benchmark simulates the execution of Behavior Trees (BTs) for a squad of AI agents.
 * Behavior Trees are a popular AI model in game development for creating complex and modular
 * behaviors. The program constructs a unique BT for each agent, with a structure determined
 * by the 'tree_depth' parameter. The simulation runs for a specified number of 'ticks'.
 * In each tick, every agent's BT is evaluated. The traversal of the tree depends on the
 * agent's current state (e.g., health) and the type of nodes in the tree (Selectors, 
 * Sequences, and Actions). Action nodes perform simple computations and modify the agent's
 * state, influencing subsequent decisions. The total computational work is a function of
 * the number of agents, the simulation duration, and the complexity of the BTs.
 *
 * The benchmark heavily exercises conditional branching, pointer/index chasing (tree traversal),
 * and integer arithmetic. It is designed to be CPU-bound, with memory access patterns
 * determined by the tree traversal logic.
 */

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

// --- Benchmark Specific Data Structures ---

#define MAX_CHILDREN 4
#define LOW_HEALTH_THRESHOLD 30
#define MAX_HEALTH 100

typedef enum {
    SUCCESS,
    FAILURE
} BTNodeStatus;

typedef enum {
    SELECTOR,       // Fallback node: tries children until one succeeds
    SEQUENCE,       // Runs children in order until one fails
    CONDITION_HEALTH_LOW, // Leaf: succeeds if health is low
    ACTION_ATTACK,  // Leaf: performs an attack action
    ACTION_FLEE,    // Leaf: performs a flee action
    ACTION_HEAL     // Leaf: performs a heal action
} BTNodeType;

typedef struct {
    BTNodeType type;
    int num_children;
    int child_indices[MAX_CHILDREN];
} BTNode;

typedef struct {
    int id;
    int health;
    int x, y; // Position
    int bt_root_idx;
} Agent;

// --- Global Variables ---

int NUM_AGENTS;
int TREE_DEPTH;
int TICKS_TO_SIMULATE;

Agent *agents;
BTNode *node_pool;
size_t node_pool_size;
size_t next_node_index;

long long final_result_accumulator;

// --- Forward declarations ---
int generate_tree_recursive(int depth);
BTNodeStatus execute_node(Agent *agent, int node_idx);

// --- Benchmark Functions ---

unsigned long long int_pow(int base, int exp) {
    unsigned long long res = 1;
    for (int i = 0; i < exp; ++i) res *= base;
    return res;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_agents> <tree_depth> <ticks_to_simulate> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_AGENTS = atoi(argv[1]);
    TREE_DEPTH = atoi(argv[2]);
    TICKS_TO_SIMULATE = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory
    agents = (Agent*)malloc(NUM_AGENTS * sizeof(Agent));
    if (!agents) { perror("Failed to allocate agents"); exit(1); }
    
    // Calculate max possible nodes to pre-allocate a single large pool.
    // Sum of geometric series for a full K-ary tree: (K^(D+1) - 1) / (K-1)
    unsigned long long max_nodes_per_tree = (int_pow(MAX_CHILDREN, TREE_DEPTH + 1) - 1) / (MAX_CHILDREN - 1);
    node_pool_size = (size_t)NUM_AGENTS * max_nodes_per_tree;
    node_pool = (BTNode*)malloc(node_pool_size * sizeof(BTNode));
    if (!node_pool) { perror("Failed to allocate node pool"); exit(1); }

    next_node_index = 0;

    // Initialize agents and their behavior trees
    for (int i = 0; i < NUM_AGENTS; ++i) {
        agents[i].id = i;
        agents[i].health = 50 + (mt_rand() % 51); // 50-100
        agents[i].x = mt_rand() % 1000;
        agents[i].y = mt_rand() % 1000;
        agents[i].bt_root_idx = generate_tree_recursive(0);
    }
}

int generate_tree_recursive(int depth) {
    if (next_node_index >= node_pool_size) {
        fprintf(stderr, "FATAL: Node pool exhausted during tree generation.\n");
        exit(1);
    }

    int current_node_idx = next_node_index++;
    BTNode* node = &node_pool[current_node_idx];

    if (depth >= TREE_DEPTH) { // Leaf node
        node->num_children = 0;
        switch (mt_rand() % 4) {
            case 0: node->type = CONDITION_HEALTH_LOW; break;
            case 1: node->type = ACTION_ATTACK; break;
            case 2: node->type = ACTION_FLEE; break;
            case 3: node->type = ACTION_HEAL; break;
        }
    } else { // Internal node
        node->type = (mt_rand() % 2 == 0) ? SELECTOR : SEQUENCE;
        node->num_children = 2 + (mt_rand() % (MAX_CHILDREN - 1)); // 2 to MAX_CHILDREN children
        for (int i = 0; i < node->num_children; ++i) {
            node->child_indices[i] = generate_tree_recursive(depth + 1);
        }
    }
    return current_node_idx;
}

BTNodeStatus execute_node(Agent *agent, int node_idx) {
    BTNode* node = &node_pool[node_idx];

    switch(node->type) {
        case SELECTOR:
            for (int i = 0; i < node->num_children; ++i) {
                if (execute_node(agent, node->child_indices[i]) == SUCCESS) {
                    return SUCCESS;
                }
            }
            return FAILURE;
        case SEQUENCE:
            for (int i = 0; i < node->num_children; ++i) {
                if (execute_node(agent, node->child_indices[i]) == FAILURE) {
                    return FAILURE;
                }
            }
            return SUCCESS;
        case CONDITION_HEALTH_LOW:
            return (agent->health < LOW_HEALTH_THRESHOLD) ? SUCCESS : FAILURE;
        case ACTION_ATTACK:
            // Simple computation to simulate action
            agent->health -= 2;
            final_result_accumulator += (agent->id ^ agent->x);
            return SUCCESS;
        case ACTION_FLEE:
            agent->x += (agent->id % 5) - 2;
            agent->y -= (agent->id % 5) - 2;
            final_result_accumulator += (agent->y * 3);
            return SUCCESS;
        case ACTION_HEAL:
            if (agent->health < MAX_HEALTH) {
                 agent->health += 5;
                 if (agent->health > MAX_HEALTH) agent->health = MAX_HEALTH;
            }
            final_result_accumulator += agent->health;
            return SUCCESS;
    }
    return FAILURE; // Should not be reached
}

void run_computation() {
    final_result_accumulator = 0;
    for (int tick = 0; tick < TICKS_TO_SIMULATE; ++tick) {
        for (int i = 0; i < NUM_AGENTS; ++i) {
            execute_node(&agents[i], agents[i].bt_root_idx);
        }
    }
}

void cleanup() {
    free(agents);
    free(node_pool);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%lld\n", final_result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
