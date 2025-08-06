#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

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

// --- BENCHMARK DATA STRUCTURES ---
typedef enum {
    OP_CONST,
    OP_VAR,
    OP_ADD,
    OP_MUL,
    OP_SIN,
    OP_COS,
    OP_EXP
} OpType;

typedef struct ExprNode {
    OpType type;
    double value;          
    int var_index;         
    struct ExprNode *left; 
    struct ExprNode *right;
} ExprNode;

// --- GLOBAL STATE ---
typedef struct {
    int num_variables;
    int expression_tree_depth;
    int series_order;

    ExprNode* base_expression;
    ExprNode** final_derivatives;
    size_t num_final_derivatives;

    long long total_nodes_result;
} BenchmarkState;

static BenchmarkState G;

// --- FORWARD DECLARATIONS ---
ExprNode* differentiate(const ExprNode* expr, int var_idx);
void free_tree(ExprNode* node);
ExprNode* copy_tree(const ExprNode* node);

// --- HELPER FUNCTIONS ---
ExprNode* create_node(OpType type, double value, int var_index, ExprNode* left, ExprNode* right) {
    ExprNode* node = (ExprNode*)malloc(sizeof(ExprNode));
    if (!node) {
        perror("Failed to allocate expression node");
        exit(1);
    }
    node->type = type;
    node->value = value;
    node->var_index = var_index;
    node->left = left;
    node->right = right;
    return node;
}

ExprNode* create_random_expr(int depth) {
    if (depth <= 0) { // Leaf node
        if (mt_rand() % 2 == 0) {
            return create_node(OP_CONST, (double)(mt_rand() % 1000) / 100.0, 0, NULL, NULL);
        } else {
            return create_node(OP_VAR, 0.0, mt_rand() % G.num_variables, NULL, NULL);
        }
    }

    uint32_t op_choice = mt_rand() % 4;
    switch (op_choice) {
        case 0: // ADD
            return create_node(OP_ADD, 0.0, 0, create_random_expr(depth - 1), create_random_expr(depth - 1));
        case 1: // MUL
            return create_node(OP_MUL, 0.0, 0, create_random_expr(depth - 1), create_random_expr(depth - 1));
        case 2: // SIN
            return create_node(OP_SIN, 0.0, 0, create_random_expr(depth - 1), NULL);
        case 3: // EXP
        default:
            return create_node(OP_EXP, 0.0, 0, create_random_expr(depth - 1), NULL);
    }
}

long long count_nodes(const ExprNode* node) {
    if (!node) return 0;
    return 1 + count_nodes(node->left) + count_nodes(node->right);
}

// --- SYMBOLIC DIFFERENTIATION CORE ---
ExprNode* differentiate(const ExprNode* expr, int var_idx) {
    if (!expr) return NULL;

    switch (expr->type) {
        case OP_CONST:
            return create_node(OP_CONST, 0.0, 0, NULL, NULL);
        case OP_VAR:
            return create_node(OP_CONST, (expr->var_index == var_idx) ? 1.0 : 0.0, 0, NULL, NULL);
        case OP_ADD:
            return create_node(OP_ADD, 0.0, 0, 
                               differentiate(expr->left, var_idx), 
                               differentiate(expr->right, var_idx));
        case OP_MUL: // Product rule: (uv)' = u'v + uv'
            return create_node(OP_ADD, 0.0, 0,
                create_node(OP_MUL, 0.0, 0, differentiate(expr->left, var_idx), copy_tree(expr->right)),
                create_node(OP_MUL, 0.0, 0, copy_tree(expr->left), differentiate(expr->right, var_idx)));
        case OP_SIN: // Chain rule: sin(u)' = cos(u) * u'
            return create_node(OP_MUL, 0.0, 0, 
                               create_node(OP_COS, 0.0, 0, copy_tree(expr->left), NULL), 
                               differentiate(expr->left, var_idx));
        case OP_COS: // Chain rule: cos(u)' = -sin(u) * u'
            return create_node(OP_MUL, 0.0, 0,
                create_node(OP_MUL, 0.0, 0, create_node(OP_CONST, -1.0, 0, NULL, NULL), create_node(OP_SIN, 0.0, 0, copy_tree(expr->left), NULL)),
                differentiate(expr->left, var_idx));
        case OP_EXP: // Chain rule: exp(u)' = exp(u) * u'
            return create_node(OP_MUL, 0.0, 0, 
                               copy_tree(expr), 
                               differentiate(expr->left, var_idx));
        default:
            return NULL;
    }
}

// --- MEMORY MANAGEMENT WRAPPERS ---
void free_tree(ExprNode* node) {
    if (!node) return;
    free_tree(node->left);
    free_tree(node->right);
    free(node);
}

ExprNode* copy_tree(const ExprNode* node) {
    if (!node) return NULL;
    return create_node(node->type, node->value, node->var_index, 
                       copy_tree(node->left), 
                       copy_tree(node->right));
}

// --- BENCHMARK FUNCTIONS --- 
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_variables expression_tree_depth series_order seed\n", argv[0]);
        exit(1);
    }

    G.num_variables = atoi(argv[1]);
    G.expression_tree_depth = atoi(argv[2]);
    G.series_order = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    G.base_expression = create_random_expr(G.expression_tree_depth);
    G.final_derivatives = NULL;
    G.num_final_derivatives = 0;
    G.total_nodes_result = 0;
}

void run_computation() {
    size_t num_current = 1;
    ExprNode** current_trees = (ExprNode**)malloc(sizeof(ExprNode*));
    if (!current_trees) { perror("malloc"); exit(1); }
    current_trees[0] = G.base_expression;

    for (int i = 0; i < G.series_order; ++i) {
        size_t num_next = num_current * G.num_variables;
        ExprNode** next_trees = (ExprNode**)malloc(num_next * sizeof(ExprNode*));
        if (!next_trees) { perror("malloc"); exit(1); }

        for (size_t j = 0; j < num_current; ++j) {
            for (int k = 0; k < G.num_variables; ++k) {
                next_trees[j * G.num_variables + k] = differentiate(current_trees[j], k);
            }
        }

        if (i > 0) {
            for (size_t j = 0; j < num_current; ++j) {
                free_tree(current_trees[j]);
            }
        }
        free(current_trees);

        current_trees = next_trees;
        num_current = num_next;
    }

    G.final_derivatives = current_trees;
    G.num_final_derivatives = num_current;

    G.total_nodes_result = 0;
    for (size_t i = 0; i < G.num_final_derivatives; ++i) {
        G.total_nodes_result += count_nodes(G.final_derivatives[i]);
    }
}

void cleanup() {
    free_tree(G.base_expression);
    if (G.final_derivatives) {
        for (size_t i = 0; i < G.num_final_derivatives; ++i) {
            free_tree(G.final_derivatives[i]);
        }
        free(G.final_derivatives);
    }
    
    // paranoia
    G.base_expression = NULL;
    G.final_derivatives = NULL;
    G.num_final_derivatives = 0;
}

// --- MAIN --- 
int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    printf("%lld\n", G.total_nodes_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
