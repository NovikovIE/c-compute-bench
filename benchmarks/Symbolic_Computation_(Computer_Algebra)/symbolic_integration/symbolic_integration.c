#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// Benchmark: Symbolic Integration
// Description: Manipulates a large expression tree by applying symbolic integration rules.
// The main computational work involves traversing an existing expression tree and building a
// new, more complex tree representing the integral. This benchmark is characterized by
// irregular memory access patterns (tree traversal) and high memory allocation overhead (building the new tree).

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Benchmark Parameters
static int NUM_VARIABLES;
static int EXPRESSION_TREE_DEPTH;
static const int INTEGRATION_VAR_IDX = 0; // Always integrate with respect to x_0

// Expression node types
typedef enum {
    CONST,   // Constant value
    VAR,     // Variable, e.g., x_i
    ADD,     // Addition operator
    MUL,     // Multiplication operator
    POW      // Power operator
} NodeType;

// Expression tree node
typedef struct Node {
    NodeType type;
    int value; // Constant value or variable index
    int has_integration_var; // Flag: 1 if this subtree depends on the integration variable
    struct Node *left;
    struct Node *right;
} Node;

// Global roots for the original and integrated expression trees
static Node* root_expression = NULL;
static Node* integrated_expression = NULL;

// Final result to prevent dead-code elimination
static long long final_result = 0;

// --- FORWARD DECLARATIONS of helper functions ---
Node* create_node(NodeType type, int value);
Node* generate_expression_recursive(int depth);
Node* copy_tree(const Node* node);
void free_tree(Node* node);
Node* integrate(const Node* node);
long long count_nodes(const Node* node);

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_variables> <expression_tree_depth> <seed>\n", argv[0]);
        exit(1);
    }
    NUM_VARIABLES = atoi(argv[1]);
    EXPRESSION_TREE_DEPTH = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    if (NUM_VARIABLES <= 0 || EXPRESSION_TREE_DEPTH < 0) {
        fprintf(stderr, "FATAL: Invalid parameters.\n");
        exit(1);
    }

    root_expression = generate_expression_recursive(EXPRESSION_TREE_DEPTH);
}

void run_computation() {
    integrated_expression = integrate(root_expression);
    // To produce a result, we count the number of nodes in the new tree.
    // This is a proxy for the complexity of the resulting expression.
    final_result = count_nodes(integrated_expression);
}

void cleanup() {
    free_tree(root_expression);
    free_tree(integrated_expression);
    root_expression = NULL;
    integrated_expression = NULL;
}

// --- HELPER FUNCTION IMPLEMENTATIONS ---

Node* create_node(NodeType type, int value) {
    Node* node = (Node*)malloc(sizeof(Node));
    if (!node) {
        fprintf(stderr, "FATAL: malloc failed\n");
        exit(1);
    }
    node->type = type;
    node->value = value;
    node->has_integration_var = 0;
    node->left = NULL;
    node->right = NULL;
    return node;
}

Node* generate_expression_recursive(int depth) {
    if (depth == 0) { // Leaf node
        if (mt_rand() % 2 == 0) { // Constant
            Node* node = create_node(CONST, mt_rand() % 100 + 1);
            return node;
        } else { // Variable
            int var_idx = mt_rand() % NUM_VARIABLES;
            Node* node = create_node(VAR, var_idx);
            if (var_idx == INTEGRATION_VAR_IDX) {
                node->has_integration_var = 1;
            }
            return node;
        }
    } else { // Operator node
        NodeType op_type = (NodeType)(ADD + (mt_rand() % 3));
        Node* node = create_node(op_type, 0);
        node->left = generate_expression_recursive(depth - 1);
        
        if (op_type == POW) {
            // For POW, ensure exponent is a small constant to keep it simple
            node->right = create_node(CONST, (mt_rand() % 3) + 2); 
        } else {
            node->right = generate_expression_recursive(depth - 1);
        }

        node->has_integration_var = node->left->has_integration_var || node->right->has_integration_var;
        return node;
    }
}

Node* copy_tree(const Node* node) {
    if (!node) return NULL;
    Node* new_node = (Node*)malloc(sizeof(Node));
     if (!new_node) {
        fprintf(stderr, "FATAL: malloc failed\n");
        exit(1);
    }
    *new_node = *node; // Copy type, value, flag
    new_node->left = copy_tree(node->left);
    new_node->right = copy_tree(node->right);
    return new_node;
}

void free_tree(Node* node) {
    if (node == NULL) return;
    free_tree(node->left);
    free_tree(node->right);
    free(node);
}

long long count_nodes(const Node* node) {
    if (node == NULL) return 0;
    return 1 + count_nodes(node->left) + count_nodes(node->right);
}

Node* integrate(const Node* node) {
    if (!node) return NULL;

    Node* new_expr = NULL;

    switch (node->type) {
        case CONST:
            // integral(C) dx = C*x
            new_expr = create_node(MUL, 0);
            new_expr->left = create_node(CONST, node->value);
            new_expr->right = create_node(VAR, INTEGRATION_VAR_IDX);
            break;
        case VAR:
            if (node->value == INTEGRATION_VAR_IDX) {
                // integral(x) dx = x^2 (ignoring 1/2 for simplicity)
                new_expr = create_node(POW, 0);
                new_expr->left = create_node(VAR, INTEGRATION_VAR_IDX);
                new_expr->right = create_node(CONST, 2);
            } else {
                // integral(y) dx = y*x
                new_expr = create_node(MUL, 0);
                new_expr->left = create_node(VAR, node->value);
                new_expr->right = create_node(VAR, INTEGRATION_VAR_IDX);
            }
            break;
        case ADD:
            // integral(f+g) = integral(f) + integral(g)
            new_expr = create_node(ADD, 0);
            new_expr->left = integrate(node->left);
            new_expr->right = integrate(node->right);
            break;
        case MUL:
            // Simplified product rule: check if one side is constant w.r.t x
            if (!node->left->has_integration_var) {
                // f is const: integral(f*g) = f * integral(g)
                new_expr = create_node(MUL, 0);
                new_expr->left = copy_tree(node->left);
                new_expr->right = integrate(node->right);
            } else if (!node->right->has_integration_var) {
                // g is const: integral(f*g) = integral(f) * g
                new_expr = create_node(MUL, 0);
                new_expr->left = integrate(node->left);
                new_expr->right = copy_tree(node->right);
            } else {
                // Both depend on x. Dummy operation to cause expression swell.
                // a complex but incorrect rule for computational load
                Node* int_f = integrate(node->left);
                Node* g = copy_tree(node->right);
                Node* term1 = create_node(MUL, 0);
                term1->left = int_f;
                term1->right = g;
                new_expr = term1;
            }
            break;
        case POW:
            // Simplified power rule: integral(x^n) = x^(n+1)
            if (node->left->type == VAR && node->left->value == INTEGRATION_VAR_IDX && node->right->type == CONST) {
                new_expr = create_node(POW, 0);
                new_expr->left = create_node(VAR, INTEGRATION_VAR_IDX);
                new_expr->right = create_node(CONST, node->right->value + 1);
            } else {
                // Other cases too complex, return original expression
                new_expr = copy_tree(node);
            }
            break;
    }
    return new_expr;
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout AFTER computation but BEFORE cleanup frees the data's container if needed
    // (in this case, final_result is a simple long long, so order is flexible)
    printf("%lld\n", final_result);

    cleanup();

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
