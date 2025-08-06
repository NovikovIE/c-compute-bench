#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// --- BENCHMARK DATA STRUCTURES ---
typedef enum {
    CONSTANT, VARIABLE,
    ADD, SUB, MUL, DIV,
    SIN, COS, POW
} NodeType;

typedef struct ExprNode {
    NodeType type;
    union {
        double value; // For CONSTANT
        int var_index; // For VARIABLE
    } data;
    struct ExprNode *left;
    struct ExprNode *right;
} ExprNode;

// Global pointers for data shared between setup, computation, and cleanup
int g_num_variables;
int g_expression_tree_depth;

ExprNode *g_numerator_root = NULL;
ExprNode *g_denominator_root = NULL;
long long g_total_nodes_processed = 0;

// --- FORWARD DECLARATIONS for helper functions ---
ExprNode* differentiate(ExprNode* node, int var_index);
void free_tree(ExprNode* node);

// --- HELPER FUNCTIONS ---

ExprNode* create_node(NodeType type) {
    ExprNode* node = (ExprNode*)malloc(sizeof(ExprNode));
    if (!node) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    node->type = type;
    node->left = NULL;
    node->right = NULL;
    return node;
}

ExprNode* create_constant(double value) {
    ExprNode* node = create_node(CONSTANT);
    node->data.value = value;
    return node;
}

ExprNode* create_variable(int var_index) {
    ExprNode* node = create_node(VARIABLE);
    node->data.var_index = var_index;
    return node;
}

ExprNode* create_unary(NodeType type, ExprNode* child) {
    ExprNode* node = create_node(type);
    node->left = child;
    return node;
}

ExprNode* create_binary(NodeType type, ExprNode* left, ExprNode* right) {
    ExprNode* node = create_node(type);
    node->left = left;
    node->right = right;
    return node;
}

ExprNode* copy_tree(ExprNode* node) {
    if (!node) return NULL;

    ExprNode* new_node = (ExprNode*)malloc(sizeof(ExprNode));
     if (!new_node) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    *new_node = *node; // Copy type and data union

    new_node->left = copy_tree(node->left);
    new_node->right = copy_tree(node->right);
    return new_node;
}

long long count_nodes(ExprNode* node) {
    if (!node) return 0;
    return 1 + count_nodes(node->left) + count_nodes(node->right);
}

ExprNode* create_random_expr(int depth) {
    if (depth == 0 || (mt_rand() % 100) < 25) { // Base case: leaf node
        if (mt_rand() % 2 == 0) {
            return create_constant((double)(mt_rand() % 100) + 1.0);
        } else {
            return create_variable(mt_rand() % g_num_variables);
        }
    } else { // Recursive step: internal node
        uint32_t op_choice = mt_rand() % 7;
        switch (op_choice) {
            case 0: return create_binary(ADD, create_random_expr(depth - 1), create_random_expr(depth - 1));
            case 1: return create_binary(SUB, create_random_expr(depth - 1), create_random_expr(depth - 1));
            case 2: return create_binary(MUL, create_random_expr(depth - 1), create_random_expr(depth - 1));
            case 3: return create_binary(DIV, create_random_expr(depth - 1), create_random_expr(depth - 1));
            case 4: return create_unary(SIN, create_random_expr(depth - 1));
            case 5: return create_unary(COS, create_random_expr(depth - 1));
            case 6: return create_binary(POW, create_random_expr(depth - 1), create_constant((double)(mt_rand() % 4) + 2.0));
        }
    }
    return NULL; // Should not be reached
}

// Symbolic differentiation function
ExprNode* differentiate(ExprNode* node, int var_index) {
    if (!node) return NULL;

    switch (node->type) {
        case CONSTANT:
            return create_constant(0.0);
        case VARIABLE:
            return create_constant(node->data.var_index == var_index ? 1.0 : 0.0);
        case ADD:
            return create_binary(ADD, differentiate(node->left, var_index), differentiate(node->right, var_index));
        case SUB:
            return create_binary(SUB, differentiate(node->left, var_index), differentiate(node->right, var_index));
        case MUL: { // d(uv) = u'v + uv'
            ExprNode* u = node->left;
            ExprNode* v = node->right;
            return create_binary(ADD,
                create_binary(MUL, differentiate(u, var_index), copy_tree(v)),
                create_binary(MUL, copy_tree(u), differentiate(v, var_index))
            );
        }
        case DIV: { // d(u/v) = (u'v - uv') / v^2
            ExprNode* u = node->left;
            ExprNode* v = node->right;
            ExprNode* u_prime = differentiate(u, var_index);
            ExprNode* v_prime = differentiate(v, var_index);
            ExprNode* u_copy = copy_tree(u);
            ExprNode* v_copy1 = copy_tree(v);
            ExprNode* v_copy2 = copy_tree(v);

            ExprNode* numerator = create_binary(SUB, 
                create_binary(MUL, u_prime, v_copy1),
                create_binary(MUL, u_copy, v_prime)
            );
            ExprNode* denominator = create_binary(POW, v_copy2, create_constant(2.0));
            return create_binary(DIV, numerator, denominator);
        }
        case SIN: { // d(sin(u)) = cos(u) * u'
            return create_binary(MUL, 
                create_unary(COS, copy_tree(node->left)), 
                differentiate(node->left, var_index)
            );
        }
        case COS: { // d(cos(u)) = -sin(u) * u'
            return create_binary(MUL,
                create_binary(MUL, create_constant(-1.0), create_unary(SIN, copy_tree(node->left))),
                differentiate(node->left, var_index)
            );
        }
        case POW: { // d(u^c) = c*u^(c-1) * u' (assuming c is constant)
            ExprNode* base = node->left;
            ExprNode* exponent = node->right;
            if (exponent->type == CONSTANT) {
                double c = exponent->data.value;
                if (c == 0.0) return create_constant(0.0);

                return create_binary(MUL,
                    create_binary(MUL,
                        create_constant(c),
                        create_binary(POW, copy_tree(base), create_constant(c - 1.0))
                    ),
                    differentiate(base, var_index)
                );
            } else { // Power is not a constant, treat derivative as 0 for simplicity.
                return create_constant(0.0);
            }
        }
    }
    return NULL; // Should not be reached
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_variables> <expression_tree_depth> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_variables = atoi(argv[1]);
    g_expression_tree_depth = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    g_numerator_root = create_random_expr(g_expression_tree_depth);
    g_denominator_root = create_random_expr(g_expression_tree_depth);
    g_total_nodes_processed = 0;
}

void run_computation() {
    const int L_HOPITAL_ITERATIONS = 3;

    for (int i = 0; i < L_HOPITAL_ITERATIONS; ++i) {
        int var_to_diff = i % g_num_variables;

        ExprNode* new_numerator = differentiate(g_numerator_root, var_to_diff);
        ExprNode* new_denominator = differentiate(g_denominator_root, var_to_diff);

        free_tree(g_numerator_root);
        free_tree(g_denominator_root);

        g_numerator_root = new_numerator;
        g_denominator_root = new_denominator;

        g_total_nodes_processed += count_nodes(g_numerator_root);
        g_total_nodes_processed += count_nodes(g_denominator_root);
    }
}

void cleanup() {
    free_tree(g_numerator_root);
    free_tree(g_denominator_root);
}

void free_tree(ExprNode* node) {
    if (node == NULL) {
        return;
    }
    free_tree(node->left);
    free_tree(node->right);
    free(node);
}


// --- MAIN ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", g_total_nodes_processed);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
