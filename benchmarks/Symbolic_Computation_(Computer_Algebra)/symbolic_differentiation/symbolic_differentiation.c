#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (Verbatim as provided) ---
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
// --- End Mersenne Twister ---


// --- Benchmark Specific Code ---

// Expression Node Type
typedef enum {
    CONSTANT,
    VARIABLE,
    ADD,
    SUB,
    MUL,
    SIN,
    COS
} NodeType;

// Expression Tree Node
typedef struct ExprNode {
    NodeType type;
    union {
        double value; // For CONSTANT
        int var_index; // For VARIABLE
        struct {      // For binary operators
            struct ExprNode *left;
            struct ExprNode *right;
        } binary;
        struct {      // For unary operators
            struct ExprNode *child;
        } unary;
    } data;
} ExprNode;

// --- Global Data ---
static int g_num_variables;
static int g_expression_tree_depth;
static int g_differentiation_order;
static ExprNode* g_root_expression = NULL; // Used for setup and cleanup
static long long g_final_node_count = 0;

// --- Forward Declarations ---
void free_tree(ExprNode* node);
ExprNode* copy_tree(ExprNode* node);
ExprNode* differentiate(ExprNode* expr, int var_idx);

// --- Helper functions for creating nodes ---
ExprNode* create_node(NodeType type) {
    ExprNode* node = (ExprNode*)malloc(sizeof(ExprNode));
    if (!node) {
        fprintf(stderr, "FATAL: malloc failed\n");
        exit(1);
    }
    node->type = type;
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

ExprNode* create_binary_op(NodeType type, ExprNode* left, ExprNode* right) {
    ExprNode* node = create_node(type);
    node->data.binary.left = left;
    node->data.binary.right = right;
    return node;
}

ExprNode* create_unary_op(NodeType type, ExprNode* child) {
    ExprNode* node = create_node(type);
    node->data.unary.child = child;
    return node;
}

// --- Tree Manipulation Functions ---

// Recursively generate a random expression tree
ExprNode* generate_tree(int depth) {
    if (depth <= 0) {
        if (mt_rand() % 2 == 0) {
            return create_constant((double)(mt_rand() % 100) / 10.0 + 1.0);
        } else {
            return create_variable(mt_rand() % g_num_variables);
        }
    }

    int op_type_rand = mt_rand() % 100;
	if (op_type_rand < 25) { // ADD
		return create_binary_op(ADD, generate_tree(depth - 1), generate_tree(depth - 1));
	} else if (op_type_rand < 50) { // SUB
		return create_binary_op(SUB, generate_tree(depth - 1), generate_tree(depth - 1));
	} else if (op_type_rand < 80) { // MUL (Higher chance for expression swell)
		return create_binary_op(MUL, generate_tree(depth - 1), generate_tree(depth - 1));
	} else if (op_type_rand < 90) { // SIN
		return create_unary_op(SIN, generate_tree(depth - 1));
	} else { // COS
		return create_unary_op(COS, generate_tree(depth - 1));
	}
}

// Recursively copy an entire tree
ExprNode* copy_tree(ExprNode* node) {
    if (!node) return NULL;
    switch (node->type) {
        case CONSTANT:
            return create_constant(node->data.value);
        case VARIABLE:
            return create_variable(node->data.var_index);
        case ADD: case SUB: case MUL:
            return create_binary_op(node->type, copy_tree(node->data.binary.left), copy_tree(node->data.binary.right));
        case SIN: case COS:
            return create_unary_op(node->type, copy_tree(node->data.unary.child));
        default:
            return NULL;
    }
}

// Recursively free an entire tree
void free_tree(ExprNode* node) {
    if (!node) return;
    switch (node->type) {
        case CONSTANT: case VARIABLE: break;
        case ADD: case SUB: case MUL:
            free_tree(node->data.binary.left);
            free_tree(node->data.binary.right);
            break;
        case SIN: case COS:
            free_tree(node->data.unary.child);
            break;
    }
    free(node);
}

// Recursively count nodes in a tree
long long count_nodes(ExprNode* node) {
    if (!node) return 0;
    long long count = 1;
    switch (node->type) {
         case CONSTANT: case VARIABLE: break;
         case ADD: case SUB: case MUL:
            count += count_nodes(node->data.binary.left);
            count += count_nodes(node->data.binary.right);
            break;
         case SIN: case COS:
            count += count_nodes(node->data.unary.child);
            break;
    }
    return count;
}

// --- Core Symbolic Differentiation ---
ExprNode* differentiate(ExprNode* expr, int var_idx) {
    if (!expr) return NULL;

    switch (expr->type) {
        case CONSTANT: // d(c)/dx = 0
            return create_constant(0.0);
        
        case VARIABLE: // d(x)/dx = 1, d(y)/dx = 0
            return create_constant(expr->data.var_index == var_idx ? 1.0 : 0.0);
        
        case ADD: // d(u+v)/dx = du/dx + dv/dx
            return create_binary_op(ADD, 
                differentiate(expr->data.binary.left, var_idx),
                differentiate(expr->data.binary.right, var_idx));
        
        case SUB: // d(u-v)/dx = du/dx - dv/dx
            return create_binary_op(SUB, 
                differentiate(expr->data.binary.left, var_idx),
                differentiate(expr->data.binary.right, var_idx));

        case MUL: { // d(u*v)/dx = u*(dv/dx) + v*(du/dx)
            ExprNode* u = expr->data.binary.left;
            ExprNode* v = expr->data.binary.right;
            return create_binary_op(ADD,
                create_binary_op(MUL, copy_tree(u), differentiate(v, var_idx)),
                create_binary_op(MUL, copy_tree(v), differentiate(u, var_idx)));
        }

        case SIN: { // d(sin(u))/dx = cos(u) * du/dx
            ExprNode* u = expr->data.unary.child;
            return create_binary_op(MUL,
                create_unary_op(COS, copy_tree(u)),
                differentiate(u, var_idx));
        }

        case COS: { // d(cos(u))/dx = -sin(u) * du/dx
            ExprNode* u = expr->data.unary.child;
            ExprNode* neg_one = create_constant(-1.0);
            ExprNode* sin_u = create_unary_op(SIN, copy_tree(u));
            ExprNode* neg_sin_u = create_binary_op(MUL, neg_one, sin_u);
            return create_binary_op(MUL, neg_sin_u, differentiate(u, var_idx));
        }
    }
    return NULL; // Should not happen
}

// --- Main Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_variables> <expression_tree_depth> <differentiation_order> <seed>\n", argv[0]);
        exit(1);
    }
    
    g_num_variables = atoi(argv[1]);
    g_expression_tree_depth = atoi(argv[2]);
    g_differentiation_order = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);
    
    mt_seed(seed);
    
    g_root_expression = generate_tree(g_expression_tree_depth);
}

void run_computation() {
    ExprNode* current_expr = g_root_expression;
    g_root_expression = NULL; 

    for (int i = 0; i < g_differentiation_order; ++i) {
        int var_to_diff = i % g_num_variables;
        ExprNode* next_expr = differentiate(current_expr, var_to_diff);
        free_tree(current_expr);
        current_expr = next_expr;
    }
    
    g_final_node_count = count_nodes(current_expr);
    g_root_expression = current_expr;
}

void cleanup() {
    free_tree(g_root_expression);
    g_root_expression = NULL;
}

// --- Main driver ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_final_node_count);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
