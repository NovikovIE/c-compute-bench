#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA STRUCTURES ---

typedef enum {
    ATOM,       // An atomic proposition, e.g., P(x, y)
    NOT,        // Negation
    AND,        // Conjunction
    OR,         // Disjunction
    FORALL,     // Universal quantifier
    EXISTS      // Existential quantifier
} NodeType;

typedef struct Formula {
    NodeType type;
    int id;         // For ATOM: predicate_id. For Quantifier: variable_id.
    struct Formula* left;  // Unary operand or left binary operand
    struct Formula* right; // Right binary operand
} Formula;

// --- GLOBAL STATE ---

// Parameters
int FORMULA_DEPTH;
int NUM_VARIABLES;
int NUM_QUANTIFIERS;
uint32_t SEED;

// Main data structure
Formula* root_formula;

// Result accumulator
volatile long final_result;

// Memory pool for formula nodes to avoid malloc during computation
Formula* formula_pool;
int formula_pool_idx;
long max_formulas;


// --- BENCHMARK-SPECIFIC IMPLEMENTATION ---

Formula* new_node() {
    if (formula_pool_idx >= max_formulas) {
        fprintf(stderr, "FATAL: Formula pool exhausted. Increase pool size.\n");
        exit(1);
    }
    Formula* node = &formula_pool[formula_pool_idx++];
    node->left = node->right = NULL;
    return node;
}

Formula* generate_formula(int depth, int* quantifiers_left) {
    // Base case: if max depth is reached, create an atomic formula.
    if (depth <= 0) {
        Formula* atom = new_node();
        atom->type = ATOM;
        atom->id = mt_rand() % NUM_VARIABLES; // Atom refers to a variable
        return atom;
    }

    int choice = mt_rand() % 100;
    Formula* node = new_node();

    // Give quantifiers a higher chance of being chosen if available
    if (*quantifiers_left > 0 && choice < 40) {
        (*quantifiers_left)--;
        node->type = (mt_rand() % 2 == 0) ? FORALL : EXISTS;
        node->id = mt_rand() % NUM_VARIABLES; // Quantifier binds a variable
        node->left = generate_formula(depth - 1, quantifiers_left);
    } else if (choice < 80) {
        // Binary connectives: AND, OR
        node->type = (mt_rand() % 2 == 0) ? AND : OR;
        node->left = generate_formula(depth - 1, quantifiers_left);
        node->right = generate_formula(depth - 1, quantifiers_left);
    } else if (choice < 95) {
        // Unary connective: NOT
        node->type = NOT;
        node->left = generate_formula(depth - 1, quantifiers_left);
    } else {
        // Leaf node: ATOM
        node->type = ATOM;
        node->id = mt_rand() % NUM_VARIABLES;
    }
    return node;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <formula_depth> <num_variables> <num_quantifiers> <seed>\n", argv[0]);
        exit(1);
    }

    FORMULA_DEPTH = atoi(argv[1]);
    NUM_VARIABLES = atoi(argv[2]);
    NUM_QUANTIFIERS = atoi(argv[3]);
    SEED = (uint32_t)atoi(argv[4]);

    mt_seed(SEED);

    // Estimate max nodes. A complete binary tree of depth D has 2^(D+1)-1 nodes.
    // We add a large buffer for safety.
    max_formulas = (1L << (FORMULA_DEPTH + 1)) + 1000;
    formula_pool = (Formula*)malloc(max_formulas * sizeof(Formula));
    if (!formula_pool) {
        fprintf(stderr, "FATAL: Memory allocation failed for formula pool.\n");
        exit(1);
    }
    formula_pool_idx = 0;
    final_result = 0;

    int quantifiers_to_generate = NUM_QUANTIFIERS;
    root_formula = generate_formula(FORMULA_DEPTH, &quantifiers_to_generate);
}

// Simulates the proof search process by recursively exploring the formula structure.
// The core computational workload comes from quantifier expansion.
void explore_and_accumulate(Formula* f, int depth) {
    if (f == NULL || depth > FORMULA_DEPTH * 2) {
        // Guard against infinite recursion in case of cyclic structures (not possible here)
        // or excessively deep paths.
        return;
    }

    // Perform some work at each node to prevent dead code elimination
    final_result += f->type + f->id * depth;

    switch (f->type) {
        case ATOM:
            // Leaf node, do nothing further.
            break;
        case NOT:
            explore_and_accumulate(f->left, depth + 1);
            break;
        case AND:
        case OR:
            explore_and_accumulate(f->left, depth + 1);
            explore_and_accumulate(f->right, depth + 1);
            break;
        case FORALL:
        case EXISTS:
            // This simulates quantifier instantiation (Gamma/Delta rules in tableau).
            // The prover would instantiate the sub-formula for different terms.
            // We simulate this by recursively exploring the sub-formula N times,
            // where N is the number of variables, representing the set of available terms.
            // This is the primary source of computational complexity.
            for (int i = 0; i < NUM_VARIABLES; ++i) {
                explore_and_accumulate(f->left, depth + 1);
            }
            break;
    }
}

void run_computation() {
    final_result = 0;
    if (root_formula) {
        explore_and_accumulate(root_formula, 0);
    }
}

void cleanup() {
    free(formula_pool);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout to ensure computation is not optimized away.
    printf("%ld\n", final_result);

    // Print the timing information to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
