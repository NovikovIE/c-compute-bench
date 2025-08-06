#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <float.h>

// Mersenne Twister (MT19937) generator - DO NOT MODIFY
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
// End of Mersenne Twister

// --- Benchmark Specific Data Structures and Globals ---

// Problem Parameters
static int NUM_VARIABLES;
static int NUM_CONSTRAINTS;

// Problem Data
static double* c; // Objective function coefficients
static double* A; // Constraint matrix (row-major)
static double* b; // Constraint right-hand side
static int* is_integer; // Flag for integer variables

// Branch and Bound Node Structure
typedef struct {
    double* l_bounds;
    double* u_bounds;
} Node;

// Node processing stack (for depth-first search)
static Node* node_stack;
static int stack_top;
static int stack_capacity;

// Benchmark result
static int nodes_processed;
static double best_objective_value;

// A helper to generate a random double between min and max
double rand_double(double min, double max) {
    return min + (double)mt_rand() / (double)UINT32_MAX * (max - min);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_variables> <num_constraints> <integer_variable_percentage> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_VARIABLES = atoi(argv[1]);
    NUM_CONSTRAINTS = atoi(argv[2]);
    int integer_percentage = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    // Allocate problem data
    c = (double*)malloc(NUM_VARIABLES * sizeof(double));
    A = (double*)malloc(NUM_CONSTRAINTS * NUM_VARIABLES * sizeof(double));
    b = (double*)malloc(NUM_CONSTRAINTS * sizeof(double));
    is_integer = (int*)malloc(NUM_VARIABLES * sizeof(int));

    if (!c || !A || !b || !is_integer) {
        fprintf(stderr, "Failed to allocate memory for problem data.\n");
        exit(1);
    }

    // Generate objective function coefficients
    for (int i = 0; i < NUM_VARIABLES; i++) {
        c[i] = rand_double(-10.0, 10.0);
    }

    // Generate constraint matrix A
    for (int i = 0; i < NUM_CONSTRAINTS * NUM_VARIABLES; i++) {
        A[i] = rand_double(-5.0, 5.0);
    }

    // Generate constraint bounds b (ensure they are non-negative)
    for (int i = 0; i < NUM_CONSTRAINTS; i++) {
        b[i] = rand_double(10.0, 100.0);
    }

    // Determine which variables are integer
    for (int i = 0; i < NUM_VARIABLES; i++) {
        is_integer[i] = (rand_double(0.0, 100.0) < integer_percentage) ? 1 : 0;
    }

    // Initialize Branch and Bound state
    best_objective_value = -DBL_MAX;
    nodes_processed = 0;
    
    // Setup node stack
    stack_capacity = 2 * NUM_VARIABLES * 200; // Heuristic capacity
    node_stack = (Node*)malloc(stack_capacity * sizeof(Node));
    if (!node_stack) {
        fprintf(stderr, "Failed to allocate memory for node stack.\n");
        exit(1);
    }
    stack_top = -1;

    // Create and push the root node
    Node root_node;
    root_node.l_bounds = (double*)malloc(NUM_VARIABLES * sizeof(double));
    root_node.u_bounds = (double*)malloc(NUM_VARIABLES * sizeof(double));
    if (!root_node.l_bounds || !root_node.u_bounds) {
        fprintf(stderr, "Failed to allocate memory for root node bounds.\n");
        exit(1);
    }

    for (int i = 0; i < NUM_VARIABLES; i++) {
        root_node.l_bounds[i] = 0.0;
        root_node.u_bounds[i] = 20.0; // A generic upper bound
    }
    
    stack_top++;
    node_stack[stack_top] = root_node;
}

// Simulates solving the LP relaxation. This is the main workload of the benchmark.
void solve_lp_relaxation(const Node* node, double* solution, double* objective) {
    // This is a proxy for a real LP solver (like Simplex).
    // It performs matrix-vector multiplications to generate a realistic workload.
    const int SOLVER_ITERATIONS = 15;
    double* temp_vec = (double*)malloc(NUM_CONSTRAINTS * sizeof(double));
    if (!temp_vec) return; // Fail silently in computation

    // Initialize solution with midpoints of bounds
    for (int i = 0; i < NUM_VARIABLES; i++) {
        solution[i] = (node->l_bounds[i] + node->u_bounds[i]) / 2.0;
    }

    // Simulate iterative solver work
    for (int iter = 0; iter < SOLVER_ITERATIONS; iter++) {
        // temp_vec = A * solution
        for (int i = 0; i < NUM_CONSTRAINTS; i++) {
            temp_vec[i] = 0.0;
            for (int j = 0; j < NUM_VARIABLES; j++) {
                temp_vec[i] += A[i * NUM_VARIABLES + j] * solution[j];
            }
        }
        // solution = A^T * temp_vec (not exactly, but generates work)
        for (int j = 0; j < NUM_VARIABLES; j++) {
            double val = 0.0;
            for (int i = 0; i < NUM_CONSTRAINTS; i++) {
                val += A[i * NUM_VARIABLES + j] * temp_vec[i];
            }
            // Normalize and update to keep values from exploding
            solution[j] = fmod(fabs(val), node->u_bounds[j]);
        }
    }
    
    free(temp_vec);

    // Calculate objective value
    *objective = 0.0;
    for (int i = 0; i < NUM_VARIABLES; i++) {
        *objective += c[i] * solution[i];
    }
    
    // To ensure branching, find the first integer variable and make it fractional
    for (int i = 0; i < NUM_VARIABLES; i++) {
        if (is_integer[i]) {
            double val = solution[i];
            if (fabs(val - floor(val)) > 1e-6) {
                // already fractional, good.
            } else {
                // make it fractional
                solution[i] = floor(val) + 0.5;
                if (solution[i] > node->u_bounds[i]) {
                    solution[i] = node->l_bounds[i];
                }
            }
            break;
        }
    }
}

void run_computation() {
    double* current_solution = (double*)malloc(NUM_VARIABLES * sizeof(double));
    if (!current_solution) return;
    
    // Limit total nodes to ensure termination for benchmark purposes
    const int MAX_NODES = 15000; 

    while (stack_top > -1 && nodes_processed < MAX_NODES) {
        // Pop a node from the stack (DFS)
        Node current_node = node_stack[stack_top--];
        nodes_processed++;
        
        double current_objective;
        solve_lp_relaxation(&current_node, current_solution, &current_objective);

        // Pruning by objective: if the relaxed solution is worse than our best
        // integer solution, no integer solution in this branch can be better.
        if (current_objective <= best_objective_value) {
            free(current_node.l_bounds);
            free(current_node.u_bounds);
            continue;
        }

        // Check for integer feasibility
        int first_fractional_var = -1;
        int is_feasible = 1;
        for (int i = 0; i < NUM_VARIABLES; i++) {
            if (is_integer[i]) {
                double val = current_solution[i];
                if (fabs(val - round(val)) > 1e-6) {
                    is_feasible = 0;
                    first_fractional_var = i;
                    break; 
                }
            }
        }
        
        if (is_feasible) {
            // Found a new, better integer solution
            if (current_objective > best_objective_value) {
                best_objective_value = current_objective;
            }
        } else {
            // Branch on the first fractional integer variable
            if (stack_top + 2 >= stack_capacity) {
                // Out of stack space, prune this branch
                free(current_node.l_bounds);
                free(current_node.u_bounds);
                continue;
            }
            int branch_var = first_fractional_var;
            
            // Create two new child nodes
            // Child 1: x_i <= floor(x_i)
            Node child1;
            child1.l_bounds = (double*)malloc(NUM_VARIABLES * sizeof(double));
            child1.u_bounds = (double*)malloc(NUM_VARIABLES * sizeof(double));
            memcpy(child1.l_bounds, current_node.l_bounds, NUM_VARIABLES * sizeof(double));
            memcpy(child1.u_bounds, current_node.u_bounds, NUM_VARIABLES * sizeof(double));
            child1.u_bounds[branch_var] = floor(current_solution[branch_var]);
            
            // Child 2: x_i >= ceil(x_i)
            Node child2;
            child2.l_bounds = (double*)malloc(NUM_VARIABLES * sizeof(double));
            child2.u_bounds = (double*)malloc(NUM_VARIABLES * sizeof(double));
            memcpy(child2.l_bounds, current_node.l_bounds, NUM_VARIABLES * sizeof(double));
            memcpy(child2.u_bounds, current_node.u_bounds, NUM_VARIABLES * sizeof(double));
            child2.l_bounds[branch_var] = ceil(current_solution[branch_var]);

            // Push children to stack if valid
            if (child1.l_bounds[branch_var] <= child1.u_bounds[branch_var]) {
                stack_top++;
                node_stack[stack_top] = child1;
            } else {
                free(child1.l_bounds);
                free(child1.u_bounds);
            }
            if (child2.l_bounds[branch_var] <= child2.u_bounds[branch_var]) {
                stack_top++;
                node_stack[stack_top] = child2;
            } else {
                free(child2.l_bounds);
                free(child2.u_bounds);
            }
        }
        
        free(current_node.l_bounds);
        free(current_node.u_bounds);
    }
    
    // Free any remaining nodes on the stack
    while (stack_top > -1) {
        Node remaining_node = node_stack[stack_top--];
        free(remaining_node.l_bounds);
        free(remaining_node.u_bounds);
    }

    free(current_solution);
}

void cleanup() {
    free(c);
    free(A);
    free(b);
    free(is_integer);
    free(node_stack);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Use nodes_processed as the final result to prevent dead code elimination
    // and provide a deterministic integer output for verification.
    printf("%d\n", nodes_processed); 

    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
