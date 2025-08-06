/**
 * @file zigzag_persistence_calculation.c
 * @brief A benchmark simulating the calculation of zigzag persistence in computational topology.
 * 
 * @details This program simulates a core component of Topological Data Analysis (TDA), specifically
 * the computation of zigzag persistence. In TDA, we study the shape of data by constructing a sequence
 * of simplicial complexes (a filtration) and tracking topological features (like connected components, holes, voids)
 * as they appear and disappear.
 * 
 * Standard persistence involves a linear sequence of inclusions: K_0 ⊆ K_1 ⊆ ... ⊆ K_n.
 * Zigzag persistence generalizes this to allow for both inclusions and exclusions (additions and removals
 * of simplices), forming a sequence like K_0 ↔ K_1 ↔ ... ↔ K_n. This is useful for analyzing
 * dynamic or time-varying data.
 * 
 * The actual computation involves complex linear algebra (matrix reductions) over a finite field to update
 * homology groups at each step. This benchmark abstracts and simulates this computationally intensive process.
 * 
 * **Simulation Details:**
 * - **setup_benchmark:** Creates a large pool of random simplices (points, edges, triangles, etc.) up to
 *   `max_dimension`. It then generates a `filtration_sequence` of a given `filtration_length`, where each step
 *   is a random ADD or REMOVE operation on a simplex from the pool.
 * - **run_computation:** Iterates through the filtration sequence. For each operation, it performs a series
 *   of mock computations that simulate updating a boundary matrix. The computational workload for processing
 *   a simplex is proportional to its dimension (i.e., the number of its faces).
 * - **Data Structures:**
 *   - `Simplex`: Represents a geometric building block (e.g., a vertex, edge, triangle).
 *   - `FiltrationOp`: An operation in the sequence (ADD/REMOVE a simplex).
 *   - `boundary_state`: A large array representing the state of the boundary matrix being updated.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator --- 
// (Provided verbatim as per requirements)
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

// --- Benchmark Data Structures and Globals ---

// Parameters
static int filtration_length;
static int max_dimension;

// A simplex is a generalized triangle. 0-simplex=vertex, 1-simplex=edge, etc.
typedef struct {
    int dimension;
    int num_vertices;
    int* vertices; // Array of vertex indices
} Simplex;

// An operation in the zigzag filtration
typedef enum { ADD, REMOVE } OpType;
typedef struct {
    OpType type;
    int simplex_id; // ID of the simplex being operated on
} FiltrationOp;

// Global data pointers
static Simplex* simplex_pool = NULL;
static int num_simplices_in_pool;

static FiltrationOp* filtration = NULL;

// A large array simulating the boundary matrix state for computation
static unsigned int* boundary_state = NULL;
static int boundary_state_size;

// Final result accumulator
static long long final_result = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <filtration_length> <max_dimension> <seed>\n", argv[0]);
        exit(1);
    }

    filtration_length = atoi(argv[1]);
    max_dimension = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Derive pool size from parameters. Make it large enough for variety but not excessive.
    num_simplices_in_pool = filtration_length / 4 > 1000 ? filtration_length / 4 : 1000;
    if (num_simplices_in_pool > 500000) num_simplices_in_pool = 500000; // Cap memory usage
    int num_total_vertices = 50 + max_dimension * 5;

    // 1. Generate a pool of random simplices
    simplex_pool = (Simplex*)malloc(num_simplices_in_pool * sizeof(Simplex));
    if (!simplex_pool) { perror("malloc failed for simplex_pool"); exit(1); }

    int* all_vertex_indices = (int*)malloc(num_total_vertices * sizeof(int));
    if (!all_vertex_indices) { perror("malloc failed for vertex_indices"); exit(1); }
    for (int i = 0; i < num_total_vertices; ++i) {
        all_vertex_indices[i] = i;
    }

    for (int i = 0; i < num_simplices_in_pool; ++i) {
        Simplex* s = &simplex_pool[i];
        s->dimension = mt_rand() % (max_dimension + 1);
        s->num_vertices = s->dimension + 1;
        s->vertices = (int*)malloc(s->num_vertices * sizeof(int));
        if (!s->vertices) { perror("malloc failed for simplex vertices"); exit(1); }

        // Select unique random vertices by shuffling the index array and taking the first N
        for (int k = num_total_vertices - 1; k > 0; --k) {
            int j = mt_rand() % (k + 1);
            int temp = all_vertex_indices[k];
            all_vertex_indices[k] = all_vertex_indices[j];
            all_vertex_indices[j] = temp;
        }
        for (int k = 0; k < s->num_vertices; ++k) {
            s->vertices[k] = all_vertex_indices[k];
        }
    }
    free(all_vertex_indices);

    // 2. Generate the filtration sequence of ADD/REMOVE operations
    filtration = (FiltrationOp*)malloc(filtration_length * sizeof(FiltrationOp));
    if (!filtration) { perror("malloc failed for filtration"); exit(1); }

    for (int i = 0; i < filtration_length; ++i) {
        filtration[i].simplex_id = mt_rand() % num_simplices_in_pool;
        filtration[i].type = (mt_rand() % 2 == 0) ? ADD : REMOVE;
    }

    // 3. Initialize the computational state array
    boundary_state_size = num_simplices_in_pool * 16;
    boundary_state = (unsigned int*)malloc(boundary_state_size * sizeof(unsigned int));
    if (!boundary_state) { perror("malloc failed for boundary_state"); exit(1); }

    for (int i = 0; i < boundary_state_size; ++i) {
        boundary_state[i] = mt_rand();
    }
}

void run_computation() {
    final_result = 0;
    for (int i = 0; i < filtration_length; ++i) {
        FiltrationOp op = filtration[i];
        Simplex simplex = simplex_pool[op.simplex_id];

        // This loop simulates the core computation: updating a column of the boundary matrix.
        // The work is proportional to the number of faces of the simplex, which is (dimension + 1).
        for (int j = 0; j < simplex.dimension + 1; ++j) {
            unsigned int hash = (op.simplex_id * 0x9E3779B9) ^ (simplex.vertices[j] * 0xDA442D27);
            
            // Heavier artificial workload to simulate matrix reduction steps.
            for(int k=0; k < 10; ++k) {
                int idx1 = (hash + j + k) % boundary_state_size;
                int idx2 = (hash * 31 + simplex.dimension + k) % boundary_state_size;

                if (op.type == ADD) {
                    boundary_state[idx1] = boundary_state[idx1] * 1664525U + boundary_state[idx2] + 1013904223U;
                } else { // REMOVE
                    boundary_state[idx1] = (boundary_state[idx1] * 22695477U) ^ (boundary_state[idx2] >> 3);
                }
                final_result += (boundary_state[idx1] >> 16) & 0xFF;
            }
        }
    }
}

void cleanup() {
    if (simplex_pool) {
        for (int i = 0; i < num_simplices_in_pool; ++i) {
            free(simplex_pool[i].vertices);
        }
        free(simplex_pool);
    }
    free(filtration);
    free(boundary_state);
}

// --- Main Execution --- 

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout to prevent dead code elimination
    printf("%lld\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
