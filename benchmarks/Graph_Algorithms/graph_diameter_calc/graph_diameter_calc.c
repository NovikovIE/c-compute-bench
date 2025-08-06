/**
 * @file graph_diameter_calc.c
 * @brief A benchmark for calculating the diameter of a graph.
 * 
 * This program generates a random undirected graph and then computes its diameter,
 * which is the longest shortest path between any pair of vertices. The computation
 * is performed using the Floyd-Warshall algorithm, which has a time complexity
 * of O(V^3), where V is the number of vertices. This makes the benchmark's
 * runtime highly dependent on the number of vertices, allowing for precise 
 * workload scaling.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

// --- Mersenne Twister (MT19937) PRNG --- (DO NOT MODIFY) 
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
// --- End of Mersenne Twister --- 

// --- Benchmark Globals ---
int num_vertices;
int** dist_matrix;     // Adjacency matrix for Floyd-Warshall
int final_result;      // Stores the calculated graph diameter

// Represents infinity; chosen to prevent overflow during addition
#define INF (INT_MAX / 2)

/**
 * @brief Parses arguments, allocates memory, and generates the graph data.
 * 
 * Creates an adjacency matrix representing an undirected graph. Edges are added
 * randomly. Distances are initialized to 1 for direct edges, 0 for self-loops,
 * and infinity otherwise.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    int num_edges = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_vertices <= 0 || num_edges < 0) {
        fprintf(stderr, "Error: number of vertices and edges must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate the adjacency matrix
    dist_matrix = (int**)malloc(num_vertices * sizeof(int*));
    if (!dist_matrix) {
        perror("malloc failed for rows");
        exit(1);
    }
    for (int i = 0; i < num_vertices; ++i) {
        dist_matrix[i] = (int*)malloc(num_vertices * sizeof(int));
        if (!dist_matrix[i]) {
            perror("malloc failed for columns");
            // Naive cleanup, but sufficient for benchmark failure
            exit(1);
        }
    }

    // Initialize the matrix
    for (int i = 0; i < num_vertices; ++i) {
        for (int j = 0; j < num_vertices; ++j) {
            if (i == j) {
                dist_matrix[i][j] = 0;
            } else {
                dist_matrix[i][j] = INF;
            }
        }
    }

    // Add random edges (undirected graph)
    for (int i = 0; i < num_edges; ++i) {
        int u = mt_rand() % num_vertices;
        int v = mt_rand() % num_vertices;
        if (u != v) {
            dist_matrix[u][v] = 1;
            dist_matrix[v][u] = 1;
        }
    }
}

/**
 * @brief Runs the core computation: Floyd-Warshall All-Pairs Shortest Path.
 * 
 * After finding all shortest paths, it finds the maximum shortest path length,
 * which is the graph's diameter. The result is stored in a global variable to
 * prevent dead-code elimination.
 */
void run_computation() {
    // Floyd-Warshall algorithm
    for (int k = 0; k < num_vertices; ++k) {
        for (int i = 0; i < num_vertices; ++i) {
            for (int j = 0; j < num_vertices; ++j) {
                if (dist_matrix[i][k] != INF && dist_matrix[k][j] != INF) {
                     if (dist_matrix[i][k] + dist_matrix[k][j] < dist_matrix[i][j]) {
                        dist_matrix[i][j] = dist_matrix[i][k] + dist_matrix[k][j];
                    }
                }
            }
        }
    }

    // Find the diameter (maximum shortest path)
    int diameter = 0;
    for (int i = 0; i < num_vertices; ++i) {
        for (int j = 0; j < num_vertices; ++j) {
            // The diameter is the max of all finite path lengths
            if (dist_matrix[i][j] != INF && dist_matrix[i][j] > diameter) {
                diameter = dist_matrix[i][j];
            }
        }
    }

    final_result = diameter;
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    if (dist_matrix) {
        for (int i = 0; i < num_vertices; ++i) {
            free(dist_matrix[i]);
        }
        free(dist_matrix);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
