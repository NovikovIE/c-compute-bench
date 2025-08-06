#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) --- DO NOT MODIFY ---
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
typedef struct {
    int num_vertices;
    int** adj_matrix; // Adjacency matrix for the graph
    int* path;         // Current path in the backtracking search
    int* visited;      // Visited vertices for the current path
    long long cycle_count; // Result: number of Hamiltonian cycles found
} BenchmarkData;

BenchmarkData g_data;

// Recursive utility function for finding Hamiltonian cycles
static void find_cycles_recursive(int pos) {
    // Base case: If all vertices are included in the path
    if (pos == g_data.num_vertices) {
        // Check if there is an edge from the last vertex to the first vertex
        if (g_data.adj_matrix[g_data.path[pos - 1]][g_data.path[0]]) {
            g_data.cycle_count++;
        }
        return;
    }

    // Try different vertices as the next candidate in the path
    for (int v = 1; v < g_data.num_vertices; v++) {
        // Check if this vertex can be added to the path
        // 1. Edge exists from previous vertex in path
        // 2. Vertex 'v' has not been visited yet
        if (g_data.adj_matrix[g_data.path[pos - 1]][v] && !g_data.visited[v]) {
            g_data.path[pos] = v;
            g_data.visited[v] = 1;

            // Recur for the next position
            find_cycles_recursive(pos + 1);

            // Backtrack: remove vertex from path and mark as unvisited
            g_data.visited[v] = 0;
        }
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <graph_density> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_vertices = atoi(argv[1]);
    double graph_density = atof(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);
    
    if (g_data.num_vertices <= 0) {
        fprintf(stderr, "FATAL: num_vertices must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for graph and helper arrays
    g_data.adj_matrix = (int**)malloc(g_data.num_vertices * sizeof(int*));
    for (int i = 0; i < g_data.num_vertices; i++) {
        g_data.adj_matrix[i] = (int*)calloc(g_data.num_vertices, sizeof(int));
    }
    g_data.path = (int*)malloc(g_data.num_vertices * sizeof(int));
    g_data.visited = (int*)malloc(g_data.num_vertices * sizeof(int));

    // Generate a random undirected graph based on density
    uint32_t max_rand = 0xFFFFFFFFU;
    for (int i = 0; i < g_data.num_vertices; i++) {
        for (int j = i + 1; j < g_data.num_vertices; j++) {
            if ((double)mt_rand() / max_rand < graph_density) {
                g_data.adj_matrix[i][j] = 1;
                g_data.adj_matrix[j][i] = 1;
            }
        }
    }
    
    g_data.cycle_count = 0;
}

void run_computation() {
    // Initialize helper arrays
    memset(g_data.visited, 0, g_data.num_vertices * sizeof(int));

    // Start with vertex 0 as the first vertex in the path
    g_data.path[0] = 0;
    g_data.visited[0] = 1;

    // Start the recursive search from the second position
    find_cycles_recursive(1);
}

void cleanup() {
    for (int i = 0; i < g_data.num_vertices; i++) {
        free(g_data.adj_matrix[i]);
    }
    free(g_data.adj_matrix);
    free(g_data.path);
    free(g_data.visited);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    cleanup();

    // Print result to stdout
    printf("%lld\n", g_data.cycle_count);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
