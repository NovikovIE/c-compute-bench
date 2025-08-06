#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Benchmark Globals ---
int NUM_VERTICES;
int NUM_EDGES;
int NUM_COLORS;

// Adjacency matrix stored as a 1D array
int *adjacency_matrix;
long long total_colorings_result;

// --- Forward Declarations ---
long long count_colorings_recursive(int vertex_index, int* colors);

// --- Benchmark Implementation ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_vertices num_edges num_colors seed\n", argv[0]);
        exit(1);
    }

    NUM_VERTICES = atoi(argv[1]);
    NUM_EDGES = atoi(argv[2]);
    NUM_COLORS = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    if (NUM_VERTICES <= 0 || NUM_EDGES < 0 || NUM_COLORS <= 0) {
        fprintf(stderr, "ERROR: Parameters must be positive integers.\n");
        exit(1);
    }

    long long max_edges = (long long)NUM_VERTICES * (NUM_VERTICES - 1) / 2;
    if (NUM_EDGES > max_edges) {
        fprintf(stderr, "ERROR: Number of edges exceeds maximum possible for a simple graph.\n");
        exit(1);
    }

    // Allocate and initialize adjacency matrix (using calloc for zero-initialization)
    adjacency_matrix = (int*)calloc((size_t)NUM_VERTICES * NUM_VERTICES, sizeof(int));
    if (!adjacency_matrix) {
        fprintf(stderr, "ERROR: Failed to allocate memory for adjacency matrix.\n");
        exit(1);
    }

    // Generate a random simple graph
    int edges_added = 0;
    int u, v;
    while(edges_added < NUM_EDGES) {
        u = mt_rand() % NUM_VERTICES;
        v = mt_rand() % NUM_VERTICES;

        // Ensure no self-loops and the edge doesn't already exist
        if (u == v || adjacency_matrix[u * NUM_VERTICES + v] == 1) {
            continue;
        }

        // Add edge for undirected graph
        adjacency_matrix[u * NUM_VERTICES + v] = 1;
        adjacency_matrix[v * NUM_VERTICES + u] = 1;
        edges_added++;
    }
}

long long count_colorings_recursive(int v_idx, int* colors) {
    // Base case: If all vertices are colored, we found one valid coloring
    if (v_idx == NUM_VERTICES) {
        return 1;
    }

    long long count = 0;

    // Try to assign a color to the current vertex v_idx
    for (int c = 1; c <= NUM_COLORS; c++) {
        bool is_safe = true;
        // Check if this color conflicts with any adjacent, already-colored vertices
        for (int i = 0; i < v_idx; i++) {
            if (adjacency_matrix[v_idx * NUM_VERTICES + i] == 1 && colors[i] == c) {
                is_safe = false;
                break;
            }
        }

        if (is_safe) {
            colors[v_idx] = c;
            count += count_colorings_recursive(v_idx + 1, colors);
            // Backtracking is implicit: colors[v_idx] will be overwritten by next iteration
        }
    }

    return count;
}

void run_computation() {
    // This array holds the color assigned to each vertex during recursion
    int* colors = (int*)malloc(NUM_VERTICES * sizeof(int));
    if (!colors) {
        fprintf(stderr, "ERROR: Failed to allocate memory for colors array.\n");
        exit(1);
    }

    // Start the recursive counting from the first vertex (index 0)
    total_colorings_result = count_colorings_recursive(0, colors);

    free(colors);
}

void cleanup() {
    free(adjacency_matrix);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", total_colorings_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
