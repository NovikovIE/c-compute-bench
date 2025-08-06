#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

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

// --- BENCHMARK DATA AND PARAMETERS ---
int V;          // Number of vertices
int E;          // Number of edges

int** adj;      // Adjacency list for the graph
int* degree;    // Degree of each vertex
int* colors;    // Array to store the color of each vertex

int final_result; // The number of colors used

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    V = atoi(argv[1]);
    E = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Two-pass graph generation to create a clean adjacency list
    // Pass 1: Generate a random edge list and calculate degrees
    int* edge_u = (int*)malloc(E * sizeof(int));
    int* edge_v = (int*)malloc(E * sizeof(int));
    if (!edge_u || !edge_v) {
        fprintf(stderr, "Failed to allocate edge list memory\n");
        exit(1);
    }

    degree = (int*)calloc(V, sizeof(int));
    if (!degree) {
        fprintf(stderr, "Failed to allocate degree memory\n");
        exit(1);
    }

    for (int i = 0; i < E; i++) {
        int u = mt_rand() % V;
        int v = mt_rand() % V;
        if (u != v) {
            edge_u[i] = u;
            edge_v[i] = v;
            degree[u]++;
            degree[v]++;
        } else {
          // Store as a sentinel value to be skipped later
          edge_u[i] = -1;
        }
    }

    // Pass 2: Allocate and populate adjacency lists based on final degrees
    adj = (int**)malloc(V * sizeof(int*));
    if (!adj) {
        fprintf(stderr, "Failed to allocate adjacency list memory\n");
        exit(1);
    }
    for (int i = 0; i < V; i++) {
        if (degree[i] > 0) {
            adj[i] = (int*)malloc(degree[i] * sizeof(int));
            if (!adj[i]) {
                 fprintf(stderr, "Failed to allocate neighbor list memory\n");
                 exit(1);
            }
        } else {
            adj[i] = NULL; // No neighbors
        }
    }

    int* current_degree = (int*)calloc(V, sizeof(int));
    if (!current_degree) {
        fprintf(stderr, "Failed to allocate temporary degree counter\n");
        exit(1);
    }

    for (int i = 0; i < E; i++) {
        if (edge_u[i] != -1) {
            int u = edge_u[i];
            int v = edge_v[i];
            adj[u][current_degree[u]++] = v;
            adj[v][current_degree[v]++] = u;
        }
    }

    free(edge_u);
    free(edge_v);
    free(current_degree);

    colors = (int*)malloc(V * sizeof(int));
    if (!colors) {
        fprintf(stderr, "Failed to allocate colors array\n");
        exit(1);
    }
}

void run_computation() {
    for (int i = 0; i < V; i++) {
        colors[i] = -1; // -1 represents "uncolored"
    }

    // Temporary array to mark colors used by neighbors. Safely sized to V.
    bool* available = (bool*)malloc(V * sizeof(bool));
    if (!available) {
         fprintf(stderr, "Failed to allocate 'available' array\n");
         exit(1);
    }

    int max_color_used = -1;

    // Assign colors to vertices one by one
    for (int u = 0; u < V; u++) {
        // Reset available colors for the current vertex
        // The number of colors can't exceed max_color_used + 2
        for (int i = 0; i <= max_color_used + 1 && i < V; i++) {
            available[i] = true;
        }

        // For all adjacent vertices, mark their colors as unavailable
        for (int i = 0; i < degree[u]; i++) {
            int neighbor = adj[u][i];
            if (colors[neighbor] != -1) {
                if(colors[neighbor] < V) { // Safety check
                    available[colors[neighbor]] = false;
                }
            }
        }

        // Find the first available color
        int first_available_color = 0;
        while (first_available_color < V) {
            if (available[first_available_color]) {
                break;
            }
            first_available_color++;
        }

        colors[u] = first_available_color;

        if (first_available_color > max_color_used) {
            max_color_used = first_available_color;
        }
    }

    free(available);

    // The result is the total number of colors used (chromatic number guess)
    final_result = max_color_used + 1;
}

void cleanup() {
    free(colors);
    for (int i = 0; i < V; i++) {
        if (adj[i] != NULL) {
          free(adj[i]);
        }
    }
    free(adj);
    free(degree);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%d\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
