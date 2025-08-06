#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <math.h> // For fabs

// --- Mersenne Twister (DO NOT MODIFY) ---
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


// --- Benchmark Globals ---
int V; // Number of vertices
int E; // Number of edges

// Graph representation (CSR for undirected graph)
int *adj;   // Adjacency list (size 2*E)
int *xadj;  // Index pointers into adj (size V+1)

// Betweenness centrality scores
double *centrality;

// Workspace for Brandes' algorithm (allocated once in setup)
int *queue;
int *stack;
int *dist;
double *sigma;
double *delta;

// Workspace for predecessor lists
// p_count: Stores predecessor counts for a given BFS run.
// p_starts: Stores start index for each vertex's predecessor list in p_buffer.
// p_buffer: A flat buffer to store all predecessors for a given BFS run.
int *p_count;
int *p_starts;
int *p_buffer;

double final_result_sum; // To prevent dead code elimination

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

// --- Setup Function ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    V = atoi(argv[1]);
    E = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (V <= 0 || E <= 0) {
        fprintf(stderr, "Number of vertices and edges must be positive.\n");
        exit(1);
    }
    if (E < V) {
        fprintf(stderr, "Number of edges must be at least number of vertices to ensure connectivity.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate all memory
    xadj = (int*)malloc((V + 1) * sizeof(int));
    adj = (int*)malloc(2 * E * sizeof(int));
    centrality = (double*)malloc(V * sizeof(double));
    
    // Workspace allocation
    queue = (int*)malloc(V * sizeof(int));
    stack = (int*)malloc(V * sizeof(int));
    dist = (int*)malloc(V * sizeof(int));
    sigma = (double*)malloc(V * sizeof(double));
    delta = (double*)malloc(V * sizeof(double));
    p_count = (int*)malloc(V * sizeof(int));
    p_starts = (int*)malloc((V + 1) * sizeof(int));
    p_buffer = (int*)malloc(E * sizeof(int)); // Upper bound for predecessors in a single BFS DAG

    if (!xadj || !adj || !centrality || !queue || !stack || !dist || !sigma || !delta || !p_count || !p_starts || !p_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    
    // --- Graph Generation ---
    int *edge_u = (int*)malloc(E * sizeof(int));
    int *edge_v = (int*)malloc(E * sizeof(int));
    int* degree = (int*)calloc(V, sizeof(int));
    
    if (!edge_u || !edge_v || !degree) {
        fprintf(stderr, "FATAL: Temp memory allocation failed for graph generation.\n");
        exit(1);
    }

    // 1. Ensure connectivity with a cycle
    for (int i = 0; i < V; ++i) {
        edge_u[i] = i;
        edge_v[i] = (i + 1) % V;
    }

    // 2. Add remaining random edges
    for (int i = V; i < E; ++i) {
        int u, v;
        do {
            u = mt_rand() % V;
            v = mt_rand() % V;
        } while (u == v); // No self-loops
        edge_u[i] = u;
        edge_v[i] = v;
    }

    // 3. Calculate degrees to build CSR structure
    for (int i = 0; i < E; ++i) {
        degree[edge_u[i]]++;
        degree[edge_v[i]]++;
    }

    // 4. Build xadj array (prefix sum of degrees)
    xadj[0] = 0;
    for (int i = 0; i < V; ++i) {
        xadj[i+1] = xadj[i] + degree[i];
    }

    // 5. Populate adj array
    int* temp_offsets = (int*)malloc(V * sizeof(int));
    memcpy(temp_offsets, xadj, V * sizeof(int));

    for (int i = 0; i < E; ++i) {
        int u = edge_u[i];
        int v = edge_v[i];
        adj[temp_offsets[u]++] = v;
        adj[temp_offsets[v]++] = u;
    }

    free(edge_u);
    free(edge_v);
    free(degree);
    free(temp_offsets);
}


// --- Computation Function ---
void run_computation() {
    for (int i = 0; i < V; ++i) {
        centrality[i] = 0.0;
    }

    // Brandes' Algorithm: iterate over each vertex as a source 's'
    for (int s = 0; s < V; ++s) {
        // --- Initialization for this source 's' ---
        int stack_top = 0;
        for (int i = 0; i < V; ++i) {
            dist[i] = -1;
            sigma[i] = 0.0;
            p_count[i] = 0; // Predecessor count
        }
        dist[s] = 0;
        sigma[s] = 1.0;
        
        int q_head = 0, q_tail = 0;
        queue[q_tail++] = s;

        // --- BFS to find shortest paths, sigmas, and predecessor counts ---
        while (q_head < q_tail) {
            int v = queue[q_head++];
            stack[stack_top++] = v;
            
            for (int i = xadj[v]; i < xadj[v+1]; ++i) {
                int w = adj[i];
                if (dist[w] < 0) {
                    dist[w] = dist[v] + 1;
                    queue[q_tail++] = w;
                }
                if (dist[w] == dist[v] + 1) {
                    sigma[w] += sigma[v];
                    p_count[w]++;
                }
            }
        }
        
        // --- Prepare predecessor lists from counts ---
        p_starts[0] = 0;
        for (int i = 0; i < V; ++i) {
            p_starts[i+1] = p_starts[i] + p_count[i];
            p_count[i] = 0; // Reset to be reused as an insertion counter
        }

        // Populate predecessors list (p_buffer) using discovered nodes
        for (int v_idx = 0; v_idx < q_tail; ++v_idx) {
            int v = queue[v_idx]; // Iterate in BFS order
            for (int i = xadj[v]; i < xadj[v+1]; ++i) {
                int w = adj[i];
                if (dist[w] == dist[v] + 1) {
                    int p_idx = p_starts[w] + p_count[w];
                    p_buffer[p_idx] = v;
                    p_count[w]++;
                }
            }
        }

        // --- Dependency accumulation phase ---
        for (int i = 0; i < V; ++i) {
            delta[i] = 0.0;
        }

        while (stack_top > 0) {
            int w = stack[--stack_top];
            
            int p_start_idx = p_starts[w];
            int p_end_idx = p_starts[w+1];

            for (int i = p_start_idx; i < p_end_idx; ++i) {
                int v = p_buffer[i];
                if (fabs(sigma[w]) > 1e-9) { // Avoid division by zero
                  delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w]);
                }
            }
            if (w != s) {
                centrality[w] += delta[w];
            }
        }
    }

    // --- Finalize and sum results ---
    final_result_sum = 0.0;
    for (int i = 0; i < V; ++i) {
        // For undirected graphs, Brandes' algorithm counts each path twice.
        centrality[i] /= 2.0;
        final_result_sum += centrality[i];
    }
}


// --- Cleanup Function ---
void cleanup() {
    free(xadj);
    free(adj);
    free(centrality);
    free(queue);
    free(stack);
    free(dist);
    free(sigma);
    free(delta);
    free(p_count);
    free(p_starts);
    free(p_buffer);
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

    // Print result to stdout
    printf("%f\n", final_result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
