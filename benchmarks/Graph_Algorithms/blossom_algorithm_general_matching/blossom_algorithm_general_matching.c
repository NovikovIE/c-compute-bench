#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// Mersenne Twister (Do Not Modify - Included Verbatim)
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

// Benchmark parameters and data structures
int V; // Number of vertices
int E; // Number of edges

// Adjacency list representation of the graph
int **adj;
int *adj_degrees;
int *adj_data;

// Data structures for Blossom Algorithm
int *match;    // match[i] = vertex matched with i, or -1
int *q;        // queue for BFS
int *p;        // parent in the alternating forest
int *base;     // base of the blossom
int *d;        // d[i] = 0 for unvisited, 1 for even, -1 for odd
int *blossom_visited;
int final_result;

static int lca(int a, int b, int start_node) {
    static int t = 0;
    t++;
    while (1) {
        if (a != -1) {
            a = base[a];
            if (blossom_visited[a] == t) return a;
            blossom_visited[a] = t;
            a = (match[a] != -1) ? p[match[a]] : -1;
        }
        if (b != -1) {
            b = base[b];
            if (blossom_visited[b] == t) return b;
            blossom_visited[b] = t;
            b = (match[b] != -1) ? p[match[b]] : -1;
        }
    }
}

static void blossom_contract(int u, int v, int l, int* head, int* tail) {
    while (base[u] != l) {
        p[u] = v;
        v = match[u];
        base[u] = l;
        base[v] = l;
        if (d[v] != 1) {
            d[v] = 1;
            q[(*tail)++] = v;
        }
        u = p[v];
    }
}

static void find_path(int start_node) {
    for (int i = 0; i < V; ++i) {
        d[i] = 0; p[i] = -1; base[i] = i;
    }
    d[start_node] = 1;
    int head = 0, tail = 0;
    q[tail++] = start_node;

    while(head < tail) {
        int u = q[head++];
        for (int i = 0; i < adj_degrees[u]; ++i) {
            int v = adj[u][i];
            if (base[u] == base[v] || match[u] == v) continue;

            if (d[v] == 0) {
                d[v] = -1;
                p[v] = u;
                if(match[v] == -1) {
                    // Augment path
                    while(v != -1) {
                        int pu = p[v];
                        int ppv = match[pu];
                        match[v] = pu;
                        match[pu] = v;
                        v = ppv;
                    }
                    return;
                }
                int mv = match[v];
                d[mv] = 1;
                q[tail++] = mv;
            } else if (d[base[v]] == 1) {
                int b = lca(u, v, start_node);
                blossom_contract(u, v, b, &head, &tail);
                blossom_contract(v, u, b, &head, &tail);
            }
        }
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    V = atoi(argv[1]);
    E = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (V <= 0 || E < 0) {
        fprintf(stderr, "FATAL: num_vertices must be positive and num_edges non-negative.\n");
        exit(1);
    }
    
    mt_seed(seed);

    adj_degrees = (int *)calloc(V, sizeof(int));
    if (!adj_degrees) { fprintf(stderr, "Failed to allocate adj_degrees\n"); exit(1); }

    // First pass: determine degrees
    for (int i = 0; i < E; ++i) {
        int u = mt_rand() % V;
        int v = mt_rand() % V;
        if (u == v) {
            i--; continue;
        }
        adj_degrees[u]++;
        adj_degrees[v]++;
    }

    long long total_degree = 0;
    for (int i = 0; i < V; ++i) total_degree += adj_degrees[i];

    adj_data = (int *)malloc(total_degree * sizeof(int));
    adj = (int **)malloc(V * sizeof(int *));
    int *current_degrees = (int *)calloc(V, sizeof(int));

    if (!adj_data || !adj || !current_degrees) {
        fprintf(stderr, "FATAL: Memory allocation failed for graph.\n");
        exit(1);
    }

    int *ptr = adj_data;
    for (int i = 0; i < V; ++i) {
        adj[i] = ptr;
        ptr += adj_degrees[i];
    }
    
    mt_seed(seed); // Re-seed for same random sequence
    for (int i = 0; i < E; ++i) {
        int u = mt_rand() % V;
        int v = mt_rand() % V;
        if (u == v) { i--; continue; }
        adj[u][current_degrees[u]++] = v;
        adj[v][current_degrees[v]++] = u;
    }
    free(current_degrees);

    match = (int *)malloc(V * sizeof(int));
    q = (int *)malloc(V * sizeof(int));
    p = (int *)malloc(V * sizeof(int));
    base = (int *)malloc(V * sizeof(int));
    d = (int *)malloc(V * sizeof(int));
    blossom_visited = (int *)calloc(V, sizeof(int));
    if (!match || !q || !p || !base || !d || !blossom_visited) {
        fprintf(stderr, "FATAL: Memory allocation failed for algorithm data.\n");
        exit(1);
    }
    
    for (int i = 0; i < V; ++i) {
        match[i] = -1;
    }
    final_result = 0;
}

void run_computation() {
    for (int i = 0; i < V; ++i) {
        if (match[i] == -1) {
            find_path(i);
        }
    }
    
    int matching_size = 0;
    for (int i = 0; i < V; ++i) {
        if (match[i] > i) { // Count each edge once
            matching_size++;
        }
    }
    final_result = matching_size;
}

void cleanup() {
    free(adj_data);
    free(adj);
    free(adj_degrees);
    free(match);
    free(q);
    free(p);
    free(base);
    free(d);
    free(blossom_visited);
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
