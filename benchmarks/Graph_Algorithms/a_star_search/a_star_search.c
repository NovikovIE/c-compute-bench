#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Specific Code ---

// Parameters
int NUM_VERTICES;
int NUM_EDGES;
int START_VERTEX;
int END_VERTEX;

// Graph Data Structures
typedef struct {
    int to;       // Destination vertex
    double weight;  // Edge weight
} Edge;

typedef struct {
    double x, y;    // Coordinates for heuristic
    Edge* edges;
    int count;      // Number of edges
    int capacity;   // Allocated capacity for edges
} Vertex;

Vertex* graph; // Global graph representation

double final_path_length = -1.0; // Benchmark result

// Priority Queue for A* Search (Min-Heap)
typedef struct {
    int vertex_id;
    double f_score; // Priority value (g_score + heuristic)
} PQNode;

PQNode* pq_heap;
int pq_size;

void pq_swap(int i, int j) {
    PQNode temp = pq_heap[i];
    pq_heap[i] = pq_heap[j];
    pq_heap[j] = temp;
}

void pq_heapify_up(int index) {
    if (index == 0) return;
    int parent_index = (index - 1) / 2;
    if (pq_heap[index].f_score < pq_heap[parent_index].f_score) {
        pq_swap(index, parent_index);
        pq_heapify_up(parent_index);
    }
}

void pq_heapify_down(int index) {
    int left = 2 * index + 1;
    int right = 2 * index + 2;
    int smallest = index;

    if (left < pq_size && pq_heap[left].f_score < pq_heap[smallest].f_score) {
        smallest = left;
    }
    if (right < pq_size && pq_heap[right].f_score < pq_heap[smallest].f_score) {
        smallest = right;
    }

    if (smallest != index) {
        pq_swap(index, smallest);
        pq_heapify_down(smallest);
    }
}

void pq_push(int vertex_id, double f_score, int max_pq_size) {
    if (pq_size >= max_pq_size) {
        // This should not happen with reasonable allocation
        return;
    }
    pq_heap[pq_size].vertex_id = vertex_id;
    pq_heap[pq_size].f_score = f_score;
    pq_size++;
    pq_heapify_up(pq_size - 1);
}

PQNode pq_pop() {
    PQNode root = pq_heap[0];
    pq_heap[0] = pq_heap[pq_size - 1];
    pq_size--;
    pq_heapify_down(0);
    return root;
}

// Heuristic function (Euclidean distance)
double heuristic(int u, int v) {
    double dx = graph[u].x - graph[v].x;
    double dy = graph[u].y - graph[v].y;
    return sqrt(dx * dx + dy * dy);
}

void add_edge(int u, int v, double weight) {
    if (graph[u].count >= graph[u].capacity) {
        graph[u].capacity *= 2;
        graph[u].edges = (Edge*)realloc(graph[u].edges, graph[u].capacity * sizeof(Edge));
    }
    graph[u].edges[graph[u].count].to = v;
    graph[u].edges[graph[u].count].weight = weight;
    graph[u].count++;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <start_vertex> <end_vertex> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_VERTICES = atoi(argv[1]);
    NUM_EDGES = atoi(argv[2]);
    START_VERTEX = atoi(argv[3]);
    END_VERTEX = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if (START_VERTEX >= NUM_VERTICES || END_VERTEX >= NUM_VERTICES || START_VERTEX < 0 || END_VERTEX < 0) {
        fprintf(stderr, "Error: Start/End vertex out of bounds.\n");
        exit(1);
    }

    mt_seed(seed);

    graph = (Vertex*)malloc(NUM_VERTICES * sizeof(Vertex));
    for (int i = 0; i < NUM_VERTICES; ++i) {
        graph[i].x = (double)mt_rand() / (double)UINT32_MAX;
        graph[i].y = (double)mt_rand() / (double)UINT32_MAX;
        graph[i].count = 0;
        graph[i].capacity = 8; // Initial capacity
        graph[i].edges = (Edge*)malloc(graph[i].capacity * sizeof(Edge));
    }

    for (int i = 0; i < NUM_EDGES; ++i) {
        int u = mt_rand() % NUM_VERTICES;
        int v = mt_rand() % NUM_VERTICES;
        if (u == v) {
            i--; // try again
            continue;
        }
        double weight = heuristic(u, v); // Use Euclidean distance as weight
        add_edge(u, v, weight);
        add_edge(v, u, weight);
    }
}

void run_computation() {
    double* g_score = (double*)malloc(NUM_VERTICES * sizeof(double));
    int* came_from = (int*)malloc(NUM_VERTICES * sizeof(int));
    
    // The number of items in PQ can exceed NUM_VERTICES because of stale entries,
    // but a safe upper bound is roughly related to NUM_EDGES.
    int max_pq_size = NUM_EDGES + 1;
    pq_heap = (PQNode*)malloc(max_pq_size * sizeof(PQNode));
    pq_size = 0;

    for (int i = 0; i < NUM_VERTICES; ++i) {
        g_score[i] = DBL_MAX;
        came_from[i] = -1;
    }

    g_score[START_VERTEX] = 0.0;
    pq_push(START_VERTEX, heuristic(START_VERTEX, END_VERTEX), max_pq_size);

    while (pq_size > 0) {
        PQNode current_node = pq_pop();
        int u = current_node.vertex_id;

        // If we pop a node with a score worse than what we've already found for it,
        // it's a stale entry. We can ignore it.
        if (current_node.f_score > g_score[u] + heuristic(u, END_VERTEX)) {
             continue;
        }

        if (u == END_VERTEX) {
            final_path_length = g_score[u];
            break;
        }

        for (int i = 0; i < graph[u].count; ++i) {
            Edge* edge = &graph[u].edges[i];
            int v = edge->to;
            double tentative_g_score = g_score[u] + edge->weight;

            if (tentative_g_score < g_score[v]) {
                came_from[v] = u;
                g_score[v] = tentative_g_score;
                double f_score = g_score[v] + heuristic(v, END_VERTEX);
                pq_push(v, f_score, max_pq_size);
            }
        }
    }

    free(g_score);
    free(came_from);
    free(pq_heap);
}

void cleanup() {
    for (int i = 0; i < NUM_VERTICES; ++i) {
        free(graph[i].edges);
    }
    free(graph);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_path_length);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
