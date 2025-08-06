#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>
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

// --- BENCHMARK DATA & STRUCTURES ---
struct Edge {
    int to;
    int weight;
    struct Edge* next;
};

// Global pointers for benchmark data
int V;
int E;
int S; // Start Vertex
struct Edge** adj;
int* dist;
long long total_distance_sum; // Benchmark result

// --- MIN-HEAP IMPLEMENTATION FOR PRIORITY QUEUE (STATIC HELPERS) ---
typedef struct MinHeapNode {
    int v; // vertex number
    int d; // distance from source
} MinHeapNode;

typedef struct MinHeap {
    int size;
    int capacity;
    int* pos; // To check if a vertex is in the heap and for decreaseKey
    MinHeapNode** array;
} MinHeap;

static MinHeapNode* newMinHeapNode(int v, int d) {
    MinHeapNode* node = (MinHeapNode*)malloc(sizeof(MinHeapNode));
    node->v = v;
    node->d = d;
    return node;
}

static MinHeap* createMinHeap(int capacity) {
    MinHeap* minHeap = (MinHeap*)malloc(sizeof(MinHeap));
    minHeap->pos = (int*)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (MinHeapNode**)malloc(capacity * sizeof(MinHeapNode*));
    return minHeap;
}

static void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b) {
    MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

static void minHeapify(MinHeap* minHeap, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;

    if (left < minHeap->size && minHeap->array[left]->d < minHeap->array[smallest]->d)
        smallest = left;

    if (right < minHeap->size && minHeap->array[right]->d < minHeap->array[smallest]->d)
        smallest = right;

    if (smallest != idx) {
        MinHeapNode* smallestNode = minHeap->array[smallest];
        MinHeapNode* idxNode = minHeap->array[idx];

        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;

        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);
        minHeapify(minHeap, smallest);
    }
}

static MinHeapNode* extractMin(MinHeap* minHeap) {
    if (minHeap->size == 0) return NULL;

    MinHeapNode* root = minHeap->array[0];
    MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;

    minHeap->pos[root->v] = minHeap->size - 1;
    minHeap->pos[lastNode->v] = 0;

    --minHeap->size;
    minHeapify(minHeap, 0);

    return root;
}

static void decreaseKey(MinHeap* minHeap, int v, int d) {
    int i = minHeap->pos[v];
    minHeap->array[i]->d = d;

    while (i && minHeap->array[i]->d < minHeap->array[(i - 1) / 2]->d) {
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

static bool isInMinHeap(MinHeap* minHeap, int v) {
    return minHeap->pos[v] < minHeap->size;
}

// --- BENCHMARK FUNCTIONS ---

void add_edge(int u, int v, int weight) {
    struct Edge* new_edge = (struct Edge*)malloc(sizeof(struct Edge));
    new_edge->to = v;
    new_edge->weight = weight;
    new_edge->next = adj[u];
    adj[u] = new_edge;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_vertices num_edges start_vertex seed\n", argv[0]);
        exit(1);
    }

    V = atoi(argv[1]);
    E = atoi(argv[2]);
    S = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    
    if (V <= 0 || E < V || S < 0 || S >= V) {
        fprintf(stderr, "Invalid arguments: V > 0, E >= V, 0 <= S < V required.\n");
        exit(1);
    }

    mt_seed(seed);

    adj = (struct Edge**)calloc(V, sizeof(struct Edge*));

    long edges_to_add = E;

    // To guarantee a connected graph, we first generate a Hamiltonian cycle.
    // This connects all vertices and uses V edges.
    for (int i = 0; i < V - 1; ++i) {
        add_edge(i, i + 1, (mt_rand() % 100) + 1);
    }
    add_edge(V - 1, 0, (mt_rand() % 100) + 1);
    edges_to_add -= V;

    // Add remaining edges randomly
    long edges_added = 0;
    while (edges_added < edges_to_add) {
        int u = mt_rand() % V;
        int v = mt_rand() % V;
        if (u == v) continue; // Avoid self-loops
        add_edge(u, v, (mt_rand() % 100) + 1);
        edges_added++;
    }

    dist = (int*)malloc(V * sizeof(int));
}

void run_computation() {
    MinHeap* minHeap = createMinHeap(V);

    for (int v = 0; v < V; ++v) {
        dist[v] = INT_MAX;
        minHeap->array[v] = newMinHeapNode(v, dist[v]);
        minHeap->pos[v] = v;
    }
    
    dist[S] = 0;
    decreaseKey(minHeap, S, 0);

    minHeap->size = V;

    while (minHeap->size > 0) {
        MinHeapNode* minHeapNode = extractMin(minHeap);
        int u = minHeapNode->v;

        struct Edge* pCrawl = adj[u];
        while (pCrawl != NULL) {
            int v = pCrawl->to;
            if (isInMinHeap(minHeap, v) && dist[u] != INT_MAX && pCrawl->weight + dist[u] < dist[v]) {
                dist[v] = dist[u] + pCrawl->weight;
                decreaseKey(minHeap, v, dist[v]);
            }
            pCrawl = pCrawl->next;
        }
        free(minHeapNode); // Free node after extracting
    }
    
    // Free the heap structure itself.
    free(minHeap->pos);
    free(minHeap->array);
    free(minHeap);

    // Accumulate the result to prevent dead code elimination.
    total_distance_sum = 0;
    for (int i = 0; i < V; i++) {
        if (dist[i] != INT_MAX) {
            total_distance_sum += dist[i];
        }
    }
}

void cleanup() {
    for (int i = 0; i < V; ++i) {
        struct Edge* current = adj[i];
        while (current != NULL) {
            struct Edge* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(adj);
    free(dist);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", total_distance_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
