#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

// --- Mersenne Twister (Do Not Modify) ---
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

// --- Benchmark Data Structures ---
// Parameters
int NUM_NODES;
int NUM_EDGES;
int NUM_FACILITIES;
int BREAK_TIME;

// Graph representation (Adjacency List)
typedef struct AdjListNode {
    int dest;
    int weight;
    struct AdjListNode* next;
} AdjListNode;

typedef struct {
    AdjListNode** head;
} Graph;

// Dijkstra's Min-Heap
typedef struct MinHeapNode {
    int v;
    int dist;
} MinHeapNode;

typedef struct MinHeap {
    int size;
    int capacity;
    int *pos;
    MinHeapNode **array;
} MinHeap;

// Global data pointers
Graph* graph;
int* facilities;
char* reachable_nodes; 
long long final_result = 0;

// --- Helper Function Prototypes ---
MinHeap* createMinHeap(int capacity);
void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b);
void minHeapify(MinHeap* minHeap, int idx);
int isEmpty(MinHeap* minHeap);
MinHeapNode* extractMin(MinHeap* minHeap);
void decreaseKey(MinHeap* minHeap, int v, int dist);
int isInMinHeap(MinHeap *minHeap, int v);
void addEdge(Graph* graph, int src, int dest, int weight);
void dijkstra(int src, int* dist, MinHeap* minHeap);

// --- Core Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_nodes> <num_edges> <num_facilities> <break_time_or_distance> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_NODES = atoi(argv[1]);
    NUM_EDGES = atoi(argv[2]);
    NUM_FACILITIES = atoi(argv[3]);
    BREAK_TIME = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    // Allocate graph
    graph = (Graph*)malloc(sizeof(Graph));
    graph->head = (AdjListNode**)malloc(NUM_NODES * sizeof(AdjListNode*));
    for (int i = 0; i < NUM_NODES; ++i) {
        graph->head[i] = NULL;
    }

    // Generate edges for the graph
    for (int i = 0; i < NUM_EDGES; ++i) {
        int u = mt_rand() % NUM_NODES;
        int v = mt_rand() % NUM_NODES;
        if (u == v) { // Avoid self-loops
             v = (v + 1) % NUM_NODES;
        }
        int weight = 1 + (mt_rand() % 100);
        addEdge(graph, u, v, weight);
        addEdge(graph, v, u, weight); // Undirected graph
    }

    // Allocate and generate facilities
    facilities = (int*)malloc(NUM_FACILITIES * sizeof(int));
    char* is_facility = (char*)calloc(NUM_NODES, sizeof(char));
    for (int i = 0; i < NUM_FACILITIES; ++i) {
        int facility_node;
        do {
            facility_node = mt_rand() % NUM_NODES;
        } while (is_facility[facility_node]); // Ensure unique facilities
        is_facility[facility_node] = 1;
        facilities[i] = facility_node;
    }
    free(is_facility);

    // Allocate result array
    reachable_nodes = (char*)calloc(NUM_NODES, sizeof(char));
}

void run_computation() {
    // Allocate data structures for Dijkstra once to isolate computation
    int* dist = (int*)malloc(NUM_NODES * sizeof(int));
    MinHeap* minHeap = createMinHeap(NUM_NODES);
    MinHeapNode* node_pool = (MinHeapNode*)malloc(NUM_NODES * sizeof(MinHeapNode));
    for (int i = 0; i < NUM_NODES; ++i) {
        minHeap->array[i] = &node_pool[i];
    }

    // Run Dijkstra from each facility
    for (int i = 0; i < NUM_FACILITIES; ++i) {
        int start_node = facilities[i];
        dijkstra(start_node, dist, minHeap);

        // Update the global set of reachable nodes
        for (int j = 0; j < NUM_NODES; ++j) {
            if (dist[j] <= BREAK_TIME) {
                reachable_nodes[j] = 1;
            }
        }
    }

    // Free Dijkstra's temporary structures
    free(dist);
    free(minHeap->pos);
    free(minHeap->array);
    free(minHeap);
    free(node_pool);
    
    // Calculate final result to prevent dead code elimination
    final_result = 0;
    for (int i = 0; i < NUM_NODES; ++i) {
        if (reachable_nodes[i]) {
            final_result += i;
        }
    }
}

void cleanup() {
    // Free graph
    for (int i = 0; i < NUM_NODES; ++i) {
        AdjListNode* current = graph->head[i];
        while (current != NULL) {
            AdjListNode* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(graph->head);
    free(graph);

    // Free other global allocations
    free(facilities);
    free(reachable_nodes);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}

// --- Helper Function Implementations ---
void addEdge(Graph* graph, int src, int dest, int weight) {
    AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
    newNode->dest = dest;
    newNode->weight = weight;
    newNode->next = graph->head[src];
    graph->head[src] = newNode;
}

MinHeap* createMinHeap(int capacity) {
    MinHeap* minHeap = (MinHeap*)malloc(sizeof(MinHeap));
    minHeap->pos = (int*)malloc(capacity * sizeof(int));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (MinHeapNode**)malloc(capacity * sizeof(MinHeapNode*));
    return minHeap;
}

void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b) {
    MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

void minHeapify(MinHeap* minHeap, int idx) {
    int smallest, left, right;
    smallest = idx;
    left = 2 * idx + 1;
    right = 2 * idx + 2;
    if (left < minHeap->size && minHeap->array[left]->dist < minHeap->array[smallest]->dist)
        smallest = left;
    if (right < minHeap->size && minHeap->array[right]->dist < minHeap->array[smallest]->dist)
        smallest = right;
    if (smallest != idx) {
        MinHeapNode *smallestNode = minHeap->array[smallest];
        MinHeapNode *idxNode = minHeap->array[idx];
        minHeap->pos[smallestNode->v] = idx;
        minHeap->pos[idxNode->v] = smallest;
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);
        minHeapify(minHeap, smallest);
    }
}

int isEmpty(MinHeap* minHeap) {
    return minHeap->size == 0;
}

MinHeapNode* extractMin(MinHeap* minHeap) {
    if (isEmpty(minHeap))
        return NULL;
    MinHeapNode* root = minHeap->array[0];
    MinHeapNode* lastNode = minHeap->array[minHeap->size - 1];
    minHeap->array[0] = lastNode;
    minHeap->pos[root->v] = minHeap->size -1;
    minHeap->pos[lastNode->v] = 0;
    --minHeap->size;
    minHeapify(minHeap, 0);
    return root;
}

void decreaseKey(MinHeap* minHeap, int v, int dist) {
    int i = minHeap->pos[v];
    minHeap->array[i]->dist = dist;
    while (i && minHeap->array[i]->dist < minHeap->array[(i - 1) / 2]->dist) {
        minHeap->pos[minHeap->array[i]->v] = (i - 1) / 2;
        minHeap->pos[minHeap->array[(i - 1) / 2]->v] = i;
        swapMinHeapNode(&minHeap->array[i], &minHeap->array[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

int isInMinHeap(MinHeap *minHeap, int v) {
   if (minHeap->pos[v] < minHeap->size)
     return 1;
   return 0;
}

void dijkstra(int src, int* dist, MinHeap* minHeap) {
    minHeap->size = NUM_NODES;
    for (int v = 0; v < NUM_NODES; ++v) {
        dist[v] = INT_MAX;
        minHeap->array[v]->v = v;
        minHeap->array[v]->dist = INT_MAX;
        minHeap->pos[v] = v;
    }

    dist[src] = 0;
    decreaseKey(minHeap, src, 0);

    while (!isEmpty(minHeap)) {
        MinHeapNode* minHeapNode = extractMin(minHeap);
        int u = minHeapNode->v;

        if (dist[u] == INT_MAX) break;

        AdjListNode* pCrawl = graph->head[u];
        while (pCrawl != NULL) {
            int v = pCrawl->dest;
            if (isInMinHeap(minHeap, v) && pCrawl->weight + dist[u] < dist[v]) {
                dist[v] = dist[u] + pCrawl->weight;
                decreaseKey(minHeap, v, dist[v]);
            }
            pCrawl = pCrawl->next;
        }
    }
}
