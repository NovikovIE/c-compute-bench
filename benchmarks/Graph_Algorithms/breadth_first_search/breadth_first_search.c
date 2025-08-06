#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- Benchmark Data Structures ---
typedef struct AdjListNode {
    int dest;
    struct AdjListNode* next;
} AdjListNode;

typedef struct {
    AdjListNode* head;
} AdjList;

typedef struct {
    int v;
    AdjList* array;
} Graph;

// Global variables to hold benchmark data and parameters
int NUM_VERTICES;
int NUM_EDGES;
int START_VERTEX;
Graph* graph;
int* distances;
int* queue;
long long final_result; // Accumulated result

// Helper to add an edge to an undirected graph
void add_edge(int src, int dest) {
    AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
    if (!newNode) exit(1);
    newNode->dest = dest;
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <start_vertex> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_VERTICES = atoi(argv[1]);
    NUM_EDGES = atoi(argv[2]);
    START_VERTEX = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (NUM_VERTICES <= 0 || NUM_EDGES < 0 || START_VERTEX < 0 || START_VERTEX >= NUM_VERTICES) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate graph structure
    graph = (Graph*)malloc(sizeof(Graph));
    if (!graph) exit(1);
    graph->v = NUM_VERTICES;
    graph->array = (AdjList*)malloc(NUM_VERTICES * sizeof(AdjList));
    if (!graph->array) exit(1);

    // Initialize adjacency lists
    for (int i = 0; i < NUM_VERTICES; ++i) {
        graph->array[i].head = NULL;
    }

    // Add edges to create a random undirected graph
    for (int i = 0; i < NUM_EDGES; ++i) {
        int u, v;
        do {
            u = mt_rand() % NUM_VERTICES;
            v = mt_rand() % NUM_VERTICES;
        } while (u == v); // Avoid self-loops
        add_edge(u, v);
        add_edge(v, u);
    }

    // Allocate memory for BFS data
    distances = (int*)malloc(NUM_VERTICES * sizeof(int));
    if (!distances) exit(1);
    queue = (int*)malloc(NUM_VERTICES * sizeof(int));
    if (!queue) exit(1);

    final_result = 0;
}

void run_computation() {
    // Initialize distances array to -1 (representing unvisited)
    for (int i = 0; i < NUM_VERTICES; ++i) {
        distances[i] = -1;
    }

    // Initialize a simple array-based queue
    int queue_front = 0, queue_rear = 0;

    // Start BFS from the START_VERTEX
    distances[START_VERTEX] = 0;
    queue[queue_rear++] = START_VERTEX;

    while (queue_front < queue_rear) {
        // Dequeue a vertex
        int u = queue[queue_front++];

        // Get all adjacent vertices of the dequeued vertex u
        AdjListNode* pCrawl = graph->array[u].head;
        while (pCrawl) {
            int v = pCrawl->dest;
            if (distances[v] == -1) {
                distances[v] = distances[u] + 1;
                queue[queue_rear++] = v;
            }
            pCrawl = pCrawl->next;
        }
    }

    // To prevent dead code elimination, compute a checksum of all distances
    long long checksum = 0;
    for (int i = 0; i < NUM_VERTICES; ++i) {
        if (distances[i] != -1) {
            checksum += distances[i];
        }
    }
    final_result = checksum;
}

void cleanup() {
    for (int i = 0; i < NUM_VERTICES; ++i) {
        AdjListNode* pCrawl = graph->array[i].head;
        while (pCrawl != NULL) {
            AdjListNode* temp = pCrawl;
            pCrawl = pCrawl->next;
            free(temp);
        }
    }
    free(graph->array);
    free(graph);
    free(distances);
    free(queue);
}


// --- Main --- 

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
