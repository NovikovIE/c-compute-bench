#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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
// --- End of Mersenne Twister ---

// --- Benchmark Data Structures ---
typedef struct AdjListNode {
    int dest;
    struct AdjListNode* next;
} AdjListNode;

typedef struct AdjList {
    AdjListNode *head;
} AdjList;

typedef struct Graph {
    int v;
    AdjList* array;
} Graph;

// --- Global Benchmark Variables ---
static int NUM_VERTICES;
static int NUM_EDGES;
static int START_VERTEX;
static Graph* graph;
static char* visited;
static long long final_result; // To store the result of the computation

void add_edge(Graph* g, int src, int dest) {
    AdjListNode* newNode = (AdjListNode*)malloc(sizeof(AdjListNode));
    if (!newNode) {
        fprintf(stderr, "Memory allocation failed for graph node\n");
        exit(1);
    }
    newNode->dest = dest;
    newNode->next = g->array[src].head;
    g->array[src].head = newNode;
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
        fprintf(stderr, "Invalid arguments.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate graph structure
    graph = (Graph*)malloc(sizeof(Graph));
    if (!graph) {
        fprintf(stderr, "Memory allocation failed for graph\n");
        exit(1);
    }
    graph->v = NUM_VERTICES;
    graph->array = (AdjList*)malloc(NUM_VERTICES * sizeof(AdjList));
    if (!graph->array) {
        fprintf(stderr, "Memory allocation failed for adjacency list\n");
        exit(1);
    }

    // Allocate visited array
    visited = (char*)malloc(NUM_VERTICES * sizeof(char));
    if (!visited) {
        fprintf(stderr, "Memory allocation failed for visited array\n");
        exit(1);
    }

    for (int i = 0; i < NUM_VERTICES; ++i) {
        graph->array[i].head = NULL;
    }

    // Generate edges
    for (int i = 0; i < NUM_EDGES; ++i) {
        int u = mt_rand() % NUM_VERTICES;
        int v = mt_rand() % NUM_VERTICES;
        add_edge(graph, u, v);
    }

    final_result = 0;
}

void run_computation() {
    long long total_visited_sum = 0;

    // Stack for iterative DFS
    int* stack = (int*)malloc(NUM_VERTICES * sizeof(int));
    if (!stack) {
        fprintf(stderr, "Memory allocation failed for DFS stack\n");
        exit(1);
    }
    int top = -1;

    // Reset visited array for this run
    for (int i = 0; i < NUM_VERTICES; i++) {
        visited[i] = 0;
    }

    // Start DFS from the start vertex
    stack[++top] = START_VERTEX;
    visited[START_VERTEX] = 1;

    while (top != -1) {
        int u = stack[top--];
        total_visited_sum += u;

        AdjListNode* pCrawl = graph->array[u].head;
        while (pCrawl) {
            int v = pCrawl->dest;
            if (!visited[v]) {
                visited[v] = 1;
                stack[++top] = v;
            }
            pCrawl = pCrawl->next;
        }
    }

    free(stack);
    final_result = total_visited_sum;
}

void cleanup() {
    for (int i = 0; i < NUM_VERTICES; i++) {
        AdjListNode* current = graph->array[i].head;
        while (current != NULL) {
            AdjListNode* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(graph->array);
    free(graph);
    free(visited);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
