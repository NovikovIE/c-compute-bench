#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Globals ---
typedef struct {
    int u;
    int v;
} Edge;

int num_vertices;
long long num_edges;
Edge* edges;
int* degrees;
unsigned long long total_degree_sum = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    num_edges = atoll(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_vertices <= 0 || num_edges <= 0) {
        fprintf(stderr, "Number of vertices and edges must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the edge list and degree counts
    edges = (Edge*)malloc(num_edges * sizeof(Edge));
    // Use calloc to initialize degree counts to zero
    degrees = (int*)calloc(num_vertices, sizeof(int));

    if (!edges || !degrees) {
        fprintf(stderr, "Failed to allocate memory for the graph.\n");
        exit(1);
    }

    // Generate a list of random edges for an undirected graph
    for (long long i = 0; i < num_edges; ++i) {
        int u = mt_rand() % num_vertices;
        int v = mt_rand() % num_vertices;
        // Avoid self-loops for a slightly more realistic graph
        while (u == v) {
            v = mt_rand() % num_vertices;
        }
        edges[i].u = u;
        edges[i].v = v;
    }
}

void run_computation() {
    // Calculate the degree for each vertex by iterating through the edge list
    // This is the core of the degree centrality calculation
    for (long long i = 0; i < num_edges; ++i) {
        degrees[edges[i].u]++;
        degrees[edges[i].v]++;
    }

    // To prevent dead-code elimination and have a single result,
    // sum all the calculated degrees.
    // By the handshaking lemma, this sum should equal 2 * num_edges.
    unsigned long long sum = 0;
    for (int i = 0; i < num_vertices; ++i) {
        sum += degrees[i];
    }
    total_degree_sum = sum;
}

void cleanup() {
    free(edges);
    free(degrees);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%llu\n", total_degree_sum);

    // Print the timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
