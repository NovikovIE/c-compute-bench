#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- MERSENNE TWISTER (verbatim) ---
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

// --- BENCHMARK DATA AND GLOBALS ---
// Benchmark: Enumerate Spanning Trees using Kirchhoff's Matrix Tree Theorem.
// The number of spanning trees is the determinant of any cofactor of the graph's Laplacian matrix.
// This benchmark generates a random graph, constructs a Laplacian cofactor,
// and computes its determinant using Gaussian elimination.

typedef struct {
    int n; // Size of the matrix (num_vertices - 1)
    double **laplacian_cofactor;
} BenchmarkData;

static BenchmarkData g_data = {0, NULL};
static long long g_result = 0;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_edges> <seed>\n", argv[0]);
        exit(1);
    }

    int num_vertices = atoi(argv[1]);
    int num_edges = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    if (num_vertices <= 1) {
        fprintf(stderr, "Error: Number of vertices must be greater than 1.\n");
        exit(1);
    }
    long long max_edges = (long long)num_vertices * (num_vertices - 1) / 2;
    if (num_edges < num_vertices - 1 || num_edges > max_edges) {
        fprintf(stderr, "Error: Number of edges must be between %d and %lld for a connected graph.\n", num_vertices - 1, max_edges);
        exit(1);
    }

    mt_seed(seed);

    int **adj = (int **)calloc(num_vertices, sizeof(int *));
    if (!adj) { exit(1); }
    for (int i = 0; i < num_vertices; ++i) {
        adj[i] = (int *)calloc(num_vertices, sizeof(int));
        if (!adj[i]) { exit(1); }
    }

    // Generate a random connected graph
    // 1. Create a random permutation for a path to ensure connectivity
    int *perm = (int *)malloc(num_vertices * sizeof(int));
    if (!perm) { exit(1); }
    for (int i = 0; i < num_vertices; ++i) perm[i] = i;

    for (int i = num_vertices - 1; i > 0; --i) {
        int j = mt_rand() % (i + 1);
        int temp = perm[i];
        perm[i] = perm[j];
        perm[j] = temp;
    }

    int edge_count = 0;
    for (int i = 0; i < num_vertices - 1; ++i) {
       adj[perm[i]][perm[i+1]] = 1;
       adj[perm[i+1]][perm[i]] = 1;
       edge_count++;
    }
    free(perm);

    // 2. Add remaining edges randomly
    while (edge_count < num_edges) {
        int u = mt_rand() % num_vertices;
        int v = mt_rand() % num_vertices;
        if (u != v && adj[u][v] == 0) {
            adj[u][v] = 1;
            adj[v][u] = 1;
            edge_count++;
        }
    }

    // 3. Construct the Laplacian cofactor (removing row 0 and column 0)
    g_data.n = num_vertices - 1;
    g_data.laplacian_cofactor = (double **)malloc(g_data.n * sizeof(double *));
    if(!g_data.laplacian_cofactor) { exit(1); }
    for (int i = 0; i < g_data.n; ++i) {
        g_data.laplacian_cofactor[i] = (double *)malloc(g_data.n * sizeof(double));
        if(!g_data.laplacian_cofactor[i]) { exit(1); }
    }

    for (int i = 0; i < g_data.n; ++i) {
        int v_i = i + 1;
        int degree = 0;
        for (int k = 0; k < num_vertices; ++k) {
            degree += adj[v_i][k];
        }
        g_data.laplacian_cofactor[i][i] = (double)degree;
        for (int j = 0; j < g_data.n; ++j) {
            if (i == j) continue;
            int v_j = j + 1;
            g_data.laplacian_cofactor[i][j] = (double)-adj[v_i][v_j];
        }
    }

    // Free temporary adjacency matrix
    for (int i = 0; i < num_vertices; ++i) {
        free(adj[i]);
    }
    free(adj);
}

// Compute determinant using Gaussian elimination with partial pivoting.
// Modifies the matrix in-place.
static double determinant_in_place(double **matrix, int size) {
    double det = 1.0;
    for (int i = 0; i < size; i++) {
        // Find pivot row
        int pivot = i;
        for (int j = i + 1; j < size; j++) {
            if (fabs(matrix[j][i]) > fabs(matrix[pivot][i])) {
                pivot = j;
            }
        }
        // Swap rows if necessary
        if (pivot != i) {
            double *temp = matrix[i];
            matrix[i] = matrix[pivot];
            matrix[pivot] = temp;
            det *= -1.0;
        }

        double pivot_val = matrix[i][i];
        if (fabs(pivot_val) < 1e-12) {
            return 0.0; // Matrix is singular
        }

        det *= pivot_val;

        for (int j = i + 1; j < size; j++) {
            double factor = matrix[j][i] / pivot_val;
            for (int k = i; k < size; k++) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
    }
    return det;
}

void run_computation() {
    double det = determinant_in_place(g_data.laplacian_cofactor, g_data.n);
    g_result = (long long)(det + 0.5); // Round to nearest integer
}

void cleanup() {
    if (g_data.laplacian_cofactor) {
        for (int i = 0; i < g_data.n; i++) {
            free(g_data.laplacian_cofactor[i]);
        }
        free(g_data.laplacian_cofactor);
    }
    g_data.laplacian_cofactor = NULL;
    g_data.n = 0;
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
    printf("%lld\n", g_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
