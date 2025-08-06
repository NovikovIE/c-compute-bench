#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Benchmark Globals and Setup ---

typedef struct {
    int num_nodes;
    int num_edges;
    int num_iterations;

    // HITS scores
    double *auth_scores;
    double *hub_scores;
    double *next_auth_scores;
    double *next_hub_scores;

    // Graph (CSR for out-links L and in-links L^T)
    int *out_row_ptr;
    int *out_col_indices;
    int *in_row_ptr;
    int *in_col_indices;

    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_nodes num_edges num_iterations seed\n", argv[0]);
        exit(1);
    }

    g_data.num_nodes = atoi(argv[1]);
    g_data.num_edges = atoi(argv[2]);
    g_data.num_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate memory for score arrays
    g_data.auth_scores = (double*)malloc(g_data.num_nodes * sizeof(double));
    g_data.hub_scores = (double*)malloc(g_data.num_nodes * sizeof(double));
    g_data.next_auth_scores = (double*)malloc(g_data.num_nodes * sizeof(double));
    g_data.next_hub_scores = (double*)malloc(g_data.num_nodes * sizeof(double));

    // Allocate memory for graph CSR representation
    g_data.out_row_ptr = (int*)malloc((g_data.num_nodes + 1) * sizeof(int));
    g_data.out_col_indices = (int*)malloc(g_data.num_edges * sizeof(int));
    g_data.in_row_ptr = (int*)malloc((g_data.num_nodes + 1) * sizeof(int));
    g_data.in_col_indices = (int*)malloc(g_data.num_edges * sizeof(int));

    // Initialize scores to 1.0
    for (int i = 0; i < g_data.num_nodes; ++i) {
        g_data.auth_scores[i] = 1.0;
        g_data.hub_scores[i] = 1.0;
    }

    // --- Generate graph --- 
    // 1. Create a temporary edge list
    int* edge_src_list = (int*)malloc(g_data.num_edges * sizeof(int));
    int* edge_dst_list = (int*)malloc(g_data.num_edges * sizeof(int));
    for (int i = 0; i < g_data.num_edges; ++i) {
        edge_src_list[i] = mt_rand() % g_data.num_nodes;
        edge_dst_list[i] = mt_rand() % g_data.num_nodes;
    }

    // 2. Count degrees to build CSR structure
    int* out_degree = (int*)calloc(g_data.num_nodes, sizeof(int));
    int* in_degree = (int*)calloc(g_data.num_nodes, sizeof(int));
    for (int i = 0; i < g_data.num_edges; ++i) {
        out_degree[edge_src_list[i]]++;
        in_degree[edge_dst_list[i]]++;
    }

    // 3. Create row pointers (prefix sum of degrees)
    g_data.out_row_ptr[0] = 0;
    g_data.in_row_ptr[0] = 0;
    for (int i = 0; i < g_data.num_nodes; ++i) {
        g_data.out_row_ptr[i + 1] = g_data.out_row_ptr[i] + out_degree[i];
        g_data.in_row_ptr[i + 1] = g_data.in_row_ptr[i] + in_degree[i];
    }

    // 4. Fill column indices
    int* temp_out_degree = (int*)malloc(g_data.num_nodes * sizeof(int));
    int* temp_in_degree = (int*)malloc(g_data.num_nodes * sizeof(int));
    for(int i = 0; i < g_data.num_nodes; ++i) { // create copies for indexing
        temp_out_degree[i] = 0;
        temp_in_degree[i] = 0;
    }

    for (int i = 0; i < g_data.num_edges; ++i) {
        int src = edge_src_list[i];
        int dst = edge_dst_list[i];
        g_data.out_col_indices[g_data.out_row_ptr[src] + temp_out_degree[src]++] = dst;
        g_data.in_col_indices[g_data.in_row_ptr[dst] + temp_in_degree[dst]++] = src;
    }

    // Free temporary generation arrays
    free(edge_src_list);
    free(edge_dst_list);
    free(out_degree);
    free(in_degree);
    free(temp_out_degree);
    free(temp_in_degree);
}

void run_computation() {
    for (int iter = 0; iter < g_data.num_iterations; ++iter) {
        // 1. Authority Update Rule: a_new = L^T * h_old
        for (int i = 0; i < g_data.num_nodes; ++i) {
            g_data.next_auth_scores[i] = 0.0;
            for (int edge_idx = g_data.in_row_ptr[i]; edge_idx < g_data.in_row_ptr[i + 1]; ++edge_idx) {
                int incoming_node = g_data.in_col_indices[edge_idx];
                g_data.next_auth_scores[i] += g_data.hub_scores[incoming_node];
            }
        }

        // 2. Hub Update Rule: h_new = L * a_new
        for (int i = 0; i < g_data.num_nodes; ++i) {
            g_data.next_hub_scores[i] = 0.0;
            for (int edge_idx = g_data.out_row_ptr[i]; edge_idx < g_data.out_row_ptr[i + 1]; ++edge_idx) {
                int outgoing_node = g_data.out_col_indices[edge_idx];
                g_data.next_hub_scores[i] += g_data.next_auth_scores[outgoing_node];
            }
        }

        // 3. Normalization (L2-norm)
        double auth_norm_sq = 0.0;
        double hub_norm_sq = 0.0;
        for (int i = 0; i < g_data.num_nodes; ++i) {
            auth_norm_sq += g_data.next_auth_scores[i] * g_data.next_auth_scores[i];
            hub_norm_sq += g_data.next_hub_scores[i] * g_data.next_hub_scores[i];
        }
        double auth_norm = sqrt(auth_norm_sq);
        double hub_norm = sqrt(hub_norm_sq);

        // 4. Score update and copy to main arrays
        if (auth_norm > 1e-9 && hub_norm > 1e-9) {
            for (int i = 0; i < g_data.num_nodes; ++i) {
                g_data.auth_scores[i] = g_data.next_auth_scores[i] / auth_norm;
                g_data.hub_scores[i] = g_data.next_hub_scores[i] / hub_norm;
            }
        } else {
          // In case of a degenerate graph, scores may go to zero. Stop iterating.
          break;
        }
    }

    // Accumulate final result to prevent dead code elimination
    g_data.final_result = 0.0;
    for (int i = 0; i < g_data.num_nodes; i++) {
        g_data.final_result += g_data.auth_scores[i] + g_data.hub_scores[i];
    }
}

void cleanup() {
    free(g_data.auth_scores);
    free(g_data.hub_scores);
    free(g_data.next_auth_scores);
    free(g_data.next_hub_scores);
    free(g_data.out_row_ptr);
    free(g_data.out_col_indices);
    free(g_data.in_row_ptr);
    free(g_data.in_col_indices);
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
    printf("%f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
