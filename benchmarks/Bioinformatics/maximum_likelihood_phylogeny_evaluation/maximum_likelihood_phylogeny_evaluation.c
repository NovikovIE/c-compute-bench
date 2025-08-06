#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Number of states (DNA bases: A, C, G, T)
#define NUM_STATES 4

// Global benchmark data structure
typedef struct {
    // Parameters
    int num_taxa;
    int sequence_length;
    int num_iterations;
    int num_nodes;

    // Data
    char **sequences;
    double **likelihood_vectors;
    int (*children)[2]; // Tree topology: maps internal node index to child indices
    double p_matrix[NUM_STATES][NUM_STATES]; // Substitution probability matrix

    // Result
    double final_log_likelihood;
} BenchmarkData;

BenchmarkData g_data;

// Helper to convert DNA base character to an integer index
static inline int base_to_int(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 0; // Should not happen
    }
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_taxa sequence_length num_iterations seed\n", argv[0]);
        exit(1);
    }

    g_data.num_taxa = atoi(argv[1]);
    g_data.sequence_length = atoi(argv[2]);
    g_data.num_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (g_data.num_taxa < 2) {
        fprintf(stderr, "num_taxa must be at least 2.\n");
        exit(1);
    }

    g_data.num_nodes = 2 * g_data.num_taxa - 1;
    g_data.final_log_likelihood = 0.0;
    
    // Allocate memory for sequences (as a contiguous block)
    g_data.sequences = (char **)malloc(g_data.num_taxa * sizeof(char *));
    if (!g_data.sequences) { perror("malloc sequences"); exit(1); }
    char *seq_data = (char *)malloc((size_t)g_data.num_taxa * g_data.sequence_length * sizeof(char));
    if (!seq_data) { perror("malloc seq_data"); exit(1); }
    for (int i = 0; i < g_data.num_taxa; ++i) {
        g_data.sequences[i] = &seq_data[ (size_t)i * g_data.sequence_length];
    }

    // Generate random DNA sequences
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < g_data.num_taxa; ++i) {
        for (int j = 0; j < g_data.sequence_length; ++j) {
            g_data.sequences[i][j] = bases[mt_rand() % NUM_STATES];
        }
    }
    
    // Allocate memory for likelihood vectors (contiguous for cache efficiency)
    g_data.likelihood_vectors = (double **)malloc(g_data.num_nodes * sizeof(double *));
    if (!g_data.likelihood_vectors) { perror("malloc likelihood_vectors"); exit(1); }
    double *lv_data = (double *)malloc((size_t)g_data.num_nodes * NUM_STATES * sizeof(double));
    if (!lv_data) { perror("malloc lv_data"); exit(1); }
    for (int i = 0; i < g_data.num_nodes; ++i) {
        g_data.likelihood_vectors[i] = &lv_data[i * NUM_STATES];
    }
    
    // Allocate and define the tree structure (a simple caterpillar tree)
    // Internal nodes are indexed num_taxa to num_nodes-1
    // children[i] stores children of node (num_taxa + i)
    g_data.children = malloc((g_data.num_taxa - 1) * sizeof(*g_data.children));
    if(!g_data.children) { perror("malloc children"); exit(1); }
    
    // Base of the caterpillar tree
    g_data.children[0][0] = 0; // Leaf 0
    g_data.children[0][1] = 1; // Leaf 1

    // Extend the caterpillar
    for (int i = 1; i < g_data.num_taxa - 1; ++i) {
        g_data.children[i][0] = g_data.num_taxa + i - 1; // Previous internal node
        g_data.children[i][1] = i + 1; // Next leaf
    }

    // Pre-calculate Jukes-Cantor P(t) matrix for a fixed branch length t=0.1
    double t = 0.1;
    double p_diag = 0.25 + 0.75 * exp(-4.0/3.0 * t);
    double p_offdiag = 0.25 - 0.25 * exp(-4.0/3.0 * t);
    for (int i = 0; i < NUM_STATES; ++i) {
        for (int j = 0; j < NUM_STATES; ++j) {
            g_data.p_matrix[i][j] = (i == j) ? p_diag : p_offdiag;
        }
    }
}

void run_computation() {
    double total_log_likelihood = 0.0;
    const double base_freq = 0.25; // Equal base frequencies for JC69 model

    for (int iter = 0; iter < g_data.num_iterations; ++iter) {
        double iter_log_likelihood = 0.0;
        
        for (int site = 0; site < g_data.sequence_length; ++site) {
            // 1. Initialize likelihoods at the tips (leaves) of the tree
            for (int i = 0; i < g_data.num_taxa; ++i) {
                int state = base_to_int(g_data.sequences[i][site]);
                for (int j = 0; j < NUM_STATES; ++j) {
                    g_data.likelihood_vectors[i][j] = (j == state) ? 1.0 : 0.0;
                }
            }

            // 2. Compute likelihoods for internal nodes using post-order traversal
            for (int i = 0; i < g_data.num_taxa - 1; ++i) {
                int parent_node = g_data.num_taxa + i;
                int child1_node = g_data.children[i][0];
                int child2_node = g_data.children[i][1];

                for (int p_state = 0; p_state < NUM_STATES; ++p_state) {
                    double sum_c1 = 0.0;
                    double sum_c2 = 0.0;
                    for (int c_state = 0; c_state < NUM_STATES; ++c_state) {
                        sum_c1 += g_data.p_matrix[p_state][c_state] * g_data.likelihood_vectors[child1_node][c_state];
                        sum_c2 += g_data.p_matrix[p_state][c_state] * g_data.likelihood_vectors[child2_node][c_state];
                    }
                    g_data.likelihood_vectors[parent_node][p_state] = sum_c1 * sum_c2;
                }
            }

            // 3. Compute the final likelihood for the site at the root
            int root_node = g_data.num_nodes - 1;
            double site_likelihood = 0.0;
            for (int i = 0; i < NUM_STATES; ++i) {
                site_likelihood += base_freq * g_data.likelihood_vectors[root_node][i];
            }
            
            if (site_likelihood > 0) {
                iter_log_likelihood += log(site_likelihood);
            }
        }
        total_log_likelihood += iter_log_likelihood;
    }
    g_data.final_log_likelihood = total_log_likelihood;
}

void cleanup() {
    free(g_data.sequences[0]); // Free contiguous block
    free(g_data.sequences);
    
    free(g_data.likelihood_vectors[0]); // Free contiguous block
    free(g_data.likelihood_vectors);
    
    free(g_data.children);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_log_likelihood);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
