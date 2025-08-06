#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// Mersenne Twister Generator (Do Not Modify - Include This Verbatim):
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

// Benchmark-specific constants and data structures
#define DAMPING_FACTOR 0.85

typedef struct {
    int num_webpages;
    long num_links;
    int num_iterations;
    
    // Compressed Sparse Row (CSR) graph representation
    int* row_ptr;       // size: num_webpages + 1
    int* col_indices;   // size: num_links
    int* out_degree;    // size: num_webpages

    // PageRank data
    double* ranks;
    double* new_ranks;
    
    // Result
    double final_checksum;
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_webpages> <num_links> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_webpages = atoi(argv[1]);
    g_data.num_links = atol(argv[2]);
    g_data.num_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate memory
    g_data.out_degree = (int*)calloc(g_data.num_webpages, sizeof(int));
    g_data.row_ptr = (int*)malloc((g_data.num_webpages + 1) * sizeof(int));
    g_data.col_indices = (int*)malloc(g_data.num_links * sizeof(int));
    g_data.ranks = (double*)malloc(g_data.num_webpages * sizeof(double));
    g_data.new_ranks = (double*)malloc(g_data.num_webpages * sizeof(double));

    if (!g_data.out_degree || !g_data.row_ptr || !g_data.col_indices || !g_data.ranks || !g_data.new_ranks) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    
    // Generate graph in CSR format using a two-pass process.
    // Pass 1: Generate random links, store them temporarily, and count out-degrees.
    int* temp_from = (int*)malloc(g_data.num_links * sizeof(int));
    int* temp_to = (int*)malloc(g_data.num_links * sizeof(int));
    if (!temp_from || !temp_to) {
        fprintf(stderr, "FATAL: Temp memory allocation failed for link generation.\n");
        exit(1);
    }
    
    for (long i = 0; i < g_data.num_links; ++i) {
        int u = mt_rand() % g_data.num_webpages;
        int v;
        do {
            v = mt_rand() % g_data.num_webpages;
        } while (u == v); // No self-links
        
        temp_from[i] = u;
        temp_to[i] = v;
        g_data.out_degree[u]++;
    }

    // Pass 2: Populate CSR arrays (row_ptr and col_indices)
    g_data.row_ptr[0] = 0;
    for (int i = 0; i < g_data.num_webpages; ++i) {
        g_data.row_ptr[i+1] = g_data.row_ptr[i] + g_data.out_degree[i];
    }
    
    // Use a copy of row_ptr to track current insertion position for each row
    int* current_pos = (int*)malloc(g_data.num_webpages * sizeof(int));
    memcpy(current_pos, g_data.row_ptr, g_data.num_webpages * sizeof(int));

    for (long i = 0; i < g_data.num_links; ++i) {
        int u = temp_from[i];
        int v = temp_to[i];
        g_data.col_indices[current_pos[u]] = v;
        current_pos[u]++;
    }

    free(temp_from);
    free(temp_to);
    free(current_pos);

    // Initialize PageRank values
    double initial_rank = 1.0 / g_data.num_webpages;
    for (int i = 0; i < g_data.num_webpages; ++i) {
        g_data.ranks[i] = initial_rank;
    }
    
    g_data.final_checksum = 0.0;
}

void run_computation() {
    for (int iter = 0; iter < g_data.num_iterations; ++iter) {
        // Calculate contribution from dangling nodes (pages with no outgoing links)
        double dangling_sum = 0.0;
        for (int i = 0; i < g_data.num_webpages; ++i) {
            if (g_data.out_degree[i] == 0) {
                dangling_sum += g_data.ranks[i];
            }
        }

        // All pages receive a base rank and a share of the dangling rank
        double base_rank = (1.0 - DAMPING_FACTOR) / g_data.num_webpages;
        double dangling_contrib = DAMPING_FACTOR * dangling_sum / g_data.num_webpages;
        double total_base_contrib = base_rank + dangling_contrib;

        for (int i = 0; i < g_data.num_webpages; ++i) {
            g_data.new_ranks[i] = total_base_contrib;
        }
        
        // Distribute rank from non-dangling nodes to their neighbors
        for (int i = 0; i < g_data.num_webpages; ++i) {
            if (g_data.out_degree[i] > 0) {
                double rank_to_distribute = DAMPING_FACTOR * g_data.ranks[i] / g_data.out_degree[i];
                for (int j = g_data.row_ptr[i]; j < g_data.row_ptr[i+1]; ++j) {
                    int dest_page = g_data.col_indices[j];
                    g_data.new_ranks[dest_page] += rank_to_distribute;
                }
            }
        }

        // Swap rank pointers for the next iteration
        double* temp = g_data.ranks;
        g_data.ranks = g_data.new_ranks;
        g_data.new_ranks = temp;
    }

    // Calculate a checksum to prevent dead code elimination.
    // The sum should be very close to 1.0.
    double checksum = 0.0;
    for(int i = 0; i < g_data.num_webpages; ++i) {
        checksum += g_data.ranks[i];
    }
    g_data.final_checksum = checksum;
}

void cleanup() {
    free(g_data.row_ptr);
    free(g_data.col_indices);
    free(g_data.out_degree);
    free(g_data.ranks);
    free(g_data.new_ranks);
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
    printf("%f\n", g_data.final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
