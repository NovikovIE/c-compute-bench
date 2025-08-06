/*
 * BENCHMARK: Quadratic Assignment Problem (Brute-Force)
 * 
 * DESCRIPTION:
 * This benchmark solves the Quadratic Assignment Problem (QAP) using a brute-force
 * exhaustive search. The QAP is a fundamental combinatorial optimization problem
 * in operations research.
 * 
 * The problem is defined as follows: Given a set of n facilities and a set of n
 * locations, along with an n x n 'flow' matrix (flow[i][j] is the flow of materials
 * between facility i and facility j) and an n x n 'distance' matrix (dist[k][l] is
 * the distance between location k and location l), the goal is to find an assignment
 * of facilities to locations (a permutation) that minimizes the total cost.
 * 
 * The cost function for a permutation 'p' is:
 *   cost = sum_{i=0 to n-1} sum_{j=0 to n-1} flow[i][j] * dist[p[i]][p[j]]
 *   where p[i] is the location assigned to facility i.
 *
 * This implementation explores every possible permutation (n!) of assignments to 
 * find the one with the minimum cost. The complexity is O(n! * n^2), which is 
 * computationally feasible only for very small n.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <limits.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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

// --- Global Benchmark Data ---
int num_facilities;
long long min_cost;
int **flow_matrix; // Flow between facilities
int **dist_matrix; // Distance between locations

// --- Function Prototypes for Helpers ---
static void generate_permutations(int k, int* p);
static void calculate_and_update_cost(const int* p);
static void swap(int* a, int* b);

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_facilities> <seed>\n", argv[0]);
        exit(1);
    }

    num_facilities = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_facilities <= 0 || num_facilities > 12) { // 12! is huge, safeguard
        fprintf(stderr, "FATAL: num_facilities must be between 1 and 12.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate matrices
    flow_matrix = (int **)malloc(num_facilities * sizeof(int *));
    dist_matrix = (int **)malloc(num_facilities * sizeof(int *));
    if (!flow_matrix || !dist_matrix) {
        fprintf(stderr, "FATAL: Memory allocation failed for matrices.\n");
        exit(1);
    }

    for (int i = 0; i < num_facilities; i++) {
        flow_matrix[i] = (int *)malloc(num_facilities * sizeof(int));
        dist_matrix[i] = (int *)malloc(num_facilities * sizeof(int));
        if (!flow_matrix[i] || !dist_matrix[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for matrix rows.\n");
            exit(1);
        }
    }

    // Populate matrices with random symmetric data
    for (int i = 0; i < num_facilities; i++) {
        for (int j = i; j < num_facilities; j++) {
            if (i == j) {
                flow_matrix[i][j] = 0;
                dist_matrix[i][j] = 0;
            } else {
                int flow_val = mt_rand() % 100 + 1;
                int dist_val = mt_rand() % 100 + 1;
                flow_matrix[i][j] = flow_matrix[j][i] = flow_val;
                dist_matrix[i][j] = dist_matrix[j][i] = dist_val;
            }
        }
    }
}

void run_computation() {
    min_cost = LLONG_MAX;

    int* p = (int*)malloc(num_facilities * sizeof(int));
    if (p == NULL) {
        fprintf(stderr, "FATAL: Failed to allocate permutation array.\n");
        exit(1);
    }
    for (int i = 0; i < num_facilities; i++) {
        p[i] = i;
    }

    generate_permutations(num_facilities, p);

    free(p);
}

void cleanup() {
    for (int i = 0; i < num_facilities; i++) {
        free(flow_matrix[i]);
        free(dist_matrix[i]);
    }
    free(flow_matrix);
    free(dist_matrix);
}

// --- Helper Functions for Computation ---

static void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

static void calculate_and_update_cost(const int* p) {
    long long current_cost = 0;
    for (int i = 0; i < num_facilities; i++) {
        for (int j = 0; j < num_facilities; j++) {
            current_cost += (long long)flow_matrix[i][j] * dist_matrix[p[i]][p[j]];
        }
    }
    if (current_cost < min_cost) {
        min_cost = current_cost;
    }
}

// Heap's algorithm for generating permutations
static void generate_permutations(int k, int* p) {
    if (k == 1) {
        calculate_and_update_cost(p);
    } else {
        generate_permutations(k - 1, p);
        for (int i = 0; i < k - 1; i++) {
            if (k % 2 == 0) {
                swap(&p[i], &p[k - 1]);
            } else {
                swap(&p[0], &p[k - 1]);
            }
            generate_permutations(k - 1, p);
        }
    }
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
    printf("%lld\n", min_cost);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
