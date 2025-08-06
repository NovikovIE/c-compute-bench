/*
 * BENCHMARK: Computational Topology - Betti Numbers from Homology Groups
 * 
 * DESCRIPTION:
 * This benchmark simulates a process common in Topological Data Analysis (TDA), which
 * involves building a simplicial complex from a point cloud and then analyzing its
 * topological properties. Specifically, it constructs a Vietoris-Rips (VR) complex 
 * from a randomly generated set of 2D points.
 *
 * The core computation involves enumerating simplices (vertices, edges, triangles, etc.)
 * up to a specified maximum dimension. A k-simplex is a set of k+1 points where every
 * pair of points is within a certain distance (epsilon) of each other. This is equivalent to
 * finding cliques in the graph where vertices are points and edges connect points closer
 * than epsilon.
 *
 * Instead of performing the full, complex boundary matrix reduction to find Betti numbers,
 * this benchmark computes the Euler characteristic (χ), a related topological invariant.
 * The Euler characteristic is the alternating sum of the number of simplices of each 
 * dimension: χ = c₀ - c₁ + c₂ - c₃ + ..., where cₖ is the number of k-simplices.
 * This calculation still requires the computationally intensive step of counting all
 * simplices, which scales factorially with the dimension and polynomially with the
 * number of points (O(N^(d+1)) for dimension d), making it a good CPU-bound workload.
 *
 * The program is hard-coded to support up to dimension 3 (tetrahedra).
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

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

// --- Benchmark Globals ---

// Parameters
static int NUM_POINTS;
static int MAX_SIMPLICIAL_DIMENSION;

// A small fixed distance threshold for building the Vietoris-Rips complex
static const float EPSILON = 0.15f;
static float EPSILON_SQUARED;

// Data Structures
static float (*points)[2] = NULL;       // 2D point cloud
static bool* adj_matrix = NULL;        // Adjacency matrix for the 1-skeleton graph
static long long* simplex_counts = NULL; // Array to store counts of k-simplices

// Final result to prevent dead code elimination
static long long final_result = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <max_simplicial_dimension> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_POINTS = atoi(argv[1]);
    MAX_SIMPLICIAL_DIMENSION = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_POINTS <= 0 || MAX_SIMPLICIAL_DIMENSION < 0) {
        fprintf(stderr, "FATAL: Invalid parameters.\n");
        exit(1);
    }

    if (MAX_SIMPLICIAL_DIMENSION > 3) {
        fprintf(stderr, "Warning: This benchmark is hard-coded for max_simplicial_dimension <= 3. Clamping value.\n");
        MAX_SIMPLICIAL_DIMENSION = 3;
    }

    mt_seed(seed);

    EPSILON_SQUARED = EPSILON * EPSILON;

    points = malloc(NUM_POINTS * sizeof(*points));
    adj_matrix = calloc((size_t)NUM_POINTS * NUM_POINTS, sizeof(bool));
    simplex_counts = calloc(MAX_SIMPLICIAL_DIMENSION + 1, sizeof(long long));

    if (!points || !adj_matrix || !simplex_counts) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate random 2D points in the unit square [0, 1] x [0, 1]
    for (int i = 0; i < NUM_POINTS; ++i) {
        points[i][0] = (float)mt_rand() / (float)UINT32_MAX;
        points[i][1] = (float)mt_rand() / (float)UINT32_MAX;
    }
}

void run_computation() {
    // Dimension 0: Points are 0-simplices
    simplex_counts[0] = NUM_POINTS;

    // Dimension 1: Edges are 1-simplices
    // Build the adjacency matrix for the 1-skeleton (graph of edges)
    if (MAX_SIMPLICIAL_DIMENSION >= 1) {
        for (int i = 0; i < NUM_POINTS; ++i) {
            for (int j = i + 1; j < NUM_POINTS; ++j) {
                float dx = points[i][0] - points[j][0];
                float dy = points[i][1] - points[j][1];
                float dist_sq = dx * dx + dy * dy;
                if (dist_sq <= EPSILON_SQUARED) {
                    adj_matrix[i * NUM_POINTS + j] = true;
                    adj_matrix[j * NUM_POINTS + i] = true;
                    simplex_counts[1]++;
                }
            }
        }
    }

    // Dimension 2: Triangles are 2-simplices
    // A triangle {i, j, k} exists if edges {i,j}, {i,k}, and {j,k} exist.
    if (MAX_SIMPLICIAL_DIMENSION >= 2) {
        for (int i = 0; i < NUM_POINTS; ++i) {
            for (int j = i + 1; j < NUM_POINTS; ++j) {
                if (adj_matrix[i * NUM_POINTS + j]) {
                    for (int k = j + 1; k < NUM_POINTS; ++k) {
                        if (adj_matrix[i * NUM_POINTS + k] && adj_matrix[j * NUM_POINTS + k]) {
                            simplex_counts[2]++;
                        }
                    }
                }
            }
        }
    }

    // Dimension 3: Tetrahedra are 3-simplices
    // A tetrahedron {i,j,k,l} exists if all 6 of its sub-edges exist.
    if (MAX_SIMPLICIAL_DIMENSION >= 3) {
        for (int i = 0; i < NUM_POINTS; ++i) {
            for (int j = i + 1; j < NUM_POINTS; ++j) {
                if (!adj_matrix[i * NUM_POINTS + j]) continue;
                for (int k = j + 1; k < NUM_POINTS; ++k) {
                    if (!adj_matrix[i * NUM_POINTS + k] || !adj_matrix[j * NUM_POINTS + k]) continue;
                    for (int l = k + 1; l < NUM_POINTS; ++l) {
                        if (adj_matrix[i * NUM_POINTS + l] && 
                            adj_matrix[j * NUM_POINTS + l] && 
                            adj_matrix[k * NUM_POINTS + l]) {
                            simplex_counts[3]++;
                        }
                    }
                }
            }
        }
    }

    // Compute the Euler Characteristic: χ = c₀ - c₁ + c₂ - c₃ + ...
    long long euler_char = 0;
    for (int i = 0; i <= MAX_SIMPLICIAL_DIMENSION; ++i) {
        if (i % 2 == 0) {
            euler_char += simplex_counts[i];
        } else {
            euler_char -= simplex_counts[i];
        }
    }
    final_result = euler_char;
}

void cleanup() {
    free(points);
    free(adj_matrix);
    free(simplex_counts);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}