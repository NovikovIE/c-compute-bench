#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND GLOBALS ---

#define MAX_SUPPORTED_DIMENSION 3 // Hard limit for this benchmark implementation

typedef struct {
    int num_points;
    int ambient_dimension;
    float epsilon_radius;
    int max_simplicial_dimension;
    uint32_t seed;

    float* points; // Size: num_points * ambient_dimension
    long long* simplex_counts; // Size: max_simplicial_dimension + 1
    char** adj_matrix; // Adjacency matrix for 1-skeleton graph

    long long final_result; // Accumulated result to prevent dead code elimination
} BenchmarkData;

static BenchmarkData g_data;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_points ambient_dimension epsilon_radius max_simplicial_dimension seed\n", argv[0]);
        exit(1);
    }

    g_data.num_points = atoi(argv[1]);
    g_data.ambient_dimension = atoi(argv[2]);
    g_data.epsilon_radius = atof(argv[3]);
    g_data.max_simplicial_dimension = atoi(argv[4]);
    g_data.seed = (uint32_t)strtoul(argv[5], NULL, 10);

    if (g_data.num_points <= 0 || g_data.ambient_dimension <= 0 || g_data.epsilon_radius <= 0.0f || g_data.max_simplicial_dimension < 0) {
        fprintf(stderr, "FATAL: Invalid benchmark parameters.\n");
        exit(1);
    }

    if (g_data.max_simplicial_dimension > MAX_SUPPORTED_DIMENSION) {
        fprintf(stderr, "Warning: max_simplicial_dimension (%d) > max supported (%d). Clamping.\n", 
                g_data.max_simplicial_dimension, MAX_SUPPORTED_DIMENSION);
        g_data.max_simplicial_dimension = MAX_SUPPORTED_DIMENSION;
    }

    mt_seed(g_data.seed);

    // Allocate memory for points in a hypercube [0, 1]^D
    g_data.points = (float*)malloc(g_data.num_points * g_data.ambient_dimension * sizeof(float));
    if (!g_data.points) {
        fprintf(stderr, "FATAL: Failed to allocate memory for points.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_points * g_data.ambient_dimension; ++i) {
        g_data.points[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    // Allocate memory for simplex counts
    g_data.simplex_counts = (long long*)calloc(g_data.max_simplicial_dimension + 1, sizeof(long long));
    if (!g_data.simplex_counts) {
        fprintf(stderr, "FATAL: Failed to allocate memory for simplex_counts.\n");
        exit(1);
    }

    // Allocate adjacency matrix for the graph (1-skeleton)
    // This is used for efficient lookup of higher-order simplices
    if (g_data.max_simplicial_dimension > 1) {
        g_data.adj_matrix = (char**)malloc(g_data.num_points * sizeof(char*));
        if (!g_data.adj_matrix) {
             fprintf(stderr, "FATAL: Failed to allocate memory for adj_matrix rows.\n");
             exit(1);
        }
        for (int i = 0; i < g_data.num_points; ++i) {
            g_data.adj_matrix[i] = (char*)calloc(g_data.num_points, sizeof(char));
            if(!g_data.adj_matrix[i]) {
                fprintf(stderr, "FATAL: Failed to allocate memory for adj_matrix cols.\n");
                exit(1);
            }
        }
    }

    g_data.final_result = 0;
}

void run_computation() {
    // This function builds a Vietoris-Rips complex, a common approximation of the
    // Čech complex in Topological Data Analysis. It counts all k-simplices up to
    // the specified max_simplicial_dimension.
    // A k-simplex is a set of k+1 points if the distance between any two of 
    // these points is less than or equal to epsilon_radius.

    const int n = g_data.num_points;
    const int dim = g_data.ambient_dimension;
    const float epsilon_sq = g_data.epsilon_radius * g_data.epsilon_radius;

    // 0-Simplices (Vertices)
    if (g_data.max_simplicial_dimension >= 0) {
        g_data.simplex_counts[0] = n;
    }
    
    // 1-Simplices (Edges)
    if (g_data.max_simplicial_dimension >= 1) {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                float dist_sq = 0.0f;
                for (int d = 0; d < dim; ++d) {
                    float diff = g_data.points[i * dim + d] - g_data.points[j * dim + d];
                    dist_sq += diff * diff;
                }

                if (dist_sq <= epsilon_sq) {
                    g_data.simplex_counts[1]++;
                    if (g_data.max_simplicial_dimension > 1) {
                        g_data.adj_matrix[i][j] = g_data.adj_matrix[j][i] = 1;
                    }
                }
            }
        }
    }

    // 2-Simplices (Triangles)
    if (g_data.max_simplicial_dimension >= 2) {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (g_data.adj_matrix[i][j]) {
                    for (int k = j + 1; k < n; ++k) {
                        if (g_data.adj_matrix[i][k] && g_data.adj_matrix[j][k]) {
                            g_data.simplex_counts[2]++;
                        }
                    }
                }
            }
        }
    }

    // 3-Simplices (Tetrahedra)
    if (g_data.max_simplicial_dimension >= 3) {
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (!g_data.adj_matrix[i][j]) continue;
                for (int k = j + 1; k < n; ++k) {
                    if (!g_data.adj_matrix[i][k] || !g_data.adj_matrix[j][k]) continue;
                    for (int l = k + 1; l < n; ++l) {
                        if (g_data.adj_matrix[i][l] && g_data.adj_matrix[j][l] && g_data.adj_matrix[k][l]) {
                            g_data.simplex_counts[3]++;
                        }
                    }
                }
            }
        }
    }

    // Sum up all found simplices for the final result
    for (int i = 0; i <= g_data.max_simplicial_dimension; ++i) {
        g_data.final_result += g_data.simplex_counts[i];
    }
}

void cleanup() {
    free(g_data.points);
    free(g_data.simplex_counts);

    if (g_data.adj_matrix) {
        for (int i = 0; i < g_data.num_points; ++i) {
            free(g_data.adj_matrix[i]);
        }
        free(g_data.adj_matrix);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print final result to stdout
    printf("%lld\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
