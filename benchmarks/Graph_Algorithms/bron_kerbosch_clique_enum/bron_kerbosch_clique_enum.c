#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

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

// --- BENCHMARK DATA AND GLOBALS ---
int N; // Number of vertices
double DENSITY; // Graph density
bool** adj_matrix; // Adjacency matrix representation of the graph
int total_cliques; // Accumulated result

void bron_kerbosch_internal(bool* R, bool* P, bool* X);

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <graph_density> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    DENSITY = atof(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (N <= 0) {
        fprintf(stderr, "FATAL: Number of vertices must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    adj_matrix = (bool **)malloc(N * sizeof(bool *));
    if (adj_matrix == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for adjacency matrix rows.\n");
        exit(1);
    }
    for (int i = 0; i < N; i++) {
        adj_matrix[i] = (bool *)calloc(N, sizeof(bool));
        if (adj_matrix[i] == NULL) {
            fprintf(stderr, "FATAL: Memory allocation failed for adjacency matrix columns.\n");
            exit(1);
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if ((double)mt_rand() / (double)UINT32_MAX < DENSITY) {
                adj_matrix[i][j] = adj_matrix[j][i] = true;
            }
        }
    }

    total_cliques = 0;
}

void cleanup() {
    for (int i = 0; i < N; i++) {
        free(adj_matrix[i]);
    }
    free(adj_matrix);
}


// Recursive helper for Bron-Kerbosch algorithm with pivoting.
// R: The set of vertices in the current clique.
// P: The set of candidate vertices to extend the clique.
// X: The set of vertices already processed and not to be used.
void bron_kerbosch_internal(bool* R, bool* P, bool* X) {
    bool p_is_empty = true;
    for (int i = 0; i < N; ++i) {
        if (P[i]) {
            p_is_empty = false;
            break;
        }
    }

    if (p_is_empty) {
        bool x_is_empty = true;
        for (int i = 0; i < N; ++i) {
            if (X[i]) {
                x_is_empty = false;
                break;
            }
        }
        if (x_is_empty) {
            total_cliques++;
        }
        return;
    }

    // Choose a pivot 'u' from P U X to maximize |P intersect N(u)|
    int pivot = -1;
    int max_neighbors = -1;

    for (int i = 0; i < N; i++) {
        if (P[i] || X[i]) {
            int count = 0;
            for (int j = 0; j < N; j++) {
                if (P[j] && adj_matrix[i][j]) {
                    count++;
                }
            }
            if (count > max_neighbors) {
                max_neighbors = count;
                pivot = i;
            }
        }
    }
    if (pivot == -1) return; // Should not happen if P is not empty

    // Create a temporary copy of P to iterate over vertices not connected to the pivot
    bool P_copy_iter[N];
    for (int i = 0; i < N; i++) {
        P_copy_iter[i] = P[i] && !adj_matrix[pivot][i];
    }
    
    for (int v = 0; v < N; v++) {
        if (P_copy_iter[v]) {
            bool R_new[N], P_new[N], X_new[N];

            for (int i = 0; i < N; i++) R_new[i] = R[i];
            R_new[v] = true;

            for (int i = 0; i < N; i++) P_new[i] = P[i] && adj_matrix[v][i];
            for (int i = 0; i < N; i++) X_new[i] = X[i] && adj_matrix[v][i];
            
            bron_kerbosch_internal(R_new, P_new, X_new);

            P[v] = false;
            X[v] = true;
        }
    }
}

void run_computation() {
    bool* R = (bool*)calloc(N, sizeof(bool));
    bool* P = (bool*)malloc(N * sizeof(bool));
    bool* X = (bool*)calloc(N, sizeof(bool));

    if (!R || !P || !X) {
        fprintf(stderr, "FATAL: Memory allocation failed for R, P, X sets.\n");
        exit(1);
    }

    for (int i = 0; i < N; i++) {
        P[i] = true; // Initially, all vertices are candidates
    }

    bron_kerbosch_internal(R, P, X);

    free(R);
    free(P);
    free(X);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", total_cliques);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}