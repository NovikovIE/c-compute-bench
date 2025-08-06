#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

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

// --- BENCHMARK SPECIFIC CODE ---

#define MAX_SUPPORTED_DIM 8

typedef struct {
    int dim;
    uint32_t* vertices;
} Simplex;

typedef struct {
    // Parameters
    int num_initial_simplices;
    int max_dim;

    // Simplicial Complex Data
    Simplex** simplices_by_dim[MAX_SUPPORTED_DIM + 1];
    int count_by_dim[MAX_SUPPORTED_DIM + 1];
    int capacity_by_dim[MAX_SUPPORTED_DIM + 1];

    // Computational Result
    long long total_betti_sum;

} BenchmarkData;

static BenchmarkData g_data;

// Helper to compare vertices for sorting
int compare_vertices(const void* a, const void* b) {
    uint32_t va = *(const uint32_t*)a;
    uint32_t vb = *(const uint32_t*)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

// Helper to compare two simplices for duplicate checking
int compare_simplices(const Simplex* s1, const Simplex* s2) {
    if (s1->dim != s2->dim) return s1->dim - s2->dim;
    for (int i = 0; i <= s1->dim; ++i) {
        if (s1->vertices[i] != s2->vertices[i]) {
            return s1->vertices[i] - s2->vertices[i];
        }
    }
    return 0;
}

// Adds a simplex to the complex if it's not already present
void add_simplex(int dim, uint32_t* vertices) {
    Simplex* new_simplex = (Simplex*)malloc(sizeof(Simplex));
    if (!new_simplex) { perror("malloc simplex"); exit(1); }
    new_simplex->dim = dim;
    new_simplex->vertices = (uint32_t*)malloc(sizeof(uint32_t) * (dim + 1));
    if (!new_simplex->vertices) { perror("malloc vertices"); exit(1); }
    memcpy(new_simplex->vertices, vertices, sizeof(uint32_t) * (dim + 1));
    qsort(new_simplex->vertices, dim + 1, sizeof(uint32_t), compare_vertices);

    // Check for duplicates before adding
    for (int i = 0; i < g_data.count_by_dim[dim]; ++i) {
        if (compare_simplices(new_simplex, g_data.simplices_by_dim[dim][i]) == 0) {
            free(new_simplex->vertices);
            free(new_simplex);
            return; // It's a duplicate
        }
    }

    // Expand storage if necessary
    if (g_data.count_by_dim[dim] >= g_data.capacity_by_dim[dim]) {
        g_data.capacity_by_dim[dim] *= 2;
        g_data.simplices_by_dim[dim] = (Simplex**)realloc(g_data.simplices_by_dim[dim], sizeof(Simplex*) * g_data.capacity_by_dim[dim]);
        if (!g_data.simplices_by_dim[dim]) { perror("realloc"); exit(1); }
    }

    g_data.simplices_by_dim[dim][g_data.count_by_dim[dim]++] = new_simplex;
}

// Recursively generate all faces of a simplex and add them
void generate_and_add_faces(int dim, uint32_t* vertices) {
    add_simplex(dim, vertices);
    if (dim == 0) return;

    uint32_t* face_vertices = (uint32_t*)malloc(sizeof(uint32_t) * dim);
    if (!face_vertices) { perror("malloc face"); exit(1); }
    
    for (int i = 0; i <= dim; ++i) {
        int current_face_idx = 0;
        for (int j = 0; j <= dim; ++j) {
            if (i == j) continue;
            face_vertices[current_face_idx++] = vertices[j];
        }
        generate_and_add_faces(dim - 1, face_vertices);
    }
    free(face_vertices);
}

int compute_rank_gf2(uint8_t** matrix, int rows, int cols) {
    int rank = 0;
    int pivot_row = 0;
    for (int j = 0; j < cols && pivot_row < rows; ++j) {
        int i = pivot_row;
        while (i < rows && matrix[i][j] == 0) {
            i++;
        }

        if (i < rows) {
            uint8_t* temp = matrix[i];
            matrix[i] = matrix[pivot_row];
            matrix[pivot_row] = temp;

            for (int k = 0; k < rows; ++k) {
                if (k != pivot_row && matrix[k][j] == 1) {
                    for (int l = j; l < cols; ++l) {
                        matrix[k][l] ^= matrix[pivot_row][l];
                    }
                }
            }
            pivot_row++;
            rank++;
        }
    }
    return rank;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_simplices> <max_simplicial_dimension> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_initial_simplices = atoi(argv[1]);
    g_data.max_dim = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (g_data.max_dim > MAX_SUPPORTED_DIM) {
        fprintf(stderr, "Max dimension %d exceeds supported max %d\n", g_data.max_dim, MAX_SUPPORTED_DIM);
        exit(1);
    }

    for (int i = 0; i <= g_data.max_dim; ++i) {
        g_data.count_by_dim[i] = 0;
        g_data.capacity_by_dim[i] = 100;
        g_data.simplices_by_dim[i] = (Simplex**)malloc(sizeof(Simplex*) * g_data.capacity_by_dim[i]);
        if(!g_data.simplices_by_dim[i]) { perror("malloc"); exit(1); }
    }
    g_data.total_betti_sum = 0;

    int vertex_pool_size = g_data.num_initial_simplices * (g_data.max_dim + 1) / 2;
    if (vertex_pool_size < g_data.max_dim + 1) vertex_pool_size = g_data.max_dim + 1;
    
    uint32_t* vertex_pool = (uint32_t*)malloc(sizeof(uint32_t) * vertex_pool_size);
    if(!vertex_pool) { perror("malloc"); exit(1); }
    for (int i = 0; i < vertex_pool_size; ++i) vertex_pool[i] = i;

    uint32_t* current_simplex_vertices = (uint32_t*)malloc(sizeof(uint32_t) * (g_data.max_dim + 1));
    if(!current_simplex_vertices) { perror("malloc"); exit(1); }

    for (int i = 0; i < g_data.num_initial_simplices; ++i) {
        int dim = (mt_rand() % g_data.max_dim) + 1;
        
        for (int j = 0; j < vertex_pool_size; j++) { // Fisher-Yates shuffle
            int k = j + (mt_rand() % (vertex_pool_size - j));
            uint32_t temp = vertex_pool[j];
            vertex_pool[j] = vertex_pool[k];
            vertex_pool[k] = temp;
        }

        memcpy(current_simplex_vertices, vertex_pool, sizeof(uint32_t) * (dim + 1));
        generate_and_add_faces(dim, current_simplex_vertices);
    }

    free(vertex_pool);
    free(current_simplex_vertices);
}

void run_computation() {
    int* ranks = (int*)calloc(g_data.max_dim + 2, sizeof(int));
    if(!ranks) { perror("malloc ranks"); exit(1); }

    for (int k = 1; k <= g_data.max_dim + 1; ++k) {
        int dim_k_minus_1 = k - 1;
        int dim_k = k;
        if (dim_k > g_data.max_dim) continue;

        int rows = g_data.count_by_dim[dim_k_minus_1];
        int cols = g_data.count_by_dim[dim_k];

        if (rows == 0 || cols == 0) {
            ranks[k] = 0;
            continue;
        }

        // Allocate and zero-initialize boundary matrix
        uint8_t** boundary_matrix = (uint8_t**)malloc(sizeof(uint8_t*) * rows);
        if(!boundary_matrix) { perror("malloc matrix rows"); exit(1); }
        for (int i = 0; i < rows; ++i) {
            boundary_matrix[i] = (uint8_t*)calloc(cols, sizeof(uint8_t));
            if(!boundary_matrix[i]) { perror("calloc matrix cols"); exit(1); }
        }

        // Build boundary matrix D_k
        uint32_t* face_v = (uint32_t*)malloc(sizeof(uint32_t) * dim_k);
        for (int j = 0; j < cols; ++j) { // For each k-simplex (column)
            Simplex* s_k = g_data.simplices_by_dim[dim_k][j];
            for (int v_rem_idx = 0; v_rem_idx <= dim_k; ++v_rem_idx) { // Find its faces
                int current_v_idx = 0;
                for (int v_idx = 0; v_idx <= dim_k; ++v_idx) {
                    if (v_idx == v_rem_idx) continue;
                    face_v[current_v_idx++] = s_k->vertices[v_idx];
                }
                // Search for this face in the list of (k-1)-simplices
                for (int i = 0; i < rows; ++i) {
                    if (memcmp(face_v, g_data.simplices_by_dim[dim_k_minus_1][i]->vertices, sizeof(uint32_t) * dim_k) == 0) {
                        boundary_matrix[i][j] = 1; // Not using signs, working in GF(2)
                        break;
                    }
                }
            }
        }
        free(face_v);
        
        ranks[k] = compute_rank_gf2(boundary_matrix, rows, cols);

        for (int i = 0; i < rows; ++i) free(boundary_matrix[i]);
        free(boundary_matrix);
    }

    g_data.total_betti_sum = 0;
    // Betti 0
    long b0 = g_data.count_by_dim[0] - ranks[1];
    g_data.total_betti_sum += b0 > 0 ? b0 : 0;

    // Betti_k for k > 0
    for (int k = 1; k <= g_data.max_dim; ++k) {
        long bk = (long)g_data.count_by_dim[k] - ranks[k] - ranks[k + 1];
        g_data.total_betti_sum += bk > 0 ? bk : 0;
    }

    free(ranks);
}

void cleanup() {
    for (int i = 0; i <= g_data.max_dim; ++i) {
        for (int j = 0; j < g_data.count_by_dim[i]; ++j) {
            free(g_data.simplices_by_dim[i][j]->vertices);
            free(g_data.simplices_by_dim[i][j]);
        }
        free(g_data.simplices_by_dim[i]);
    }
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data.total_betti_sum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}