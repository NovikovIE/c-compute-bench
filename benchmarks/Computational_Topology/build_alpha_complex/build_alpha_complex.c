#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) KERNEL DO NOT MODIFY ---
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
// --- END MT19937 KERNEL ---

// --- BENCHMARK DATA AND GLOBALS ---
typedef struct {
    double x, y, z;
} Point3D;

typedef struct {
    int num_points;
    int max_simplicial_dimension;
    double alpha; // This is used as epsilon_sq in Vietoris-Rips

    Point3D* points;
    int* current_simplex_indices;
    long long total_simplices_count;
} BenchmarkState;

BenchmarkState g_state;

// --- HELPER FUNCTIONS ---

// Calculate squared Euclidean distance between two points
double dist_sq(int p_idx1, int p_idx2) {
    Point3D p1 = g_state.points[p_idx1];
    Point3D p2 = g_state.points[p_idx2];
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return dx * dx + dy * dy + dz * dz;
}

// Recursive function to find cliques of size k (simplices of dimension k-1)
void find_simplices_recursive(int k, int level, int start_index, long long* count) {
    if (level == k) {
        (*count)++;
        return;
    }

    for (int i = start_index; i <= g_state.num_points - (k - level); ++i) {
        g_state.current_simplex_indices[level] = i;
        
        // Check if the new point is connected to all previous points in the current path
        int is_clique = 1;
        for (int j = 0; j < level; ++j) {
            if (dist_sq(i, g_state.current_simplex_indices[j]) > g_state.alpha) {
                is_clique = 0;
                break;
            }
        }

        if (is_clique) {
            find_simplices_recursive(k, level + 1, i + 1, count);
        }
    }
}


// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_points max_simplicial_dimension alpha seed\n", argv[0]);
        exit(1);
    }

    g_state.num_points = atoi(argv[1]);
    g_state.max_simplicial_dimension = atoi(argv[2]);
    g_state.alpha = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (g_state.max_simplicial_dimension < 0) {
        fprintf(stderr, "FATAL: max_simplicial_dimension must be non-negative.\n");
        exit(1);
    }
    
    // Allocate memory for points and helper array
    g_state.points = (Point3D*)malloc(g_state.num_points * sizeof(Point3D));
    if (!g_state.points) {
        fprintf(stderr, "FATAL: Memory allocation failed for points.\n");
        exit(1);
    }
    
    // The helper array needs to store indices for the largest possible simplex
    g_state.current_simplex_indices = (int*)malloc((g_state.max_simplicial_dimension + 1) * sizeof(int));
    if(!g_state.current_simplex_indices) {
        fprintf(stderr, "FATAL: Memory allocation failed for simplex indices.\n");
        free(g_state.points);
        exit(1);
    }

    // Generate random 3D points in the [0, 1] cube
    for (int i = 0; i < g_state.num_points; ++i) {
        g_state.points[i].x = (double)mt_rand() / UINT32_MAX;
        g_state.points[i].y = (double)mt_rand() / UINT32_MAX;
        g_state.points[i].z = (double)mt_rand() / UINT32_MAX;
    }

    g_state.total_simplices_count = 0;
}

void run_computation() {
    // 0-simplices are the points themselves
    g_state.total_simplices_count = g_state.num_points;

    // Find simplices for each dimension from 1 up to max_simplicial_dimension
    for (int d = 1; d <= g_state.max_simplicial_dimension; ++d) {
        int k = d + 1; // A d-simplex has k=d+1 vertices
        long long simplex_count_for_dim = 0;
        find_simplices_recursive(k, 0, 0, &simplex_count_for_dim);
        g_state.total_simplices_count += simplex_count_for_dim;
    }
}

void cleanup() {
    free(g_state.points);
    free(g_state.current_simplex_indices);
}

// --- MAIN FUNCTION ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", g_state.total_simplices_count);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
