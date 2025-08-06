#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator --- DO NOT MODIFY ---
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

// --- Benchmark Data Structures and Globals ---

#define K_DIM 2

typedef struct {
    float x, y;
} Point;

typedef struct {
    int point_index; // Index into the total_points array
    int left, right; // Index into the kd_nodes array. -1 for null.
} KDNode;

// Benchmark parameters
int num_total_points;
int num_query_points;
float query_radius;

// Data arrays
Point* total_points;    // All points to build the tree from
Point* query_points;    // Points to perform range searches on
KDNode* kd_nodes;       // Pre-allocated array for tree nodes
int* point_indices;     // Array of indices to be sorted during build

// Result accumulator
long long final_result = 0;

// Globals for K-D tree build
static int current_axis; // Used by the qsort comparator
static int next_node_idx; // Node allocator index

// --- Helper Functions ---

// Utility to generate a random float between 0.0 and 1.0
float rand_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

// Calculate squared Euclidean distance between two points
float dist_sq(const Point* p1, const Point* p2) {
    float dx = p1->x - p2->x;
    float dy = p1->y - p2->y;
    return dx * dx + dy * dy;
}

// qsort comparator for sorting point indices based on the current axis
int point_comparator(const void* a, const void* b) {
    int idx1 = *(const int*)a;
    int idx2 = *(const int*)b;
    float diff;
    if (current_axis == 0) { // Compare x-coordinates
        diff = total_points[idx1].x - total_points[idx2].x;
    } else { // Compare y-coordinates
        diff = total_points[idx1].y - total_points[idx2].y;
    }
    if (diff < 0.0f) return -1;
    if (diff > 0.0f) return 1;
    return 0;
}

// Recursive function to build the K-D tree
int build_kdtree_recursive(int start, int count, int depth) {
    if (count <= 0) {
        return -1;
    }

    current_axis = depth % K_DIM;
    qsort(point_indices + start, count, sizeof(int), point_comparator);

    int median_offset = (count - 1) / 2;
    int node_idx = next_node_idx++;
    kd_nodes[node_idx].point_index = point_indices[start + median_offset];

    kd_nodes[node_idx].left = build_kdtree_recursive(start, median_offset, depth + 1);
    kd_nodes[node_idx].right = build_kdtree_recursive(start + median_offset + 1, count - median_offset - 1, depth + 1);

    return node_idx;
}

// Recursive function for range search
void search_kdtree_recursive(int node_idx, const Point* query, float radius_sq, int depth, long long* count) {
    if (node_idx == -1) {
        return;
    }

    const Point* current_point = &total_points[kd_nodes[node_idx].point_index];

    if (dist_sq(current_point, query) <= radius_sq) {
        (*count)++;
    }

    int axis = depth % K_DIM;
    float delta; // Distance from query point to the splitting plane
    if (axis == 0) {
        delta = query->x - current_point->x;
    } else {
        delta = query->y - current_point->y;
    }

    // Search the near child first
    int near_child = (delta < 0) ? kd_nodes[node_idx].left : kd_nodes[node_idx].right;
    int far_child = (delta < 0) ? kd_nodes[node_idx].right : kd_nodes[node_idx].left;

    search_kdtree_recursive(near_child, query, radius_sq, depth + 1, count);

    // Pruning step: only search the far child if the query radius crosses the splitting plane
    if (delta * delta <= radius_sq) {
        search_kdtree_recursive(far_child, query, radius_sq, depth + 1, count);
    }
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_total_points num_query_points query_radius seed\n", argv[0]);
        exit(1);
    }

    num_total_points = atoi(argv[1]);
    num_query_points = atoi(argv[2]);
    query_radius = atof(argv[3]);
    uint32_t seed = atoi(argv[4]);
    mt_seed(seed);

    total_points = (Point*)malloc(num_total_points * sizeof(Point));
    query_points = (Point*)malloc(num_query_points * sizeof(Point));
    
    if (!total_points || !query_points) {
        fprintf(stderr, "Memory allocation failed for points.\n");
        exit(1);
    }

    for (int i = 0; i < num_total_points; ++i) {
        total_points[i].x = rand_float();
        total_points[i].y = rand_float();
    }

    for (int i = 0; i < num_query_points; ++i) {
        query_points[i].x = rand_float();
        query_points[i].y = rand_float();
    }
}

void run_computation() {
    // Allocate memory for tree structure and indices inside computation
    kd_nodes = (KDNode*)malloc(num_total_points * sizeof(KDNode));
    point_indices = (int*)malloc(num_total_points * sizeof(int));
    if (!kd_nodes || !point_indices) {
        fprintf(stderr, "Memory allocation failed for tree structures.\n");
        exit(1);
    }
    
    for (int i = 0; i < num_total_points; ++i) {
        point_indices[i] = i;
    }

    // 1. Build the K-D Tree
    next_node_idx = 0;
    int root_idx = build_kdtree_recursive(0, num_total_points, 0);

    // 2. Perform Range Searches
    long long total_found = 0;
    float radius_sq = query_radius * query_radius;

    for (int i = 0; i < num_query_points; ++i) {
        search_kdtree_recursive(root_idx, &query_points[i], radius_sq, 0, &total_found);
    }

    final_result = total_found;
    
    // Free memory allocated within this function
    free(kd_nodes);
    kd_nodes = NULL;
    free(point_indices);
    point_indices = NULL;
}

void cleanup() {
    free(total_points);
    total_points = NULL;
    free(query_points);
    query_points = NULL;
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
