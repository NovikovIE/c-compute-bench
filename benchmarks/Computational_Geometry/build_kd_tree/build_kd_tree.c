#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER (MT19937) ---
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

// Benchmark-specific data structures
typedef struct KdNode {
    int point_index;
    int axis;
    struct KdNode *left;
    struct KdNode *right;
} KdNode;

// Global variables to share data between setup, computation, and cleanup
static int NUM_POINTS;
static int DIMENSIONS;

static double** points;       // Array of D-dimensional points
static int* point_indices;   // Indices used for building the tree
static KdNode* root_node = NULL; // Root of the k-d tree
static double final_result = 0.0; // Result to prevent dead code elimination

// Forward declarations of helper functions
static void swap_indices(int* a, int* b);
static int partition(int* indices, int low, int high, int axis);
static void quickselect(int* indices, int low, int high, int k, int axis);
static KdNode* build_kdtree_recursive(int* indices, int count, int depth);
static void free_kdtree(KdNode* node);
static double traverse_and_sum(KdNode* node);

// setup_benchmark: Parses arguments, allocates memory, and generates data.
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <dimensions> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_POINTS = atoi(argv[1]);
    DIMENSIONS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_POINTS <= 0 || DIMENSIONS <= 0) {
        fprintf(stderr, "FATAL: num_points and dimensions must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    points = (double**)malloc(NUM_POINTS * sizeof(double*));
    if (!points) {
        fprintf(stderr, "FATAL: Memory allocation failed for points array.\n");
        exit(1);
    }
    for (int i = 0; i < NUM_POINTS; ++i) {
        points[i] = (double*)malloc(DIMENSIONS * sizeof(double));
        if (!points[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for a point's coordinates.\n");
            exit(1);
        }
    }

    // Generate random points in the range [0.0, 1.0]
    for (int i = 0; i < NUM_POINTS; ++i) {
        for (int j = 0; j < DIMENSIONS; ++j) {
            points[i][j] = (double)mt_rand() / (double)UINT32_MAX;
        }
    }

    // Allocate and initialize point indices for the tree building process
    point_indices = (int*)malloc(NUM_POINTS * sizeof(int));
    if (!point_indices) {
        fprintf(stderr, "FATAL: Memory allocation failed for point indices.\n");
        exit(1);
    }
    for (int i = 0; i < NUM_POINTS; ++i) {
        point_indices[i] = i;
    }
}

// Helper: swap two integer pointers
static void swap_indices(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Helper: Partition for quickselect using Lomuto partition scheme
static int partition(int* indices, int low, int high, int axis) {
    double pivot_val = points[indices[high]][axis];
    int i = low - 1;
    for (int j = low; j < high; ++j) {
        if (points[indices[j]][axis] < pivot_val) {
            i++;
            swap_indices(&indices[i], &indices[j]);
        }
    }
    swap_indices(&indices[i + 1], &indices[high]);
    return i + 1;
}

// Helper: Iterative quickselect to find the k-th smallest element and partition the array
static void quickselect(int* indices, int low, int high, int k, int axis) {
    while (low < high) {
        int pi = partition(indices, low, high, axis);
        if (pi == k) {
            return;
        }
        if (pi > k) {
            high = pi - 1;
        } else {
            low = pi + 1;
        }
    }
}

// Recursive function to build the k-d tree
static KdNode* build_kdtree_recursive(int* indices, int count, int depth) {
    if (count <= 0) {
        return NULL;
    }

    int axis = depth % DIMENSIONS;
    int median_offset = count / 2;

    // Use quickselect to find the median and partition the sub-array of indices
    quickselect(indices, 0, count - 1, median_offset, axis);

    int median_point_idx = indices[median_offset];

    KdNode* node = (KdNode*)malloc(sizeof(KdNode));
    if (!node) {
        fprintf(stderr, "FATAL: Memory allocation failed for KdNode.\n");
        exit(1);
    }
    node->point_index = median_point_idx;
    node->axis = axis;

    // Recursively build left and right subtrees
    node->left = build_kdtree_recursive(indices, median_offset, depth + 1);
    node->right = build_kdtree_recursive(indices + median_offset + 1, count - median_offset - 1, depth + 1);

    return node;
}

// Helper: Traverse the tree and compute a sum to prevent dead code elimination
static double traverse_and_sum(KdNode* node) {
    if (node == NULL) {
        return 0.0;
    }
    // Sum the first coordinate of each point in the tree
    double sum = points[node->point_index][0];
    sum += traverse_and_sum(node->left);
    sum += traverse_and_sum(node->right);
    return sum;
}

// run_computation: Executes the core k-d tree building algorithm
void run_computation() {
    root_node = build_kdtree_recursive(point_indices, NUM_POINTS, 0);

    // To prevent dead code elimination, traverse the tree and calculate a sum.
    if (root_node) {
        final_result = traverse_and_sum(root_node);
    }
}

// Helper: Free k-d tree memory using post-order traversal
static void free_kdtree(KdNode* node) {
    if (node == NULL) {
        return;
    }
    free_kdtree(node->left);
    free_kdtree(node->right);
    free(node);
}

// cleanup: Frees all memory allocated in setup_benchmark
void cleanup() {
    if (root_node) {
        free_kdtree(root_node);
    }

    if (point_indices) {
        free(point_indices);
    }

    if (points) {
        for (int i = 0; i < NUM_POINTS; ++i) {
            if(points[i]) free(points[i]);
        }
        free(points);
    }
    
}

// main: orchestrates the benchmark execution and timing
int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout to be checked for correctness
    // and to prevent dead code elimination.
    printf("%.1f\n", final_result);

    // Print the time taken to stderr.
    fprintf(stderr, "%.6f", time_taken);
    
    // Cleanup must be after printing results.
    cleanup();

    return 0;
}
