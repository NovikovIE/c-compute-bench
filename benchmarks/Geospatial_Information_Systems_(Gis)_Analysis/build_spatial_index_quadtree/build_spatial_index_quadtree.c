#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) ---
// Do not modify this section
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
// --- end Mersenne Twister ---

// --- Benchmark Data Structures ---
#define NODE_CAPACITY 4

typedef struct {
    double x, y;
} Point;

// Represents a rectangular region by its top-left corner and dimensions.
typedef struct {
    double x, y, width, height;
} BoundingBox;

typedef struct QuadTreeNode {
    BoundingBox boundary;
    Point* points[NODE_CAPACITY];
    int point_count;
    int depth;
    struct QuadTreeNode* children[4]; // NW, NE, SW, SE
} QuadTreeNode;

// --- Global Benchmark State ---
static int g_num_points;
static int g_max_tree_depth;
static Point* g_points;
static QuadTreeNode* g_root;
static long long g_total_nodes_created; // The final result for stdout

// --- Forward Declarations ---
void quadtree_insert(QuadTreeNode* node, Point* p);
void quadtree_subdivide(QuadTreeNode* node);
QuadTreeNode* create_quadtree_node(BoundingBox boundary, int depth);

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <max_tree_depth> <seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    g_num_points = atoi(argv[1]);
    g_max_tree_depth = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_num_points <= 0 || g_max_tree_depth <= 0) {
        fprintf(stderr, "Error: num_points and max_tree_depth must be positive.\n");
        exit(EXIT_FAILURE);
    }

    mt_seed(seed);

    g_points = (Point*)malloc(g_num_points * sizeof(Point));
    if (!g_points) {
        perror("Failed to allocate points array");
        exit(EXIT_FAILURE);
    }

    // Generate points in a unit square [0, 1] x [0, 1]
    const double max_val = (double)UINT32_MAX;
    for (int i = 0; i < g_num_points; ++i) {
        g_points[i].x = (double)mt_rand() / max_val;
        g_points[i].y = (double)mt_rand() / max_val;
    }
}

QuadTreeNode* create_quadtree_node(BoundingBox boundary, int depth) {
    QuadTreeNode* node = (QuadTreeNode*)malloc(sizeof(QuadTreeNode));
    if (!node) {
        perror("Failed to allocate quadtree node");
        exit(EXIT_FAILURE);
    }
    node->boundary = boundary;
    node->point_count = 0;
    node->depth = depth;
    memset(node->children, 0, sizeof(node->children));
    g_total_nodes_created++;
    return node;
}

int get_quadrant(BoundingBox* boundary, Point* p) {
    double mid_x = boundary->x + boundary->width / 2.0;
    double mid_y = boundary->y + boundary->height / 2.0;

    if (p->y < mid_y) {
        return (p->x < mid_x) ? 0 : 1; // Top half: 0=NW, 1=NE
    } else {
        return (p->x < mid_x) ? 2 : 3; // Bottom half: 2=SW, 3=SE
    }
}

void quadtree_subdivide(QuadTreeNode* node) {
    double hw = node->boundary.width / 2.0;
    double hh = node->boundary.height / 2.0;
    double x = node->boundary.x;
    double y = node->boundary.y;
    int next_depth = node->depth + 1;

    node->children[0] = create_quadtree_node((BoundingBox){x, y, hw, hh}, next_depth);       // NW
    node->children[1] = create_quadtree_node((BoundingBox){x + hw, y, hw, hh}, next_depth); // NE
    node->children[2] = create_quadtree_node((BoundingBox){x, y + hh, hw, hh}, next_depth); // SW
    node->children[3] = create_quadtree_node((BoundingBox){x + hw, y + hh, hw, hh}, next_depth); // SE
}

void quadtree_insert(QuadTreeNode* node, Point* p) {
    // If it's an internal node, recurse into the correct child
    if (node->children[0] != NULL) {
        int quadrant = get_quadrant(&node->boundary, p);
        quadtree_insert(node->children[quadrant], p);
        return;
    }

    // It's a leaf node. Add the point if there is space.
    if (node->point_count < NODE_CAPACITY) {
        node->points[node->point_count++] = p;
        return;
    }

    // Leaf is full. Subdivide if max depth not reached.
    if (node->depth < g_max_tree_depth) {
        quadtree_subdivide(node);

        // Re-insert existing points into the new children
        for (int i = 0; i < node->point_count; ++i) {
            int quadrant = get_quadrant(&node->boundary, node->points[i]);
            quadtree_insert(node->children[quadrant], node->points[i]);
        }
        node->point_count = 0; // Parent no longer stores points directly

        // Insert the new point into the correct new child
        int new_quadrant = get_quadrant(&node->boundary, p);
        quadtree_insert(node->children[new_quadrant], p);
    }
    // If leaf is full and at max depth, we cannot insert the point. It is dropped.
}

void run_computation() {
    g_total_nodes_created = 0;
    BoundingBox root_boundary = {0.0, 0.0, 1.0, 1.0};
    g_root = create_quadtree_node(root_boundary, 0);

    for (int i = 0; i < g_num_points; ++i) {
        quadtree_insert(g_root, &g_points[i]);
    }
}

void free_quadtree(QuadTreeNode* node) {
    if (node == NULL) {
        return;
    }
    if (node->children[0] != NULL) {
        for (int i = 0; i < 4; ++i) {
            free_quadtree(node->children[i]);
        }
    }
    free(node);
}

void cleanup() {
    if (g_points) {
        free(g_points);
    }
    if (g_root) {
        free_quadtree(g_root);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start_time, end_time;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start_time);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end_time);

    cleanup();

    double elapsed_time = (end_time.tv_sec - start_time.tv_sec) + 
                          (end_time.tv_nsec - start_time.tv_nsec) / 1e9;

    // Print result (total nodes) to stdout
    printf("%lld\n", g_total_nodes_created);

    // Print time to stderr
    fprintf(stderr, "%.6f", elapsed_time);

    return 0;
}
