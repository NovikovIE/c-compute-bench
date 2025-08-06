#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

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
// --- End of MT19937 --- 

// --- Benchmark Data Structures and Globals ---

// Maximum and minimum number of children in a node
#define MAX_CHILDREN 8
#define MIN_CHILDREN (MAX_CHILDREN / 2)

typedef struct {
    float min_x, min_y, max_x, max_y;
} Rectangle;

typedef struct RTreeNode {
    int is_leaf;
    int num_children;
    Rectangle mbr;
    struct RTreeNode* parent;
    // Leaf nodes point to Rectangles, internal nodes point to other RTreeNodes
    void* children[MAX_CHILDREN];
} RTreeNode;

// Parameters
static int NUM_RECTANGLES;

// Input data
static Rectangle* rectangles;

// R-Tree root
static RTreeNode* root = NULL;

// Final result to prevent dead code elimination
static long long total_nodes_created = 0;

// --- R-Tree Helper Functions ---

float rect_area(Rectangle r) {
    if (r.min_x > r.max_x || r.min_y > r.max_y) return 0.0f;
    return (r.max_x - r.min_x) * (r.max_y - r.min_y);
}

Rectangle combine_rects(Rectangle r1, Rectangle r2) {
    Rectangle r = {
        .min_x = fminf(r1.min_x, r2.min_x),
        .min_y = fminf(r1.min_y, r2.min_y),
        .max_x = fmaxf(r1.max_x, r2.max_x),
        .max_y = fmaxf(r1.max_y, r2.max_y)
    };
    return r;
}

RTreeNode* create_node(RTreeNode* parent, int is_leaf) {
    RTreeNode* node = (RTreeNode*)malloc(sizeof(RTreeNode));
    if (!node) {
        fprintf(stderr, "FATAL: Memory allocation failed for RTreeNode.\n");
        exit(1);
    }
    node->is_leaf = is_leaf;
    node->num_children = 0;
    node->parent = parent;
    node->mbr = (Rectangle){FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX};
    total_nodes_created++;
    return node;
}

void update_mbr(RTreeNode* node) {
    if (node->num_children == 0) {
        node->mbr = (Rectangle){FLT_MAX, FLT_MAX, -FLT_MAX, -FLT_MAX};
        return;
    }

    if (node->is_leaf) {
        node->mbr = *((Rectangle*)node->children[0]);
        for (int i = 1; i < node->num_children; ++i) {
            node->mbr = combine_rects(node->mbr, *((Rectangle*)node->children[i]));
        }
    } else {
        node->mbr = ((RTreeNode*)node->children[0])->mbr;
        for (int i = 1; i < node->num_children; ++i) {
            node->mbr = combine_rects(node->mbr, ((RTreeNode*)node->children[i])->mbr);
        }
    }
}

RTreeNode* choose_leaf(RTreeNode* node, Rectangle* rect) {
    if (node->is_leaf) {
        return node;
    }

    float min_enlargement = FLT_MAX;
    RTreeNode* best_child = NULL;

    for (int i = 0; i < node->num_children; ++i) {
        RTreeNode* child = (RTreeNode*)node->children[i];
        float old_area = rect_area(child->mbr);
        Rectangle enlarged_mbr = combine_rects(child->mbr, *rect);
        float enlargement = rect_area(enlarged_mbr) - old_area;

        if (enlargement < min_enlargement) {
            min_enlargement = enlargement;
            best_child = child;
        } else if (enlargement == min_enlargement) {
            if (rect_area(child->mbr) < rect_area(best_child->mbr)) {
                best_child = child;
            }
        }
    }
    return choose_leaf(best_child, rect);
}

void quadratic_split(RTreeNode* node, void* new_entry, RTreeNode** new_node_ptr);
void adjust_tree(RTreeNode* node, RTreeNode* split_node);

void add_child_to_node(RTreeNode* node, void* child) {
    if (node->is_leaf) {
        node->children[node->num_children++] = child;
    } else { // internal node
        RTreeNode* child_node = (RTreeNode*)child;
        node->children[node->num_children++] = child_node;
        child_node->parent = node;
    }
    update_mbr(node);
}

void insert(Rectangle* rect) {
    RTreeNode* leaf = choose_leaf(root, rect);

    if (leaf->num_children < MAX_CHILDREN) {
        add_child_to_node(leaf, rect);
        adjust_tree(leaf, NULL);
    } else { // Leaf needs to be split
        RTreeNode* new_leaf = NULL;
        quadratic_split(leaf, rect, &new_leaf);
        adjust_tree(leaf, new_leaf);
    }
}

void adjust_tree(RTreeNode* node, RTreeNode* split_node) {
    while (node != NULL) {
        RTreeNode* parent = node->parent;
        update_mbr(node);
        if (split_node != NULL) {
            if (parent != NULL) {
                if (parent->num_children < MAX_CHILDREN) {
                    add_child_to_node(parent, split_node);
                    split_node = NULL; // Split is absorbed
                } else { // Parent must also split
                    RTreeNode* new_parent_split = NULL;
                    quadratic_split(parent, split_node, &new_parent_split);
                    split_node = new_parent_split;
                }
            } else { // Node is the root
                RTreeNode* new_root = create_node(NULL, 0);
                add_child_to_node(new_root, node);
                add_child_to_node(new_root, split_node);
                root = new_root;
                split_node = NULL; // Tree has a new root
            }
        }
        node = parent;
    }
}


void quadratic_split(RTreeNode* node, void* new_entry, RTreeNode** new_node_ptr) {
    void* all_children[MAX_CHILDREN + 1];
    for(int i=0; i < MAX_CHILDREN; ++i) all_children[i] = node->children[i];
    all_children[MAX_CHILDREN] = new_entry;

    int seed1 = 0, seed2 = 1;
    float max_waste = -1.0f;

    // PickSeeds: O(n^2) search for the worst pair
    for (int i = 0; i <= MAX_CHILDREN; ++i) {
        for (int j = i + 1; j <= MAX_CHILDREN; ++j) {
            Rectangle r1 = node->is_leaf ? *((Rectangle*)all_children[i]) : ((RTreeNode*)all_children[i])->mbr;
            Rectangle r2 = node->is_leaf ? *((Rectangle*)all_children[j]) : ((RTreeNode*)all_children[j])->mbr;
            float waste = rect_area(combine_rects(r1, r2)) - rect_area(r1) - rect_area(r2);
            if (waste > max_waste) {
                max_waste = waste;
                seed1 = i;
                seed2 = j;
            }
        }
    }

    RTreeNode* new_node = create_node(node->parent, node->is_leaf);
    *new_node_ptr = new_node;

    node->num_children = 0;
    
    add_child_to_node(node, all_children[seed1]);
    add_child_to_node(new_node, all_children[seed2]);

    int used[MAX_CHILDREN + 1] = {0};
    used[seed1] = 1;
    used[seed2] = 1;

    // Distribute remaining children
    for (int i = 0; i <= MAX_CHILDREN; ++i) {
        if (used[i]) continue;

        if (node->num_children + (MAX_CHILDREN + 1 - (node->num_children + new_node->num_children)) <= MIN_CHILDREN){
            add_child_to_node(node, all_children[i]);
            continue;
        }
        if (new_node->num_children + (MAX_CHILDREN + 1 - (node->num_children + new_node->num_children)) <= MIN_CHILDREN){
            add_child_to_node(new_node, all_children[i]);
            continue;
        }

        Rectangle r = node->is_leaf ? *((Rectangle*)all_children[i]) : ((RTreeNode*)all_children[i])->mbr;
        
        float enlargement1 = rect_area(combine_rects(node->mbr, r)) - rect_area(node->mbr);
        float enlargement2 = rect_area(combine_rects(new_node->mbr, r)) - rect_area(new_node->mbr);

        if (enlargement1 < enlargement2) {
            add_child_to_node(node, all_children[i]);
        } else if (enlargement2 < enlargement1) {
            add_child_to_node(new_node, all_children[i]);
        } else { // Tie break: smaller area, then fewer children
            if (rect_area(node->mbr) < rect_area(new_node->mbr)) {
                 add_child_to_node(node, all_children[i]);
            } else if (rect_area(new_node->mbr) < rect_area(node->mbr)){
                add_child_to_node(new_node, all_children[i]);
            } else {
                if(node->num_children <= new_node->num_children)
                    add_child_to_node(node, all_children[i]);
                else
                    add_child_to_node(new_node, all_children[i]);
            }
        }
    }
}

// --- Benchmark-specific Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_rectangles> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_RECTANGLES = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);

    if (NUM_RECTANGLES <= 0) {
        fprintf(stderr, "FATAL: num_rectangles must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    rectangles = (Rectangle*)malloc(NUM_RECTANGLES * sizeof(Rectangle));
    if (rectangles == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for rectangles.\n");
        exit(1);
    }

    const float MAX_COORD = 10000.0f;
    const float MAX_SIZE = 100.0f;

    for (int i = 0; i < NUM_RECTANGLES; i++) {
        float x1 = (mt_rand() / (float)UINT32_MAX) * MAX_COORD;
        float y1 = (mt_rand() / (float)UINT32_MAX) * MAX_COORD;
        float w = (mt_rand() / (float)UINT32_MAX) * MAX_SIZE + 1.0f;
        float h = (mt_rand() / (float)UINT32_MAX) * MAX_SIZE + 1.0f;
        rectangles[i] = (Rectangle){.min_x = x1, .min_y = y1, .max_x = x1 + w, .max_y = y1 + h};
    }

    root = create_node(NULL, 1); // Create root node, is a leaf initially
}

void run_computation() {
    for (int i = 0; i < NUM_RECTANGLES; ++i) {
        insert(&rectangles[i]);
    }
}

void free_rtree_nodes(RTreeNode* node) {
    if (node == NULL) return;
    if (!node->is_leaf) {
        for (int i = 0; i < node->num_children; ++i) {
            free_rtree_nodes((RTreeNode*)node->children[i]);
        }
    }
    free(node);
}

void cleanup() {
    if (root != NULL) {
        free_rtree_nodes(root);
    }
    free(rectangles);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%lld\n", total_nodes_created);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
