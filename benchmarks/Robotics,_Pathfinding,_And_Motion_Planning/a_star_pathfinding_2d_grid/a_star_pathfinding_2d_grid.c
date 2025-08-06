#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

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

// --- Global Data Structures ---

// Benchmark parameters
static int grid_width;
static int grid_height;

// Input data
static unsigned char *grid; // 0 for walkable, 1 for obstacle
static int start_node_idx;
static int end_node_idx;

// Output result
static volatile int final_path_length;

// A* Node structure
typedef struct {
    int g; // Cost from start
    int h; // Heuristic cost to end
    int parent_idx;
} NodeData;

// Min-heap for the open set
typedef struct {
    int node_idx;
    int f; // Total cost (g + h)
} HeapNode;

// --- Utility Functions ---

static inline int manhattan_distance(int x1, int y1, int x2, int y2) {
    return abs(x1 - x2) + abs(y1 - y2);
}

// --- Min-Heap for Priority Queue ---

static void sift_down(HeapNode* heap, int size, int index) {
    int smallest = index;
    int left = 2 * index + 1;
    int right = 2 * index + 2;

    if (left < size && heap[left].f < heap[smallest].f) {
        smallest = left;
    }
    if (right < size && heap[right].f < heap[smallest].f) {
        smallest = right;
    }

    if (smallest != index) {
        HeapNode temp = heap[index];
        heap[index] = heap[smallest];
        heap[smallest] = temp;
        sift_down(heap, size, smallest);
    }
}

static void sift_up(HeapNode* heap, int index) {
    while (index != 0 && heap[(index - 1) / 2].f > heap[index].f) {
        HeapNode temp = heap[index];
        heap[index] = heap[(index - 1) / 2];
        heap[(index - 1) / 2] = temp;
        index = (index - 1) / 2;
    }
}

static void heap_push(HeapNode** heap, int* size, int* capacity, HeapNode new_node) {
    if (*size == *capacity) {
        *capacity *= 2;
        *heap = (HeapNode*)realloc(*heap, *capacity * sizeof(HeapNode));
    }
    (*heap)[*size] = new_node;
    sift_up(*heap, *size);
    (*size)++;
}

static HeapNode heap_pop(HeapNode* heap, int* size) {
    HeapNode root = heap[0];
    heap[0] = heap[*size - 1];
    (*size)--;
    sift_down(heap, *size, 0);
    return root;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_width grid_height obstacle_percentage seed\n", argv[0]);
        exit(1);
    }

    grid_width = atoi(argv[1]);
    grid_height = atoi(argv[2]);
    double obstacle_percentage = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    grid = (unsigned char *)malloc(grid_width * grid_height * sizeof(unsigned char));
    if (!grid) {
        fprintf(stderr, "Failed to allocate memory for grid.\n");
        exit(1);
    }

    double obs_threshold = obstacle_percentage * UINT32_MAX;
    for (int i = 0; i < grid_width * grid_height; ++i) {
        grid[i] = (mt_rand() < obs_threshold) ? 1 : 0;
    }

    // Select valid start and end points
    do {
        start_node_idx = mt_rand() % (grid_width * grid_height);
    } while (grid[start_node_idx] == 1);

    do {
        end_node_idx = mt_rand() % (grid_width * grid_height);
    } while (grid[end_node_idx] == 1 || end_node_idx == start_node_idx);

    final_path_length = -1; // Default to -1 (no path found)
}

void run_computation() {
    int grid_size = grid_width * grid_height;
    NodeData *nodes = (NodeData *)malloc(grid_size * sizeof(NodeData));
    char *closed_set = (char *)calloc(grid_size, sizeof(char));

    int heap_capacity = 1024;
    int heap_size = 0;
    HeapNode* open_set = (HeapNode*)malloc(heap_capacity*sizeof(HeapNode));

    for (int i = 0; i < grid_size; ++i) {
        nodes[i].g = -1; // Using -1 as infinity
    }
 
    int end_x = end_node_idx % grid_width;
    int end_y = end_node_idx / grid_width;

    nodes[start_node_idx].g = 0;
    nodes[start_node_idx].h = manhattan_distance(start_node_idx % grid_width, start_node_idx / grid_width, end_x, end_y);
    nodes[start_node_idx].parent_idx = -1;
    HeapNode start_heap_node = {start_node_idx, nodes[start_node_idx].h};
    heap_push(&open_set, &heap_size, &heap_capacity, start_heap_node);

    int dx[] = {0, 0, 1, -1}; // 4-directional movement
    int dy[] = {1, -1, 0, 0};

    while (heap_size > 0) {
        HeapNode current_heap_node = heap_pop(open_set, &heap_size);
        int current_idx = current_heap_node.node_idx;

        if (current_idx == end_node_idx) {
            int path_len = 0;
            int curr = end_node_idx;
            while (curr != -1) {
                path_len++;
                curr = nodes[curr].parent_idx;
            }
            final_path_length = path_len - 1; // Number of steps
            goto cleanup_run;
        }

        if (closed_set[current_idx]) continue;
        closed_set[current_idx] = 1;

        int current_x = current_idx % grid_width;
        int current_y = current_idx / grid_width;

        for (int i = 0; i < 4; i++) {
            int neighbor_x = current_x + dx[i];
            int neighbor_y = current_y + dy[i];

            if (neighbor_x >= 0 && neighbor_x < grid_width && neighbor_y >= 0 && neighbor_y < grid_height) {
                int neighbor_idx = neighbor_y * grid_width + neighbor_x;
                if (grid[neighbor_idx] == 1 || closed_set[neighbor_idx]) {
                    continue;
                }

                int tentative_g = nodes[current_idx].g + 1;
                if (nodes[neighbor_idx].g == -1 || tentative_g < nodes[neighbor_idx].g) {
                    nodes[neighbor_idx].parent_idx = current_idx;
                    nodes[neighbor_idx].g = tentative_g;
                    nodes[neighbor_idx].h = manhattan_distance(neighbor_x, neighbor_y, end_x, end_y);
                    int f = tentative_g + nodes[neighbor_idx].h;

                    HeapNode neighbor_heap_node = {neighbor_idx, f};
                    heap_push(&open_set, &heap_size, &heap_capacity, neighbor_heap_node);
                }
            }
        }
    }

cleanup_run:
    free(nodes);
    free(closed_set);
    free(open_set);
}

void cleanup() {
    free(grid);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", final_path_length);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
