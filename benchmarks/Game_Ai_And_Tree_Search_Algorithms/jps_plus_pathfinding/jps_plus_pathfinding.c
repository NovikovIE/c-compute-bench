#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- Benchmark Data Structures and Globals ---

typedef struct {
    int x, y;
} Point;

typedef struct Node {
    Point pos;
    double g, f;
    struct Node *parent;
} Node;

typedef struct {
    Node **data;
    int size;
    int capacity;
} PriorityQueue;

static struct {
    int grid_width;
    int grid_height;
    int obstacle_percentage;
    int num_paths;
    bool *grid;
    Point *starts;
    Point *ends;
    long long final_result;

    // JPS Search-specific data reused for each run
    Node *node_pool;
    bool *closed_set;
    PriorityQueue open_list;
} g_state;

// --- Forward Declarations for JPS ---
long long jps_find_path(Point start_pos, Point end_pos);

// --- Priority Queue (Min-Heap) Implementation ---
void pq_init(PriorityQueue *pq, int capacity) {
    pq->data = (Node **)malloc(capacity * sizeof(Node *));
    pq->size = 0;
    pq->capacity = capacity;
}

void pq_swap(Node **a, Node **b) {
    Node *temp = *a;
    *a = *b;
    *b = temp;
}

void pq_heapify_up(PriorityQueue *pq, int index) {
    if (index == 0) return;
    int parent_index = (index - 1) / 2;
    if (pq->data[index]->f < pq->data[parent_index]->f) {
        pq_swap(&pq->data[index], &pq->data[parent_index]);
        pq_heapify_up(pq, parent_index);
    }
}

void pq_push(PriorityQueue *pq, Node *node) {
    if (pq->size >= pq->capacity) return; // Should not happen with proper capacity
    pq->data[pq->size] = node;
    pq_heapify_up(pq, pq->size);
    pq->size++;
}

void pq_heapify_down(PriorityQueue *pq, int index) {
    int left = 2 * index + 1;
    int right = 2 * index + 2;
    int smallest = index;
    if (left < pq->size && pq->data[left]->f < pq->data[smallest]->f) {
        smallest = left;
    }
    if (right < pq->size && pq->data[right]->f < pq->data[smallest]->f) {
        smallest = right;
    }
    if (smallest != index) {
        pq_swap(&pq->data[index], &pq->data[smallest]);
        pq_heapify_down(pq, smallest);
    }
}

Node* pq_pop(PriorityQueue *pq) {
    if (pq->size == 0) return NULL;
    Node *top = pq->data[0];
    pq->data[0] = pq->data[pq->size - 1];
    pq->size--;
    pq_heapify_down(pq, 0);
    return top;
}

void pq_clear(PriorityQueue *pq) {
    pq->size = 0;
}

void pq_destroy(PriorityQueue *pq) {
    free(pq->data);
}

// --- JPS Algorithm Implementation ---

inline int get_grid_index(Point p) {
    return p.y * g_state.grid_width + p.x;
}

inline bool is_walkable(Point p) {
    return p.x >= 0 && p.x < g_state.grid_width && 
           p.y >= 0 && p.y < g_state.grid_height && 
           !g_state.grid[get_grid_index(p)];
}

inline double heuristic(Point a, Point b) {
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2)); // Euclidean distance
}

Node* jump(Point current, Point parent, const Point end_pos) {
    int dx = current.x - parent.x;
    int dy = current.y - parent.y;

    if (!is_walkable(current)) return NULL;
    if (current.x == end_pos.x && current.y == end_pos.y) return &g_state.node_pool[get_grid_index(current)];

    // Check for forced neighbors
    if (dx != 0 && dy != 0) { // Diagonal move
        if ((is_walkable((Point){current.x - dx, current.y + dy}) && !is_walkable((Point){current.x - dx, current.y})) ||
            (is_walkable((Point){current.x + dx, current.y - dy}) && !is_walkable((Point){current.x, current.y - dy}))) {
            return &g_state.node_pool[get_grid_index(current)];
        }
         // Check for sub-path jumps
        if (jump((Point){current.x + dx, current.y}, current, end_pos) || jump((Point){current.x, current.y + dy}, current, end_pos)) {
            return &g_state.node_pool[get_grid_index(current)];
        }
    } else if (dx != 0) { // Horizontal move
        if ((is_walkable((Point){current.x + dx, current.y + 1}) && !is_walkable((Point){current.x, current.y + 1})) ||
            (is_walkable((Point){current.x + dx, current.y - 1}) && !is_walkable((Point){current.x, current.y - 1}))) {
            return &g_state.node_pool[get_grid_index(current)];
        }
    } else { // Vertical move
        if ((is_walkable((Point){current.x + 1, current.y + dy}) && !is_walkable((Point){current.x + 1, current.y})) ||
            (is_walkable((Point){current.x - 1, current.y + dy}) && !is_walkable((Point){current.x - 1, current.y}))) {
            return &g_state.node_pool[get_grid_index(current)];
        }
    }

    // Continue jumping
    Point next = {current.x + dx, current.y + dy};
    if (dx != 0 && dy != 0) {
         if (is_walkable(next)){
             return jump(next, current, end_pos);
         }
    } else {
         if (is_walkable(current)) {
            return jump(next, current, end_pos);
         }
    }
    
    return NULL;
}

void identify_successors(Node *node, const Point end_pos) {
    int_fast8_t neighbors[8][2];
    int num_neighbors = 0;

    Point p = node->pos;
    Node *parent_node = node->parent;
    
    if (parent_node == NULL) { // Starting node
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                if (i == 0 && j == 0) continue;
                neighbors[num_neighbors][0] = i;
                neighbors[num_neighbors][1] = j;
                num_neighbors++;
            }
        }
    } else {
        Point parent_pos = parent_node->pos;
        int dx = (p.x - parent_pos.x) / fmax(1, abs(p.x - parent_pos.x));
        int dy = (p.y - parent_pos.y) / fmax(1, abs(p.y - parent_pos.y));
        
        if (dx != 0 && dy != 0) { // Diagonal
            if (is_walkable((Point){p.x, p.y + dy})) { neighbors[num_neighbors][0] = 0; neighbors[num_neighbors][1] = dy; num_neighbors++; }
            if (is_walkable((Point){p.x + dx, p.y})) { neighbors[num_neighbors][0] = dx; neighbors[num_neighbors][1] = 0; num_neighbors++; }
            if ((is_walkable((Point){p.x, p.y + dy}) || is_walkable((Point){p.x + dx, p.y})) && is_walkable((Point){p.x + dx, p.y + dy})) {
                neighbors[num_neighbors][0] = dx; neighbors[num_neighbors][1] = dy; num_neighbors++;
            }
        } else { // Cardinal
            if (dx == 0) { // Vertical
                if(is_walkable((Point){p.x, p.y+dy})){ neighbors[num_neighbors][0] = 0; neighbors[num_neighbors][1] = dy; num_neighbors++; }
                if(is_walkable((Point){p.x+1, p.y}) && !is_walkable((Point){p.x+1, p.y-dy})) { neighbors[num_neighbors][0] = 1; neighbors[num_neighbors][1] = 0; num_neighbors++; }
                if(is_walkable((Point){p.x-1, p.y}) && !is_walkable((Point){p.x-1, p.y-dy})) { neighbors[num_neighbors][0] = -1; neighbors[num_neighbors][1] = 0; num_neighbors++;}
            } else { // Horizontal
                if(is_walkable((Point){p.x+dx, p.y})){ neighbors[num_neighbors][0] = dx; neighbors[num_neighbors][1] = 0; num_neighbors++; }
                if(is_walkable((Point){p.x, p.y+1}) && !is_walkable((Point){p.x-dx, p.y+1})) { neighbors[num_neighbors][0] = 0; neighbors[num_neighbors][1] = 1; num_neighbors++; }
                if(is_walkable((Point){p.x, p.y-1}) && !is_walkable((Point){p.x-dx, p.y-1})) { neighbors[num_neighbors][0] = 0; neighbors[num_neighbors][1] = -1; num_neighbors++; }
            }
        }
    }

    for (int i = 0; i < num_neighbors; ++i) {
        Node *jp = jump((Point){p.x + neighbors[i][0], p.y + neighbors[i][1]}, p, end_pos);
        if (jp) {
             int jp_idx = get_grid_index(jp->pos);
             double new_g = node->g + heuristic(node->pos, jp->pos);
             if (!g_state.closed_set[jp_idx] || new_g < jp->g) {
                jp->g = new_g;
                jp->f = new_g + heuristic(jp->pos, end_pos);
                jp->parent = node;
                if (!g_state.closed_set[jp_idx]) {
                     pq_push(&g_state.open_list, jp);
                     g_state.closed_set[jp_idx] = true;
                }
             } 
        }
    }
}

long long reconstruct_path_and_sum(Node *end_node) {
    long long sum = 0;
    Node *current = end_node;
    while (current != NULL) {
        sum += current->pos.x + current->pos.y;
        current = current->parent;
    }
    return sum;
}

long long jps_find_path(Point start_pos, Point end_pos) {
    int grid_size = g_state.grid_width * g_state.grid_height;
    memset(g_state.closed_set, 0, grid_size * sizeof(bool));
    pq_clear(&g_state.open_list);

    int start_idx = get_grid_index(start_pos);
    Node *start_node = &g_state.node_pool[start_idx];
    start_node->pos = start_pos;
    start_node->g = 0;
    start_node->f = heuristic(start_pos, end_pos);
    start_node->parent = NULL;

    pq_push(&g_state.open_list, start_node);
    g_state.closed_set[start_idx] = true;

    while (g_state.open_list.size > 0) {
        Node *current = pq_pop(&g_state.open_list);
        if (current->pos.x == end_pos.x && current->pos.y == end_pos.y) {
            return reconstruct_path_and_sum(current);
        }
        identify_successors(current, end_pos);
    }
    return 0; // No path found
}

// --- Benchmark Setup, Computation, and Cleanup ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_width grid_height obstacle_percentage seed\n", argv[0]);
        exit(1);
    }
    g_state.grid_width = atoi(argv[1]);
    g_state.grid_height = atoi(argv[2]);
    g_state.obstacle_percentage = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    mt_seed(seed);

    g_state.num_paths = 50; // Fixed number of pathfinding runs
    int grid_size = g_state.grid_width * g_state.grid_height;
    
    g_state.grid = (bool *)malloc(grid_size * sizeof(bool));
    g_state.starts = (Point *)malloc(g_state.num_paths * sizeof(Point));
    g_state.ends = (Point *)malloc(g_state.num_paths * sizeof(Point));
    g_state.node_pool = (Node *)malloc(grid_size * sizeof(Node));
    g_state.closed_set = (bool *)malloc(grid_size * sizeof(bool));

    if (!g_state.grid || !g_state.starts || !g_state.ends || !g_state.node_pool || !g_state.closed_set) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize Node Pool
     for (int i = 0; i < grid_size; ++i) {
        g_state.node_pool[i].pos.y = i / g_state.grid_width;
        g_state.node_pool[i].pos.x = i % g_state.grid_width;
        g_state.node_pool[i].g = -1; // -1 indicates not visited yet
    }
    
    // Generate grid
    for (int i = 0; i < grid_size; i++) {
        g_state.grid[i] = (mt_rand() % 100) < g_state.obstacle_percentage;
    }

    // Generate pathfinding tasks
    for (int i = 0; i < g_state.num_paths; i++) {
        do {
            g_state.starts[i].x = mt_rand() % g_state.grid_width;
            g_state.starts[i].y = mt_rand() % g_state.grid_height;
        } while (g_state.grid[get_grid_index(g_state.starts[i])]);

        do {
            g_state.ends[i].x = mt_rand() % g_state.grid_width;
            g_state.ends[i].y = mt_rand() % g_state.grid_height;
        } while (g_state.grid[get_grid_index(g_state.ends[i])] || 
                 (g_state.starts[i].x == g_state.ends[i].x && g_state.starts[i].y == g_state.ends[i].y));
    }

    // Initialize priority queue
    pq_init(&g_state.open_list, grid_size);
}

void run_computation() {
    g_state.final_result = 0;
    for (int i = 0; i < g_state.num_paths; ++i) {
        long long path_sum = jps_find_path(g_state.starts[i], g_state.ends[i]);
        g_state.final_result += path_sum;
    }
}

void cleanup() {
    free(g_state.grid);
    free(g_state.starts);
    free(g_state.ends);
    free(g_state.node_pool);
    free(g_state.closed_set);
    pq_destroy(&g_state.open_list);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", g_state.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
