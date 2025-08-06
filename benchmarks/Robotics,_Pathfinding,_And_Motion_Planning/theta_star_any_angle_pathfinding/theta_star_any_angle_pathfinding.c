#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

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

// Benchmark parameters
int GRID_WIDTH;
int GRID_HEIGHT;
int OBSTACLE_PERCENTAGE;

// Node structure for pathfinding
typedef struct Node {
    int x, y;
    double g_score;       // Cost from start to current node
    double h_score;       // Heuristic cost from a node to the goal
    double f_score;       // g_score + h_score
    struct Node* parent;
    bool in_open_set;
    bool in_closed_set;
} Node;

// Global data structures
Node* grid_nodes;     // 1D array representing the 2D grid of nodes
bool* obstacle_map;   // 1D array for obstacle data
Node** open_set;      // Array of pointers to nodes in the open set
int open_set_count = 0;

Node* start_node = NULL;
Node* goal_node = NULL;

long long final_result = 0; // Accumulated result to prevent dead code elimination

// Helper functions
static inline int get_index(int x, int y) {
    return y * GRID_WIDTH + x;
}

static inline double euclidean_distance(Node* a, Node* b) {
    double dx = a->x - b->x;
    double dy = a->y - b->y;
    return sqrt(dx * dx + dy * dy);
}

bool line_of_sight(Node* a, Node* b) {
    int x0 = a->x, y0 = a->y;
    int x1 = b->x, y1 = b->y;
    int dx = abs(x1 - x0);
    int dy = -abs(y1 - y0);
    int sx = (x0 < x1) ? 1 : -1;
    int sy = (y0 < y1) ? 1 : -1;
    int err = dx + dy;

    while (true) {
        if (obstacle_map[get_index(x0, y0)]) {
            return false;
        }
        if (x0 == x1 && y0 == y1) {
            break;
        }
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
        }
    }
    return true;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_width grid_height obstacle_percentage seed\n", argv[0]);
        exit(1);
    }

    GRID_WIDTH = atoi(argv[1]);
    GRID_HEIGHT = atoi(argv[2]);
    OBSTACLE_PERCENTAGE = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate memory
    grid_nodes = (Node*)malloc(GRID_WIDTH * GRID_HEIGHT * sizeof(Node));
    obstacle_map = (bool*)malloc(GRID_WIDTH * GRID_HEIGHT * sizeof(bool));
    open_set = (Node**)malloc(GRID_WIDTH * GRID_HEIGHT * sizeof(Node*));

    if (!grid_nodes || !obstacle_map || !open_set) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    start_node = &grid_nodes[get_index(1, 1)];
    goal_node = &grid_nodes[get_index(GRID_WIDTH - 2, GRID_HEIGHT - 2)];

    // Initialize grid and obstacles
    for (int y = 0; y < GRID_HEIGHT; ++y) {
        for (int x = 0; x < GRID_WIDTH; ++x) {
            int idx = get_index(x, y);
            grid_nodes[idx] = (Node){
                .x = x, .y = y,
                .g_score = DBL_MAX, .h_score = 0.0,
                .f_score = DBL_MAX, .parent = NULL,
                .in_open_set = false, .in_closed_set = false
            };
            obstacle_map[idx] = (mt_rand() % 100) < OBSTACLE_PERCENTAGE;
        }
    }
    
    // Ensure borders are obstacles
    for (int y = 0; y < GRID_HEIGHT; ++y) {
        obstacle_map[get_index(0, y)] = true;
        obstacle_map[get_index(GRID_WIDTH - 1, y)] = true;
    }
    for (int x = 0; x < GRID_WIDTH; ++x) {
        obstacle_map[get_index(x, 0)] = true;
        obstacle_map[get_index(x, GRID_HEIGHT - 1)] = true;
    }

    // Ensure start and goal are not obstacles
    obstacle_map[get_index(start_node->x, start_node->y)] = false;
    obstacle_map[get_index(goal_node->x, goal_node->y)] = false;

    // Initialize all h_scores (heuristic to goal)
    for (int i = 0; i < GRID_WIDTH * GRID_HEIGHT; ++i) {
        grid_nodes[i].h_score = euclidean_distance(&grid_nodes[i], goal_node);
    }
}

void run_computation() {
    start_node->g_score = 0.0;
    start_node->f_score = start_node->h_score;
    
    open_set[open_set_count++] = start_node;
    start_node->in_open_set = true;

    while (open_set_count > 0) {
        // Find node with the lowest f_score in the open set
        int current_idx_in_set = 0;
        for (int i = 1; i < open_set_count; ++i) {
            if (open_set[i]->f_score < open_set[current_idx_in_set]->f_score) {
                current_idx_in_set = i;
            }
        }
        Node* current = open_set[current_idx_in_set];

        if (current == goal_node) {
            Node* path_node = goal_node;
            while (path_node != NULL) {
                final_result += path_node->x + path_node->y;
                path_node = path_node->parent;
            }
            return;
        }

        // Remove current from open set
        open_set[current_idx_in_set] = open_set[--open_set_count];
        current->in_open_set = false;
        current->in_closed_set = true;

        // Process neighbors
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0) continue;

                int nx = current->x + dx;
                int ny = current->y + dy;

                Node* neighbor = &grid_nodes[get_index(nx, ny)];
                if (neighbor->in_closed_set || obstacle_map[get_index(nx, ny)]) {
                    continue;
                }

                Node* parent = current;
                // Theta* modification: check line of sight from parent's parent
                if (current->parent != NULL && line_of_sight(current->parent, neighbor)) {
                    parent = current->parent;
                }
                
                double new_g_score = parent->g_score + euclidean_distance(parent, neighbor);

                if (new_g_score < neighbor->g_score) {
                    neighbor->parent = parent;
                    neighbor->g_score = new_g_score;
                    neighbor->f_score = new_g_score + neighbor->h_score;

                    if (!neighbor->in_open_set) {
                        open_set[open_set_count++] = neighbor;
                        neighbor->in_open_set = true;
                    }
                }
            }
        }
    }
    // If no path is found, final_result remains 0
}

void cleanup() {
    free(grid_nodes);
    free(obstacle_map);
    free(open_set);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
