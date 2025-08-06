#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// For portability, define M_SQRT2 if not available
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

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
// --- END MERSENNE TWISTER ---

// Benchmark parameters and data
static int LOOKAHEAD_DEPTH;
static int NUM_MOVES_TO_SIMULATE;
static int MAP_SIZE;

static int** map;       // 0 for open, 1 for obstacle
static double** H;      // Heuristic table for LRTA*

static int agent_x, agent_y;
static int goal_x, goal_y;

static unsigned long long final_result;

// Node for A* search, used within run_computation
typedef struct {
    int x, y;
    double g; // cost from start of lookahead
    double f; // g + h
} AStarNode;

// Helper to calculate Manhattan distance
double heuristic_cost_estimate(int x1, int y1, int x2, int y2) {
    return (double)(abs(x1 - x2) + abs(y1 - y2));
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <lookahead_depth> <num_moves_to_simulate> <map_size> <seed>\n", argv[0]);
        exit(1);
    }
    LOOKAHEAD_DEPTH = atoi(argv[1]);
    NUM_MOVES_TO_SIMULATE = atoi(argv[2]);
    MAP_SIZE = atoi(argv[3]);
    unsigned int seed = (unsigned int)strtoul(argv[4], NULL, 10);
    mt_seed(seed);

    map = (int**)malloc(MAP_SIZE * sizeof(int*));
    H = (double**)malloc(MAP_SIZE * sizeof(double*));
    for (int i = 0; i < MAP_SIZE; i++) {
        map[i] = (int*)malloc(MAP_SIZE * sizeof(int));
        H[i] = (double*)malloc(MAP_SIZE * sizeof(double));
    }

    goal_x = mt_rand() % MAP_SIZE;
    goal_y = mt_rand() % MAP_SIZE;

    for (int i = 0; i < MAP_SIZE; i++) {
        for (int j = 0; j < MAP_SIZE; j++) {
            // 20% chance of obstacle
            map[i][j] = (mt_rand() % 100 < 20) ? 1 : 0;
            H[i][j] = heuristic_cost_estimate(j, i, goal_x, goal_y);
        }
    }

    map[goal_y][goal_x] = 0; // Ensure goal is not an obstacle

    do {
        agent_x = mt_rand() % MAP_SIZE;
        agent_y = mt_rand() % MAP_SIZE;
    } while (map[agent_y][agent_x] == 1 || (agent_x == goal_x && agent_y == goal_y));

    final_result = 0;
}

void cleanup() {
    for (int i = 0; i < MAP_SIZE; i++) {
        free(map[i]);
        free(H[i]);
    }
    free(map);
    free(H);
}

void run_computation() {
    size_t open_list_capacity = (LOOKAHEAD_DEPTH + 1) * 8;
    AStarNode* open_list = (AStarNode*)malloc(open_list_capacity * sizeof(AStarNode));
    int* closed_list_map = (int*)calloc((size_t)MAP_SIZE * MAP_SIZE, sizeof(int));
    #define IS_CLOSED(x, y, id) (closed_list_map[(y) * MAP_SIZE + (x)] == id)
    #define SET_CLOSED(x, y, id) (closed_list_map[(y) * MAP_SIZE + (x)] = id)
    int visited_id = 1;

    int dx[] = {0, 0, 1, -1, 1, -1, 1, -1};
    int dy[] = {1, -1, 0, 0, 1, 1, -1, -1};
    double move_cost[] = {1.0, 1.0, 1.0, 1.0, M_SQRT2, M_SQRT2, M_SQRT2, M_SQRT2};

    for (int move = 0; move < NUM_MOVES_TO_SIMULATE; ++move) {
        if (agent_x == goal_x && agent_y == goal_y) {
           // Goal reached, continue for consistent benchmark workload
        }

        double best_successor_f_cost = -1.0;
        int best_next_x = -1, best_next_y = -1;

        // --- One-step lookahead to find best immediate move ---
        for(int i = 0; i < 8; ++i) {
            int nx = agent_x + dx[i];
            int ny = agent_y + dy[i];
            if (nx >= 0 && nx < MAP_SIZE && ny >= 0 && ny < MAP_SIZE && map[ny][nx] == 0) {
                double successor_f_cost = move_cost[i] + H[ny][nx];
                if (best_next_x == -1 || successor_f_cost < best_successor_f_cost) {
                    best_successor_f_cost = successor_f_cost;
                    best_next_x = nx;
                    best_next_y = ny;
                }
            }
        }

        // --- Deeper search for computational workload ---
        int open_list_size = 0;
        int open_list_head = 0;
        SET_CLOSED(agent_x, agent_y, visited_id);
        open_list[open_list_size++] = (AStarNode){agent_x, agent_y, 0.0, H[agent_y][agent_x]};
        int expansions = 0;
        while (open_list_head < open_list_size && expansions < LOOKAHEAD_DEPTH) {
            int best_idx = -1;
            for (int i = open_list_head; i < open_list_size; ++i) {
                if (best_idx == -1 || open_list[i].f < open_list[best_idx].f) {
                    best_idx = i;
                }
            }
            AStarNode current_node = open_list[best_idx];
            open_list[best_idx] = open_list[open_list_head];
            open_list_head++;
            expansions++;
            for (int i = 0; i < 8; ++i) {
                int nx = current_node.x + dx[i];
                int ny = current_node.y + dy[i];
                if (nx >= 0 && nx < MAP_SIZE && ny >= 0 && ny < MAP_SIZE && map[ny][nx] == 0 && !IS_CLOSED(nx, ny, visited_id)) {
                    double new_g = current_node.g + move_cost[i];
                    SET_CLOSED(nx, ny, visited_id);
                    if (open_list_size < open_list_capacity) {
                        open_list[open_list_size++] = (AStarNode){nx, ny, new_g, new_g + H[ny][nx]};
                    }
                }
            }
        }

        if (best_next_x != -1) {
            H[agent_y][agent_x] = best_successor_f_cost;
            agent_x = best_next_x;
            agent_y = best_next_y;
        }

        final_result += agent_x + agent_y;
        visited_id++; 
        if (visited_id == 0) {
             for(size_t i = 0; i < (size_t)MAP_SIZE * MAP_SIZE; ++i) closed_list_map[i] = 0;
             visited_id = 1;
        }
    }

    free(open_list);
    free(closed_list_map);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%llu\n", final_result);

    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
