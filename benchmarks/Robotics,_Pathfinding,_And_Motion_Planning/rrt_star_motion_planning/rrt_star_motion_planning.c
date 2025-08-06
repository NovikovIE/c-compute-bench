#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>

// --- Mersenne Twister (Verbatim) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Data Structures and Globals ---
typedef struct {
    double* coords;
} Point;

typedef struct Node {
    Point state;
    struct Node* parent;
    double cost;
} Node;

typedef struct {
    Point center;
    double radius;
} Obstacle;

// Parameters
static int CONFIG_SPACE_DIMS;
static int NUM_OBSTACLES;
static int NUM_ITERATIONS;
static double REWIRE_RADIUS;

// Data
static Node* tree;
static int tree_size;
static Obstacle* obstacles;
static Point start_point;
static Point goal_point;

// Result
static double final_path_cost;

// Constants
#define STEP_SIZE 0.1
#define GOAL_BIAS 0.05
#define GOAL_RADIUS 0.1

// --- Helper Functions ---
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

double dist_sq(const Point* p1, const Point* p2) {
    double d_sq = 0.0;
    for (int i = 0; i < CONFIG_SPACE_DIMS; i++) {
        double diff = p1->coords[i] - p2->coords[i];
        d_sq += diff * diff;
    }
    return d_sq;
}

// Simple collision check by sampling points along the segment
int is_collision(const Point* p1, const Point* p2) {
    if (CONFIG_SPACE_DIMS > 256) { // Safety check for stack usage
      fprintf(stderr, "FATAL: CONFIG_SPACE_DIMS too large for stack allocation in collision check.\n");
      exit(1);
    }
    double p_check_coords[CONFIG_SPACE_DIMS];
    Point p_check = { .coords = p_check_coords };
    int num_checks = 10;

    for (int i = 1; i <= num_checks; i++) {
        double t = (double)i / num_checks;
        for (int d = 0; d < CONFIG_SPACE_DIMS; d++) {
            p_check.coords[d] = p1->coords[d] + t * (p2->coords[d] - p1->coords[d]);
        }

        for (int j = 0; j < NUM_OBSTACLES; j++) {
            if (dist_sq(&p_check, &obstacles[j].center) <= obstacles[j].radius * obstacles[j].radius) {
                return 1; // Collision
            }
        }
    }
    return 0; // No collision
}


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s config_space_dims num_obstacles num_iterations rewire_radius seed\n", argv[0]);
        exit(1);
    }

    CONFIG_SPACE_DIMS = atoi(argv[1]);
    NUM_OBSTACLES = atoi(argv[2]);
    NUM_ITERATIONS = atoi(argv[3]);
    REWIRE_RADIUS = atof(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    tree = (Node*)malloc((NUM_ITERATIONS + 2) * sizeof(Node));
    if (!tree) { perror("Failed to allocate tree"); exit(1); }

    obstacles = (Obstacle*)malloc(NUM_OBSTACLES * sizeof(Obstacle));
    if (!obstacles) { perror("Failed to allocate obstacles"); exit(1); }

    start_point.coords = (double*)malloc(CONFIG_SPACE_DIMS * sizeof(double));
    goal_point.coords = (double*)malloc(CONFIG_SPACE_DIMS * sizeof(double));
    for (int i = 0; i < CONFIG_SPACE_DIMS; i++) {
        start_point.coords[i] = 0.05;
        goal_point.coords[i] = 0.95;
    }

    tree[0].state.coords = (double*)malloc(CONFIG_SPACE_DIMS * sizeof(double));
    memcpy(tree[0].state.coords, start_point.coords, CONFIG_SPACE_DIMS * sizeof(double));
    tree[0].parent = NULL;
    tree[0].cost = 0.0;
    tree_size = 1;

    for (int i = 0; i < NUM_OBSTACLES; i++) {
        obstacles[i].center.coords = (double*)malloc(CONFIG_SPACE_DIMS * sizeof(double));
        obstacles[i].radius = random_double() * 0.05 + 0.05;
        
        int too_close;
        do {
            too_close = 0;
            for (int d = 0; d < CONFIG_SPACE_DIMS; d++) {
                obstacles[i].center.coords[d] = random_double() * 0.8 + 0.1;
            }
            if (dist_sq(&obstacles[i].center, &start_point) < pow(obstacles[i].radius + 0.1, 2)) too_close = 1;
            if (dist_sq(&obstacles[i].center, &goal_point) < pow(obstacles[i].radius + 0.1, 2)) too_close = 1;
        } while (too_close);
    }

    final_path_cost = -1.0;
}


void run_computation() {
    double rewire_radius_sq = REWIRE_RADIUS * REWIRE_RADIUS;

    for (int k = 0; k < NUM_ITERATIONS && tree_size <= NUM_ITERATIONS; k++) {
        double q_rand_coords[CONFIG_SPACE_DIMS];
        Point q_rand = { .coords = q_rand_coords };

        if (random_double() < GOAL_BIAS) {
            memcpy(q_rand.coords, goal_point.coords, CONFIG_SPACE_DIMS * sizeof(double));
        } else {
            for (int i = 0; i < CONFIG_SPACE_DIMS; i++) {
                q_rand.coords[i] = random_double();
            }
        }

        Node* q_near = &tree[0];
        double min_dist_sq = DBL_MAX;
        for (int i = 0; i < tree_size; i++) {
            double d_sq = dist_sq(&tree[i].state, &q_rand);
            if (d_sq < min_dist_sq) {
                min_dist_sq = d_sq;
                q_near = &tree[i];
            }
        }

        double q_new_coords[CONFIG_SPACE_DIMS];
        Point q_new = { .coords = q_new_coords };
        double dist = sqrt(min_dist_sq);

        if (dist < STEP_SIZE) {
            memcpy(q_new.coords, q_rand.coords, CONFIG_SPACE_DIMS * sizeof(double));
        } else {
            for (int i = 0; i < CONFIG_SPACE_DIMS; i++) {
                q_new.coords[i] = q_near->state.coords[i] + (q_rand.coords[i] - q_near->state.coords[i]) * STEP_SIZE / dist;
            }
        }

        if (!is_collision(&q_near->state, &q_new)) {
            Node* q_parent = q_near;
            double min_cost = q_near->cost + sqrt(dist_sq(&q_near->state, &q_new));
            
            for (int i = 0; i < tree_size; i++) {
                 if (dist_sq(&tree[i].state, &q_new) < rewire_radius_sq) {
                    double cost_via_neighbor = tree[i].cost + sqrt(dist_sq(&tree[i].state, &q_new));
                    if (cost_via_neighbor < min_cost && !is_collision(&tree[i].state, &q_new)) {
                        min_cost = cost_via_neighbor;
                        q_parent = &tree[i];
                    }
                }
            }

            Node* new_node = &tree[tree_size];
            new_node->state.coords = (double*)malloc(CONFIG_SPACE_DIMS * sizeof(double));
            memcpy(new_node->state.coords, q_new.coords, CONFIG_SPACE_DIMS * sizeof(double));
            new_node->parent = q_parent;
            new_node->cost = min_cost;
            
            for (int i = 0; i < tree_size; i++) {
                if (dist_sq(&tree[i].state, &new_node->state) < rewire_radius_sq) {
                    double cost_via_new = new_node->cost + sqrt(dist_sq(&tree[i].state, &new_node->state));
                    if (cost_via_new < tree[i].cost && !is_collision(&tree[i].state, &new_node->state)) {
                        tree[i].parent = new_node;
                        tree[i].cost = cost_via_new;
                    }
                }
            }
            tree_size++;
        }
    }

    Node* goal_node = NULL;
    double min_goal_dist_sq = GOAL_RADIUS * GOAL_RADIUS;
    for (int i = 0; i < tree_size; i++) {
        double d_sq = dist_sq(&tree[i].state, &goal_point);
        if (d_sq < min_goal_dist_sq) {
            final_path_cost = tree[i].cost + sqrt(d_sq);
            min_goal_dist_sq = d_sq;
        }
    }
}


void cleanup() {
    for (int i = 0; i < tree_size; i++) {
        free(tree[i].state.coords);
    }
    free(tree);

    for (int i = 0; i < NUM_OBSTACLES; i++) {
        free(obstacles[i].center.coords);
    }
    free(obstacles);

    free(start_point.coords);
    free(goal_point.coords);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    printf("%f\n", final_path_cost);

    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
