#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
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

// Benchmark Configuration
#define STEP_SIZE 0.1

// Global variables for benchmark data
int CONFIG_SPACE_DIMS;
int NUM_OBSTACLES;
int NUM_ITERATIONS;

typedef struct {
    double *center;
    double radius_sq;
} Obstacle;

typedef struct {
    double *config;
    int parent_index;
} Node;

Obstacle *obstacles;
Node *tree;
int tree_size;
double total_tree_length = 0.0;

// Utility function to generate a random double between 0 and 1
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Utility function to calculate squared Euclidean distance
double dist_sq(const double *p1, const double *p2, int dims) {
    double d_sq = 0.0;
    for (int i = 0; i < dims; ++i) {
        double diff = p1[i] - p2[i];
        d_sq += diff * diff;
    }
    return d_sq;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <config_space_dims> <num_obstacles> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    CONFIG_SPACE_DIMS = atoi(argv[1]);
    NUM_OBSTACLES = atoi(argv[2]);
    NUM_ITERATIONS = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate obstacles
    obstacles = (Obstacle *)malloc(NUM_OBSTACLES * sizeof(Obstacle));
    for (int i = 0; i < NUM_OBSTACLES; ++i) {
        obstacles[i].center = (double *)malloc(CONFIG_SPACE_DIMS * sizeof(double));
        for (int d = 0; d < CONFIG_SPACE_DIMS; ++d) {
            obstacles[i].center[d] = rand_double();
        }
        double radius = 0.02 + rand_double() * 0.08; // Radius between 0.02 and 0.1
        obstacles[i].radius_sq = radius * radius;
    }

    // Allocate the RRT tree structure and node configurations
    tree = (Node *)malloc((NUM_ITERATIONS + 1) * sizeof(Node));
    for (int i = 0; i < NUM_ITERATIONS + 1; ++i) {
        tree[i].config = (double *)malloc(CONFIG_SPACE_DIMS * sizeof(double));
    }

    // Initialize the tree with a root node
    for (int d = 0; d < CONFIG_SPACE_DIMS; ++d) {
        tree[0].config[d] = 0.05; // Starting point
    }
    tree[0].parent_index = -1; // Root has no parent
    tree_size = 1;
}

void run_computation() {
    double q_rand[CONFIG_SPACE_DIMS];
    double q_new[CONFIG_SPACE_DIMS];

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        // 1. Sample a random point in the configuration space
        for (int d = 0; d < CONFIG_SPACE_DIMS; ++d) {
            q_rand[d] = rand_double();
        }

        // 2. Find the nearest node in the tree ('q_near')
        int nearest_node_idx = -1;
        double min_dist_sq = -1.0;
        for (int j = 0; j < tree_size; ++j) {
            double d_sq = dist_sq(tree[j].config, q_rand, CONFIG_SPACE_DIMS);
            if (nearest_node_idx == -1 || d_sq < min_dist_sq) {
                min_dist_sq = d_sq;
                nearest_node_idx = j;
            }
        }
        Node *q_near = &tree[nearest_node_idx];

        // 3. Steer from q_near towards q_rand to get a new point 'q_new'
        double *near_config = q_near->config;
        double direction_len = sqrt(dist_sq(near_config, q_rand, CONFIG_SPACE_DIMS));
        if (direction_len < 1e-9) continue; // Avoid division by zero

        for (int d = 0; d < CONFIG_SPACE_DIMS; ++d) {
            q_new[d] = near_config[d] + (q_rand[d] - near_config[d]) * STEP_SIZE / direction_len;
        }

        // 4. Check for collision of q_new with any obstacle
        int collision = 0;
        for (int j = 0; j < NUM_OBSTACLES; ++j) {
            if (dist_sq(q_new, obstacles[j].center, CONFIG_SPACE_DIMS) <= obstacles[j].radius_sq) {
                collision = 1;
                break;
            }
        }

        if (collision) {
            continue;
        }

        // 5. Add q_new to the tree if no collision was found
        if (tree_size <= NUM_ITERATIONS) {
            for (int d = 0; d < CONFIG_SPACE_DIMS; ++d) {
                tree[tree_size].config[d] = q_new[d];
            }
            tree[tree_size].parent_index = nearest_node_idx;
            total_tree_length += sqrt(dist_sq(q_new, near_config, CONFIG_SPACE_DIMS));
            tree_size++;
        }
    }
}

void cleanup() {
    for (int i = 0; i < NUM_ITERATIONS + 1; ++i) {
        free(tree[i].config);
    }
    free(tree);

    for (int i = 0; i < NUM_OBSTACLES; ++i) {
        free(obstacles[i].center);
    }
    free(obstacles);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%.6f\n", total_tree_length);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}