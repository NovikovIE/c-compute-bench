/*
 * Program: prm_construction_phase
 * Theme: Robotics, Pathfinding, and Motion Planning
 * Description: Simulates the construction phase of a Probabilistic Roadmap (PRM).
 *              This involves sampling random points (nodes) in a high-dimensional
 *              configuration space and checking if they are collision-free with
 *              respect to a set of obstacles.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- START MERSENNE TWISTER (MT19937) --- Do Not Modify ---
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

// --- BENCHMARK DATA STRUCTURES ---
typedef struct {
    float* center;    // Center point in N-dimensional space
    float radius_sq;  // Squared radius for efficient collision checking
} Obstacle;

// --- GLOBAL STATE ---
static int g_config_space_dims;
static int g_num_obstacles;
static int g_num_nodes_to_sample;

static Obstacle* g_obstacles;
static float* g_sample_point_buffer; // Reusable buffer for sampled points

static int g_final_result; // The computed result to prevent dead code elimination

// --- BENCHMARK FUNCTIONS ---

// Allocates memory and generates input data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s config_space_dims num_obstacles num_nodes_to_sample seed\n", argv[0]);
        exit(1);
    }

    g_config_space_dims = atoi(argv[1]);
    g_num_obstacles = atoi(argv[2]);
    g_num_nodes_to_sample = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Allocate obstacles
    g_obstacles = (Obstacle*)malloc(g_num_obstacles * sizeof(Obstacle));
    if (!g_obstacles) {
        fprintf(stderr, "FATAL: Memory allocation failed for obstacles.\n");
        exit(1);
    }

    // Generate random obstacles (hyper-spheres)
    for (int i = 0; i < g_num_obstacles; i++) {
        g_obstacles[i].center = (float*)malloc(g_config_space_dims * sizeof(float));
        if (!g_obstacles[i].center) {
            fprintf(stderr, "FATAL: Memory allocation failed for obstacle center.\n");
            exit(1);
        }
        for (int d = 0; d < g_config_space_dims; d++) {
            // Each coordinate is a float between 0.0 and 1.0
            g_obstacles[i].center[d] = (float)mt_rand() / (float)UINT32_MAX;
        }
        // Radius is a small value, so space isn't completely blocked
        float radius = ((float)mt_rand() / (float)UINT32_MAX) * 0.05f + 0.01f;
        g_obstacles[i].radius_sq = radius * radius; // Store squared radius
    }

    // Allocate a reusable buffer for sampling points
    g_sample_point_buffer = (float*)malloc(g_config_space_dims * sizeof(float));
    if (!g_sample_point_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed for sample point buffer.\n");
        exit(1);
    }

    g_final_result = 0;
}

// Runs the core computation
void run_computation() {
    int valid_nodes_count = 0;
    for (int i = 0; i < g_num_nodes_to_sample; i++) {
        // 1. Sample a random point in the configuration space
        for (int d = 0; d < g_config_space_dims; d++) {
            g_sample_point_buffer[d] = (float)mt_rand() / (float)UINT32_MAX;
        }

        // 2. Check for collision with all obstacles
        int is_valid = 1; // Assume valid until a collision is found
        for (int j = 0; j < g_num_obstacles; j++) {
            float dist_sq = 0.0f;
            for (int d = 0; d < g_config_space_dims; d++) {
                float diff = g_sample_point_buffer[d] - g_obstacles[j].center[d];
                dist_sq += diff * diff;
            }

            if (dist_sq < g_obstacles[j].radius_sq) {
                is_valid = 0; // Collision detected
                break;        // No need to check other obstacles
            }
        }

        if (is_valid) {
            valid_nodes_count++;
        }
    }
    g_final_result = valid_nodes_count;
}

// Frees all allocated memory
void cleanup() {
    if (g_obstacles) {
        for (int i = 0; i < g_num_obstacles; i++) {
            free(g_obstacles[i].center);
        }
        free(g_obstacles);
    }
    free(g_sample_point_buffer);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (number of valid nodes) to stdout
    printf("%d\n", g_final_result);

    // Print the timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
