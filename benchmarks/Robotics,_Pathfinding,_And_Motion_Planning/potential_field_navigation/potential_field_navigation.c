#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- START: Mersenne Twister (MT19937) --- Do Not Modify ---
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
// --- END: Mersenne Twister ---

// Benchmark data structure
typedef struct {
    int grid_width;
    int grid_height;
    int num_obstacles;
    int max_iterations;

    float* potential_field;
    
    int start_x;
    int start_y;
    
    long long final_result;
} BenchmarkData;

// Global pointer to benchmark data
static BenchmarkData* g_data = NULL;

// Helper to calculate squared distance
static inline float dist_sq(int x1, int y1, int x2, int y2) {
    float dx = (float)(x1 - x2);
    float dy = (float)(y1 - y2);
    return dx * dx + dy * dy;
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <grid_width> <grid_height> <num_obstacles> <max_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        fprintf(stderr, "FATAL: Could not allocate memory for benchmark data.\n");
        exit(1);
    }

    g_data->grid_width = atoi(argv[1]);
    g_data->grid_height = atoi(argv[2]);
    g_data->num_obstacles = atoi(argv[3]);
    g_data->max_iterations = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate grid
    size_t grid_size = (size_t)g_data->grid_width * g_data->grid_height;
    g_data->potential_field = (float*)malloc(grid_size * sizeof(float));
    if (!g_data->potential_field) {
        fprintf(stderr, "FATAL: Could not allocate memory for potential field.\n");
        free(g_data);
        exit(1);
    }
    
    // Define goal and obstacles
    const int goal_x = g_data->grid_width - 10;
    const int goal_y = g_data->grid_height - 10;

    typedef struct { int x, y; } Obstacle;
    Obstacle* obstacles = (Obstacle*)malloc(g_data->num_obstacles * sizeof(Obstacle));
    if (!obstacles) {
        fprintf(stderr, "FATAL: Could not allocate memory for obstacles.\n");
        free(g_data->potential_field);
        free(g_data);
        exit(1);
    }

    for (int i = 0; i < g_data->num_obstacles; ++i) {
        obstacles[i].x = mt_rand() % g_data->grid_width;
        obstacles[i].y = mt_rand() % g_data->grid_height;
    }

    const float attractive_k = 0.005f;
    const float repulsive_influence_sq = 100.0f; // radius of 10

    // Calculate potential for each cell in the grid
    for (int y = 0; y < g_data->grid_height; ++y) {
        for (int x = 0; x < g_data->grid_width; ++x) {
            size_t index = (size_t)y * g_data->grid_width + x;
            
            float attractive_potential = attractive_k * dist_sq(x, y, goal_x, goal_y);
            float repulsive_potential = 0.0f;

            for (int i = 0; i < g_data->num_obstacles; ++i) {
                float d2 = dist_sq(x, y, obstacles[i].x, obstacles[i].y);
                if (d2 < 0.1f) { // on or very close to an obstacle
                    repulsive_potential = INFINITY;
                    break;
                } else if (d2 < repulsive_influence_sq) {
                    repulsive_potential += (1.0f / d2 - 1.0f / repulsive_influence_sq);
                }
            }
            g_data->potential_field[index] = attractive_potential + (50000.0f * repulsive_potential);
        }
    }

    free(obstacles);

    // Set start position and final result accumulator
    g_data->start_x = 10;
    g_data->start_y = 10;
    g_data->final_result = 0;
}

void run_computation() {
    int current_x = g_data->start_x;
    int current_y = g_data->start_y;
    const int width = g_data->grid_width;
    const int height = g_data->grid_height;
    float* field = g_data->potential_field;

    for (int i = 0; i < g_data->max_iterations; ++i) {
        float min_potential = field[(size_t)current_y * width + current_x];
        int next_x = current_x;
        int next_y = current_y;

        // Check 8 neighbors for the steepest descent path
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0) continue;

                int nx = current_x + dx;
                int ny = current_y + dy;

                if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                    size_t index = (size_t)ny * width + nx;
                    if (field[index] < min_potential) {
                        min_potential = field[index];
                        next_x = nx;
                        next_y = ny;
                    }
                }
            }
        }
        current_x = next_x;
        current_y = next_y;
    }

    g_data->final_result = (long long)current_x + current_y;
}

void cleanup() {
    if (g_data) {
        free(g_data->potential_field);
        free(g_data);
        g_data = NULL;
    }
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    long long final_result = g_data->final_result;

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
