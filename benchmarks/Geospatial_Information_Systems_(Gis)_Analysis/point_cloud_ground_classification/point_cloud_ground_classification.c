#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <float.h>

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

// --- BENCHMARK DATA AND PARAMETERS ---
typedef struct {
    float x, y, z;
} Point3D;

// Global pointers for benchmark data
Point3D* g_points = NULL;
int* g_classification = NULL; // 0 for ground, 1 for non-ground
float* g_min_z_grid = NULL;

// Parameters from command line
int P_NUM_POINTS;
int P_GRID_RESOLUTION;
float P_HEIGHT_THRESHOLD;

// Global result accumulator to prevent dead-code elimination
int g_ground_point_count = 0;

// --- UTILITY FUNCTIONS ---
float rand_float(float min, float max) {
    return min + ((float)mt_rand() / (float)UINT32_MAX) * (max - min);
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_points grid_resolution height_threshold seed\n", argv[0]);
        exit(1);
    }

    P_NUM_POINTS = atoi(argv[1]);
    P_GRID_RESOLUTION = atoi(argv[2]);
    P_HEIGHT_THRESHOLD = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    g_points = (Point3D*)malloc(P_NUM_POINTS * sizeof(Point3D));
    if (!g_points) { fprintf(stderr, "Error: Could not allocate memory for points.\n"); exit(1); }
    
    g_classification = (int*)malloc(P_NUM_POINTS * sizeof(int));
    if (!g_classification) { fprintf(stderr, "Error: Could not allocate memory for classification.\n"); exit(1); }

    size_t grid_size = (size_t)P_GRID_RESOLUTION * P_GRID_RESOLUTION;
    g_min_z_grid = (float*)malloc(grid_size * sizeof(float));
    if (!g_min_z_grid) { fprintf(stderr, "Error: Could not allocate memory for Z grid.\n"); exit(1); }

    // Generate random point cloud data
    float max_coord = (float)P_GRID_RESOLUTION;
    for (int i = 0; i < P_NUM_POINTS; i++) {
        g_points[i].x = rand_float(0.0f, max_coord - 0.001f);
        g_points[i].y = rand_float(0.0f, max_coord - 0.001f);
        // Create a sloped plane with random noise for Z to simulate terrain
        float slope_component = 0.05f * g_points[i].x + 0.02f * g_points[i].y;
        g_points[i].z = slope_component + rand_float(0.0f, 10.0f); 
    }

    // Initialize the minimum Z grid to a very large value
    for (size_t i = 0; i < grid_size; i++) {
        g_min_z_grid[i] = FLT_MAX;
    }
}

void run_computation() {
    size_t grid_dim = (size_t)P_GRID_RESOLUTION;

    // Phase 1: Populate the minimum Z grid by iterating through all points
    for (int i = 0; i < P_NUM_POINTS; i++) {
        Point3D p = g_points[i];
        int grid_x = (int)p.x;
        int grid_y = (int)p.y;
        size_t index = (size_t)grid_y * grid_dim + grid_x;

        if (p.z < g_min_z_grid[index]) {
            g_min_z_grid[index] = p.z;
        }
    }

    // Phase 2: Classify each point based on its height relative to the grid minimum
    int ground_count = 0;
    for (int i = 0; i < P_NUM_POINTS; i++) {
        Point3D p = g_points[i];
        int grid_x = (int)p.x;
        int grid_y = (int)p.y;
        size_t index = (size_t)grid_y * grid_dim + grid_x;

        float min_z_in_cell = g_min_z_grid[index];

        if (p.z - min_z_in_cell <= P_HEIGHT_THRESHOLD) {
            g_classification[i] = 0; // Ground
            ground_count++;
        } else {
            g_classification[i] = 1; // Non-ground
        }
    }
    g_ground_point_count = ground_count;
}

void cleanup() {
    free(g_points);
    free(g_classification);
    free(g_min_z_grid);
    g_points = NULL;
    g_classification = NULL;
    g_min_z_grid = NULL;
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
    printf("%d\n", g_ground_point_count);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
