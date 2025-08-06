#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// Mersenne Twister (Do Not Modify - Include This Verbatim):
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

// ------------ Benchmark-specific code starts here ------------

// Helper for random floats
double mt_rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Global parameters
static int map_width;
static int map_height;
static int num_sensor_readings;
static int readings_per_scan;

// Data structures
typedef struct {
    float x, y, theta;
} Pose;

typedef struct {
    float distance, angle;
} SensorReading;

// Global pointers for data
static float *occupancy_grid;
static Pose *robot_poses;
static SensorReading *all_sensor_readings;
static int global_result;

// Benchmark constants
#define PI 3.14159265358979323846f
#define MAX_SENSOR_RANGE 50.0f // meters
#define MAP_SCALE 0.1f // meters per cell
#define LOG_ODDS_FREE -0.4f
#define LOG_ODDS_OCC 0.85f
#define LOG_ODDS_CLAMP_MIN -10.0f
#define LOG_ODDS_CLAMP_MAX 10.0f

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s map_width map_height num_sensor_readings readings_per_scan seed\n", argv[0]);
        exit(1);
    }

    map_width = atoi(argv[1]);
    map_height = atoi(argv[2]);
    num_sensor_readings = atoi(argv[3]);
    readings_per_scan = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);
    
    mt_seed(seed);

    // Allocate memory
    occupancy_grid = (float*)malloc((size_t)map_width * map_height * sizeof(float));
    if (!occupancy_grid) { fprintf(stderr, "Memory allocation failed for grid\n"); exit(1); }

    robot_poses = (Pose*)malloc((size_t)num_sensor_readings * sizeof(Pose));
    if (!robot_poses) { fprintf(stderr, "Memory allocation failed for poses\n"); exit(1); }

    const size_t total_readings = (size_t)num_sensor_readings * readings_per_scan;
    all_sensor_readings = (SensorReading*)malloc(total_readings * sizeof(SensorReading));
    if (!all_sensor_readings) { fprintf(stderr, "Memory allocation failed for readings\n"); exit(1); }

    // Initialize map (log-odds = 0 means 50% probability)
    for (int i = 0; i < map_width * map_height; ++i) {
        occupancy_grid[i] = 0.0f;
    }

    // Generate robot poses and sensor readings
    for (int i = 0; i < num_sensor_readings; ++i) {
        // Random pose within the map (with a margin)
        robot_poses[i].x = (float)(mt_rand_double() * (map_width * 0.8) + (map_width * 0.1));
        robot_poses[i].y = (float)(mt_rand_double() * (map_height * 0.8) + (map_height * 0.1));
        robot_poses[i].theta = (float)(mt_rand_double() * 2.0 * PI);
    }

    for (size_t i = 0; i < total_readings; ++i) {
        // Simulate a 180-degree scan
        all_sensor_readings[i].angle = (float)((mt_rand_double() * PI) - (PI / 2.0));
        // Distance in cells
        all_sensor_readings[i].distance = (float)(mt_rand_double() * MAX_SENSOR_RANGE / MAP_SCALE);
    }
}

void run_computation() {
    for (int i = 0; i < num_sensor_readings; ++i) {
        Pose current_pose = robot_poses[i];
        int start_idx_scan = i * readings_per_scan;

        for (int j = 0; j < readings_per_scan; ++j) {
            SensorReading reading = all_sensor_readings[start_idx_scan + j];
            
            float total_angle = current_pose.theta + reading.angle;
            float end_x_f = current_pose.x + reading.distance * cosf(total_angle);
            float end_y_f = current_pose.y + reading.distance * sinf(total_angle);
            
            // Convert to grid coordinates
            int x0 = (int)current_pose.x;
            int y0 = (int)current_pose.y;
            int x1 = (int)end_x_f;
            int y1 = (int)end_y_f;

            // Bresenham's line algorithm to trace the beam
            int dx = abs(x1 - x0);
            int sx = x0 < x1 ? 1 : -1;
            int dy = -abs(y1 - y0);
            int sy = y0 < y1 ? 1 : -1;
            int err = dx + dy;
            int e2;
            
            int cur_x = x0;
            int cur_y = y0;

            while (1) {
                // Mark cells along the ray as free, up to the endpoint
                if (cur_x >= 0 && cur_x < map_width && cur_y >= 0 && cur_y < map_height) {
                    int grid_idx = cur_y * map_width + cur_x;
                    occupancy_grid[grid_idx] += LOG_ODDS_FREE;
                    if (occupancy_grid[grid_idx] < LOG_ODDS_CLAMP_MIN) {
                       occupancy_grid[grid_idx] = LOG_ODDS_CLAMP_MIN;
                    }
                }
                
                if (cur_x == x1 && cur_y == y1) break;
                
                e2 = 2 * err;
                if (e2 >= dy) { err += dy; cur_x += sx; }
                if (e2 <= dx) { err += dx; cur_y += sy; }
            }

            // Update the endpoint as occupied, if within bounds
            if (x1 >= 0 && x1 < map_width && y1 >= 0 && y1 < map_height) {
                int grid_idx = y1 * map_width + x1;
                // Revert the 'free' update and add the 'occupied' update
                occupancy_grid[grid_idx] -= LOG_ODDS_FREE;
                occupancy_grid[grid_idx] += LOG_ODDS_OCC;
                if (occupancy_grid[grid_idx] > LOG_ODDS_CLAMP_MAX) {
                    occupancy_grid[grid_idx] = LOG_ODDS_CLAMP_MAX;
                }
            }
        }
    }

    // Accumulate a result to prevent dead code elimination
    int occupied_cells = 0;
    for (int i = 0; i < map_width * map_height; ++i) {
        if (occupancy_grid[i] > 0.0f) {
            occupied_cells++;
        }
    }
    global_result = occupied_cells;
}

void cleanup() {
    free(occupancy_grid);
    free(robot_poses);
    free(all_sensor_readings);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", global_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
