#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// To compile: gcc -o inverse_distance_weighting_interpolation inverse_distance_weighting_interpolation.c -lm -O3

// --- Mersenne Twister (Do Not Modify) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Data Structures ---
typedef struct {
    float x;
    float y;
    float value;
} KnownPoint;

// Global structure to hold all benchmark data
struct {
    int num_known_points;
    int output_grid_width;
    int output_grid_height;
    float power_parameter;
    
    KnownPoint *known_points;
    float *output_grid;

    long long final_result;
} g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_known_points output_grid_width output_grid_height power_parameter seed\n", argv[0]);
        exit(1);
    }
    
    g_data.num_known_points = atoi(argv[1]);
    g_data.output_grid_width = atoi(argv[2]);
    g_data.output_grid_height = atoi(argv[3]);
    g_data.power_parameter = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);
    
    g_data.known_points = (KnownPoint *)malloc(g_data.num_known_points * sizeof(KnownPoint));
    if (!g_data.known_points) {
        perror("Failed to allocate memory for known_points");
        exit(1);
    }

    size_t grid_size = (size_t)g_data.output_grid_width * g_data.output_grid_height;
    g_data.output_grid = (float *)malloc(grid_size * sizeof(float));
    if (!g_data.output_grid) {
        perror("Failed to allocate memory for output_grid");
        free(g_data.known_points);
        exit(1);
    }

    // Initialize known points with random data in a 1.0 x 1.0 space
    for (int i = 0; i < g_data.num_known_points; ++i) {
        g_data.known_points[i].x = (float)mt_rand() / (float)UINT32_MAX;
        g_data.known_points[i].y = (float)mt_rand() / (float)UINT32_MAX;
        g_data.known_points[i].value = (float)mt_rand() / (float)UINT32_MAX * 100.0f;
    }
}

void run_computation() {
    const float epsilon = 1e-6f;
    double checksum = 0.0;
    
    const float inv_width_minus_1 = (g_data.output_grid_width > 1) ? 1.0f / (g_data.output_grid_width - 1) : 0.0f;
    const float inv_height_minus_1 = (g_data.output_grid_height > 1) ? 1.0f / (g_data.output_grid_height - 1) : 0.0f;

    for (int gy = 0; gy < g_data.output_grid_height; ++gy) {
        float py = (float)gy * inv_height_minus_1;
        for (int gx = 0; gx < g_data.output_grid_width; ++gx) {
            float px = (float)gx * inv_width_minus_1;
            
            float numerator = 0.0f;
            float denominator = 0.0f;
            float interpolated_value = 0.0f;
            int found_exact = 0;

            for (int i = 0; i < g_data.num_known_points; ++i) {
                float dx = px - g_data.known_points[i].x;
                float dy = py - g_data.known_points[i].y;
                float distance_sq = dx * dx + dy * dy;

                if (distance_sq < epsilon * epsilon) {
                    interpolated_value = g_data.known_points[i].value;
                    found_exact = 1;
                    break;
                }
                
                float distance = sqrtf(distance_sq);
                float weight = 1.0f / powf(distance, g_data.power_parameter);
                numerator += weight * g_data.known_points[i].value;
                denominator += weight;
            }
            
            if (!found_exact) {
                if (denominator > epsilon) {
                    interpolated_value = numerator / denominator;
                } else {
                    interpolated_value = 0.0f; // Default value if all points are too far/weights are zero
                }
            }
            
            int grid_index = gy * g_data.output_grid_width + gx;
            g_data.output_grid[grid_index] = interpolated_value;
            checksum += interpolated_value;
        }
    }
    
    g_data.final_result = (long long)checksum;
}

void cleanup() {
    free(g_data.known_points);
    free(g_data.output_grid);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%lld\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
