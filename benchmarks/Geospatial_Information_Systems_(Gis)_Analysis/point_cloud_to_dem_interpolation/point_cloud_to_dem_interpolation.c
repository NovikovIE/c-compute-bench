#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Data and Structures ---
typedef struct {
    double x, y, z;
} Point3D;

enum InterpolationMethod {
    NEAREST, // 0
    IDW      // 1
};

struct {
    // Parameters
    int num_points;
    int grid_width;
    int grid_height;
    enum InterpolationMethod method;

    // Data
    Point3D* points;
    double* dem_grid;
    double final_result;
} benchmark_data;

// Helper to generate a random double in a range
double random_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_points grid_width grid_height interpolation_method seed\n", argv[0]);
        fprintf(stderr, "interpolation_method can be 'NEAREST' or 'IDW'\n");
        exit(1);
    }

    benchmark_data.num_points = atoi(argv[1]);
    benchmark_data.grid_width = atoi(argv[2]);
    benchmark_data.grid_height = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if (strcmp(argv[4], "NEAREST") == 0) {
        benchmark_data.method = NEAREST;
    } else if (strcmp(argv[4], "IDW") == 0) {
        benchmark_data.method = IDW;
    } else {
        fprintf(stderr, "FATAL: Unknown interpolation method '%s'\n", argv[4]);
        exit(1);
    }

    mt_seed(seed);

    benchmark_data.points = (Point3D*)malloc(benchmark_data.num_points * sizeof(Point3D));
    if (!benchmark_data.points) {
        fprintf(stderr, "FATAL: Memory allocation failed for points.\n");
        exit(1);
    }

    for (int i = 0; i < benchmark_data.num_points; ++i) {
        benchmark_data.points[i].x = random_double(0.0, benchmark_data.grid_width);
        benchmark_data.points[i].y = random_double(0.0, benchmark_data.grid_height);
        benchmark_data.points[i].z = random_double(0.0, 1000.0); // Elevation
    }

    size_t grid_size = (size_t)benchmark_data.grid_width * benchmark_data.grid_height;
    benchmark_data.dem_grid = (double*)malloc(grid_size * sizeof(double));
    if (!benchmark_data.dem_grid) {
        fprintf(stderr, "FATAL: Memory allocation failed for DEM grid.\n");
        free(benchmark_data.points);
        exit(1);
    }
    
    benchmark_data.final_result = 0.0;
}

void run_computation() {
    double total_elevation_sum = 0.0;
    const int width = benchmark_data.grid_width;
    const int height = benchmark_data.grid_height;
    const int num_pts = benchmark_data.num_points;
    const double epsilon = 1e-9; // To handle grid point being on a data point

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double interpolated_z = 0.0;

            if (benchmark_data.method == NEAREST) {
                double min_dist_sq = DBL_MAX;
                for (int i = 0; i < num_pts; ++i) {
                    double dx = benchmark_data.points[i].x - (double)x;
                    double dy = benchmark_data.points[i].y - (double)y;
                    double dist_sq = dx * dx + dy * dy;
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                        interpolated_z = benchmark_data.points[i].z;
                    }
                }
            } else { // IDW
                double numerator = 0.0;
                double denominator = 0.0;
                int exact_match = 0;

                for (int i = 0; i < num_pts; ++i) {
                    double dx = benchmark_data.points[i].x - (double)x;
                    double dy = benchmark_data.points[i].y - (double)y;
                    double dist_sq = dx * dx + dy * dy;

                    if (dist_sq < epsilon) {
                        interpolated_z = benchmark_data.points[i].z;
                        exact_match = 1;
                        break;
                    }

                    double weight = 1.0 / dist_sq; // IDW with power=2
                    numerator += weight * benchmark_data.points[i].z;
                    denominator += weight;
                }
                if (!exact_match && denominator > 0) {
                    interpolated_z = numerator / denominator;
                }
            }
            benchmark_data.dem_grid[y * width + x] = interpolated_z;
            total_elevation_sum += interpolated_z;
        }
    }
    benchmark_data.final_result = total_elevation_sum;
}

void cleanup() {
    free(benchmark_data.points);
    free(benchmark_data.dem_grid);
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print result to stdout
    printf("%f\n", benchmark_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
