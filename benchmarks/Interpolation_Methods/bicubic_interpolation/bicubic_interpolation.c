#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// Benchmark parameters
static int grid_width;
static int grid_height;
static int num_queries;

// Data arrays
static float *grid_data = NULL;
static float *query_x = NULL;
static float *query_y = NULL;

// Final result
static double final_result = 0.0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_width grid_height num_queries seed\n", argv[0]);
        exit(1);
    }

    grid_width = atoi(argv[1]);
    grid_height = atoi(argv[2]);
    num_queries = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    grid_data = (float *)malloc((size_t)grid_width * grid_height * sizeof(float));
    query_x = (float *)malloc((size_t)num_queries * sizeof(float));
    query_y = (float *)malloc((size_t)num_queries * sizeof(float));

    if (!grid_data || !query_x || !query_y) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Populate grid with random values
    for (int i = 0; i < grid_width * grid_height; ++i) {
        grid_data[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    // Populate query points
    // Queries are generated in a sub-region to avoid edge effects in interpolation
    float w_range = (float)grid_width - 4.0f;
    float h_range = (float)grid_height - 4.0f;
    for (int i = 0; i < num_queries; ++i) {
        query_x[i] = 1.0f + ((float)mt_rand() / (float)UINT32_MAX) * w_range;
        query_y[i] = 1.0f + ((float)mt_rand() / (float)UINT32_MAX) * h_range;
    }
}

// 1D cubic interpolation (Catmull-Rom)
static inline float cubic_interpolate(float p[4], float x) {
    return p[1] + 0.5f * x * (p[2] - p[0] + x * (2.0f * p[0] - 5.0f * p[1] + 4.0f * p[2] - p[3] + x * (3.0f * (p[1] - p[2]) + p[3] - p[0])));
}

void run_computation() {
    double total_sum = 0.0;
    float p[4][4];

    for (int i = 0; i < num_queries; ++i) {
        float qx = query_x[i];
        float qy = query_y[i];

        int x_int = (int)qx;
        int y_int = (int)qy;

        float x_frac = qx - x_int;
        float y_frac = qy - y_int;

        // Gather the 4x4 grid of points surrounding the query point
        for (int j = 0; j < 4; ++j) {
            for (int k = 0; k < 4; ++k) {
                p[j][k] = grid_data[(y_int + j - 1) * grid_width + (x_int + k - 1)];
            }
        }

        // Interpolate along x for each of the 4 rows
        float col[4];
        col[0] = cubic_interpolate(p[0], x_frac);
        col[1] = cubic_interpolate(p[1], x_frac);
        col[2] = cubic_interpolate(p[2], x_frac);
        col[3] = cubic_interpolate(p[3], x_frac);

        // Interpolate along y using the results from the x-interpolations
        total_sum += cubic_interpolate(col, y_frac);
    }
    final_result = total_sum;
}

void cleanup() {
    free(grid_data);
    free(query_x);
    free(query_y);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
