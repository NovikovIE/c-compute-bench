#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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

// --- Benchmark Globals ---
int num_points;
int num_eval_points;

// Control points (y-values). The x-coordinates are implicit (0, 1, 2, ...).
float *control_points_y;

// 't' values at which to evaluate the spline.
float *eval_times;

// Array to store the results of the interpolation.
float *results;

// Final accumulated result to prevent dead-code elimination.
double final_result;

// --- Benchmark Functions ---

// Utility to generate a random float between 0.0 and 1.0
float rand_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <num_eval_points> <seed>\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    num_eval_points = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (num_points < 4) {
        fprintf(stderr, "Error: num_points must be at least 4 for Catmull-Rom splines.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    control_points_y = (float *)malloc(num_points * sizeof(float));
    eval_times = (float *)malloc(num_eval_points * sizeof(float));
    results = (float*)malloc(num_eval_points * sizeof(float));

    if (!control_points_y || !eval_times || !results) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate random control points
    for (int i = 0; i < num_points; ++i) {
        control_points_y[i] = rand_float() * 100.0f; // y-values between 0 and 100
    }

    // Generate random evaluation times within the valid range [0, num_points - 1)
    float max_t = (float)(num_points - 1);
    for (int i = 0; i < num_eval_points; ++i) {
        eval_times[i] = rand_float() * max_t;
    }
}

void run_computation() {
    double accumulator = 0.0;

    for (int i = 0; i < num_eval_points; ++i) {
        float global_t = eval_times[i];

        // The segment is between point P1 (at p1_idx) and P2 (at p1_idx + 1)
        int p1_idx = (int)global_t;

        // The parameter 't' for the formula, from 0 to 1
        float t = global_t - (float)p1_idx;

        // Clamp p1_idx to ensure p2_idx is always in bounds.
        // If p1_idx is at the last point, we treat it as the end of the previous segment.
        if (p1_idx >= num_points - 1) {
            p1_idx = num_points - 2;
            t = 1.0f;
        }

        // Determine indices for the 4 required control points (P0, P1, P2, P3)
        // We clamp at the boundaries to create phantom points.
        int p0_idx = (p1_idx > 0) ? p1_idx - 1 : p1_idx;
        int p2_idx = p1_idx + 1;
        int p3_idx = (p1_idx < num_points - 2) ? p1_idx + 2 : p2_idx;

        // Fetch the y-values of the control points
        float y0 = control_points_y[p0_idx];
        float y1 = control_points_y[p1_idx];
        float y2 = control_points_y[p2_idx];
        float y3 = control_points_y[p3_idx];

        float t2 = t * t;
        float t3 = t2 * t;

        // Catmull-Rom spline formula
        float interpolated_y = 0.5f * (
            (2.0f * y1) +
            (-y0 + y2) * t +
            (2.0f * y0 - 5.0f * y1 + 4.0f * y2 - y3) * t2 +
            (-y0 + 3.0f * y1 - 3.0f * y2 + y3) * t3
        );

        results[i] = interpolated_y;
    }

    // Accumulate results to prevent dead-code elimination by the compiler
    for (int i = 0; i < num_eval_points; ++i) {
        accumulator += results[i];
    }
    final_result = accumulator;
}

void cleanup() {
    free(control_points_y);
    free(eval_times);
    free(results);
    control_points_y = NULL;
    eval_times = NULL;
    results = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
