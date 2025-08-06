/*
 * BENCHMARK: Natural Neighbor Interpolation
 * THEME: Interpolation Methods
 * DESCRIPTION: A computationally-intensive routine that mimics Natural Neighbor
 *              Interpolation. For each query point, it calculates weights based on
 *              the geometric properties of triangles formed with pairs of known data
 *              points. This O(N*M^2) approach avoids complex data structures like
 *              Delaunay triangulations but provides a significant computational load,
 *              making it suitable for a CPU benchmark.
 */

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
static int NUM_POINTS;
static int NUM_QUERIES;

// Data structures
typedef struct {
    double x, y, value;
} Point;

static Point *points = NULL;
static Point *queries = NULL;
static double *interpolated_values = NULL;
static double final_result = 0.0;

// Helper to generate a random double in [0, 1]
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Helper for squared Euclidean distance
double dist_sq(Point p1, Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return dx * dx + dy * dy;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_points num_queries seed\n", argv[0]);
        exit(1);
    }

    NUM_POINTS = atoi(argv[1]);
    NUM_QUERIES = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    points = (Point *)malloc(NUM_POINTS * sizeof(Point));
    queries = (Point *)malloc(NUM_QUERIES * sizeof(Point));
    interpolated_values = (double *)malloc(NUM_QUERIES * sizeof(double));

    if (!points || !queries || !interpolated_values) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate known data points
    for (int i = 0; i < NUM_POINTS; i++) {
        points[i].x = random_double() * 1000.0;
        points[i].y = random_double() * 1000.0;
        points[i].value = random_double() * 100.0;
    }

    // Generate query points
    for (int i = 0; i < NUM_QUERIES; i++) {
        queries[i].x = random_double() * 1000.0;
        queries[i].y = random_double() * 1000.0;
        queries[i].value = 0; // To be computed
    }
}

void run_computation() {
    for (int i = 0; i < NUM_QUERIES; i++) {
        Point q = queries[i];
        double total_weight = 0.0;
        double weighted_sum = 0.0;

        // This is a computationally intensive proxy for NNI.
        // For each known point p_j, we calculate a weight based on its
        // geometric relationship with the query point q and all other known points p_k.
        // The weight for p_j is an aggregation of values derived from triangles (q, p_j, p_k).
        for (int j = 0; j < NUM_POINTS; j++) {
            double current_point_weight = 0.0;
            for (int k = 0; k < NUM_POINTS; k++) {
                if (j == k) continue;

                double d_q_j_sq = dist_sq(q, points[j]);
                double d_q_k_sq = dist_sq(q, points[k]);
                double d_j_k_sq = dist_sq(points[j], points[k]);

                // Using a formulation that is sensitive to the geometry of the triangle (q, pj, pk).
                // A smaller perimeter implies a more 'compact' triangle, which we give a higher weight.
                if (d_q_j_sq > 1e-12 && d_q_k_sq > 1e-12 && d_j_k_sq > 1e-12) {
                    double perimeter = sqrt(d_q_j_sq) + sqrt(d_q_k_sq) + sqrt(d_j_k_sq);
                    current_point_weight += 1.0 / (perimeter * perimeter * perimeter); // O(N^2) weight calculation
                }
            }

            weighted_sum += points[j].value * current_point_weight;
            total_weight += current_point_weight;
        }

        if (total_weight > 1e-12) {
            interpolated_values[i] = weighted_sum / total_weight;
        } else {
            // Fallback to nearest neighbor if weights are negligible (e.g., collinear points)
            double min_dist_sq = 1e30;
            int nn_idx = -1;
            for (int j = 0; j < NUM_POINTS; j++) {
                double d_sq = dist_sq(q, points[j]);
                if (d_sq < min_dist_sq) {
                    min_dist_sq = d_sq;
                    nn_idx = j;
                }
            }
            interpolated_values[i] = (nn_idx != -1) ? points[nn_idx].value : 0.0;
        }
    }

    // Accumulate result to prevent dead code elimination
    double sum = 0.0;
    for (int i = 0; i < NUM_QUERIES; i++) {
        sum += interpolated_values[i];
    }
    final_result = sum;
}

void cleanup() {
    free(points);
    free(queries);
    free(interpolated_values);
    points = NULL;
    queries = NULL;
    interpolated_values = NULL;
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
