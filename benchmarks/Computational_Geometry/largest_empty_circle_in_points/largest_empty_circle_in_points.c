#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of MT19937 ---

// --- Benchmark Globals ---
typedef struct {
    double x, y;
} Point;

int num_points;
Point *points;
double final_result_radius; // To store the final result

// --- Benchmark Functions ---

// Calculates the squared Euclidean distance between two points
static inline double dist_sq(Point p1, Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return dx * dx + dy * dy;
}

// Calculates the circumcenter of a triangle defined by three points.
// Returns false if the points are collinear, true otherwise.
static bool get_circumcenter(Point p1, Point p2, Point p3, Point* center) {
    double D = 2.0 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y));

    if (fabs(D) < 1e-12) {
        return false; // Points are collinear
    }

    double p1_sq = p1.x * p1.x + p1.y * p1.y;
    double p2_sq = p2.x * p2.x + p2.y * p2.y;
    double p3_sq = p3.x * p3.x + p3.y * p3.y;

    center->x = (p1_sq * (p2.y - p3.y) + p2_sq * (p3.y - p1.y) + p3_sq * (p1.y - p2.y)) / D;
    center->y = (p1_sq * (p3.x - p2.x) + p2_sq * (p1.x - p3.x) + p3_sq * (p2.x - p1.x)) / D;

    return true;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_points> <seed>\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    mt_seed(seed);

    points = (Point*)malloc(num_points * sizeof(Point));
    if (!points) {
        fprintf(stderr, "Failed to allocate memory for points.\n");
        exit(1);
    }

    for (int i = 0; i < num_points; ++i) {
        points[i].x = (double)mt_rand() / (double)UINT32_MAX;
        points[i].y = (double)mt_rand() / (double)UINT32_MAX;
    }
    
    final_result_radius = 0.0;
}

void run_computation() {
    double max_radius_sq = 0.0;
    const double epsilon = 1e-9; // For floating point comparisons

    // This benchmark uses a brute-force approach to find the largest empty circle.
    // A circle is a candidate if its boundary is defined by 2 or 3 of the input points.
    // For each candidate circle, we check if it is empty (contains no other points).
    // Total complexity is dominated by the triplet check, resulting in O(N^4).

    // Case 1: Circle defined by 2 points (as diameter)
    // Complexity: O(N^3)
    for (int i = 0; i < num_points; ++i) {
        for (int j = i + 1; j < num_points; ++j) {
            Point center = {(points[i].x + points[j].x) / 2.0, (points[i].y + points[j].y) / 2.0};
            double radius_sq = dist_sq(points[i], center);

            // Optimization: if this can't be a new max, skip the expensive check
            if (radius_sq < max_radius_sq) continue;

            bool is_empty = true;
            for (int k = 0; k < num_points; ++k) {
                if (dist_sq(center, points[k]) < radius_sq - epsilon) {
                    is_empty = false;
                    break;
                }
            }
            if (is_empty) {
                max_radius_sq = radius_sq;
            }
        }
    }

    // Case 2: Circle defined by 3 points (circumcircle)
    // Complexity: O(N^4)
    Point center;
    for (int i = 0; i < num_points; ++i) {
        for (int j = i + 1; j < num_points; ++j) {
            for (int k = j + 1; k < num_points; ++k) {
                if (get_circumcenter(points[i], points[j], points[k], &center)) {
                    double radius_sq = dist_sq(center, points[i]);

                    if (radius_sq < max_radius_sq) continue;
                    
                    // Heuristic: check if center is within the initial bounding box [0,1]x[0,1]
                    if (center.x < 0.0 || center.x > 1.0 || center.y < 0.0 || center.y > 1.0) {
                        continue;
                    }

                    bool is_empty = true;
                    for (int l = 0; l < num_points; ++l) {
                        if (dist_sq(center, points[l]) < radius_sq - epsilon) {
                            is_empty = false;
                            break;
                        }
                    }
                    if (is_empty) {
                        max_radius_sq = radius_sq;
                    }
                }
            }
        }
    }

    final_result_radius = sqrt(max_radius_sq);
}

void cleanup() {
    free(points);
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
    printf("%f\n", final_result_radius);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
