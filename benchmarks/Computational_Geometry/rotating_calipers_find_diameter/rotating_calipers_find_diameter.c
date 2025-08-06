#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Mersenne Twister (DO NOT MODIFY)
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

// Benchmark-specific data structures and globals
typedef struct {
    double x;
    double y;
} Point;

int N; // num_convex_hull_points
Point *convex_hull;
double max_dist_sq; // Final result: squared diameter

// Helper to sort angles for convex polygon generation
int compare_doubles(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da > db) return 1;
    if (da < db) return -1;
    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_convex_hull_points> <seed>\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    mt_seed(seed);

    if (N <= 1) {
        convex_hull = NULL;
        return;
    }

    convex_hull = (Point *)malloc(N * sizeof(Point));
    if (!convex_hull) {
        fprintf(stderr, "Failed to allocate memory for convex hull points.\n");
        exit(1);
    }

    // Generate points for a convex polygon by creating points on a circle/ellipse
    // with varying radii and sorting them by angle.
    double *angles = (double *)malloc(N * sizeof(double));
    if (!angles) {
        fprintf(stderr, "Failed to allocate memory for angles.\n");
        free(convex_hull);
        exit(1);
    }

    // Generate N random angles
    for (int i = 0; i < N; i++) {
        angles[i] = (mt_rand() / (double)UINT32_MAX) * 2.0 * M_PI;
    }

    // Sort angles to ensure points form a convex shape when connected in order
    qsort(angles, N, sizeof(double), compare_doubles);

    // Create points based on sorted angles and random radii
    for (int i = 0; i < N; i++) {
        double radius = 1000.0 + (mt_rand() / (double)UINT32_MAX) * 100.0;
        convex_hull[i].x = radius * cos(angles[i]);
        convex_hull[i].y = radius * sin(angles[i]);
    }

    free(angles);
}

// Helper to calculate squared distance between two points
double dist_sq(Point *p1, Point *p2) {
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    return dx * dx + dy * dy;
}

// Helper to compute 2 * signed area of triangle(a, b, c)
// This is the 2D cross-product (b-a) x (c-a)
double cross_product_area(Point *a, Point *b, Point *c) {
    return (b->x - a->x) * (c->y - a->y) - (b->y - a->y) * (c->x - a->x);
}

void run_computation() {
    if (N < 2) {
        max_dist_sq = 0.0;
        return;
    }

    max_dist_sq = 0.0;
    int j = 1;

    // The rotating calipers algorithm
    // i and j are indices of the two 'calipers' (antipodal points/edges)
    // They both sweep over the polygon, each making one full turn.
    for (int i = 0; i < N; i++) {
        int next_i = (i + 1) % N;
        
        // Move caliper j forward while the area of triangle(i, next_i, j) increases
        while (cross_product_area(&convex_hull[i], &convex_hull[next_i], &convex_hull[(j + 1) % N]) > 
               cross_product_area(&convex_hull[i], &convex_hull[next_i], &convex_hull[j])) {
            j = (j + 1) % N;
        }

        // The pair of points with the maximal distance must include a vertex.
        // The candidates are the current antipodal pair (i, j) and (next_i, j).
        double d1_sq = dist_sq(&convex_hull[i], &convex_hull[j]);
        if (d1_sq > max_dist_sq) {
            max_dist_sq = d1_sq;
        }

        double d2_sq = dist_sq(&convex_hull[next_i], &convex_hull[j]);
        if (d2_sq > max_dist_sq) {
            max_dist_sq = d2_sq;
        }
    }
}

void cleanup() {
    if (convex_hull) {
        free(convex_hull);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (squared diameter) to stdout
    printf("%f\n", max_dist_sq);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
