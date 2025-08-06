#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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
// --- End of Mersenne Twister ---

// Benchmark-specific definitions
typedef struct {
    double x, y;
} Point;

// Global variables for data
int num_points;
Point* points = NULL;
Point* hull = NULL;
int result_hull_size = 0;

// Helper function to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// qsort comparison function for points
int compare_points(const void* a, const void* b) {
    Point* p1 = (Point*)a;
    Point* p2 = (Point*)b;
    if (p1->x < p2->x) return -1;
    if (p1->x > p2->x) return 1;
    if (p1->y < p2->y) return -1;
    if (p1->y > p2->y) return 1;
    return 0;
}

// Cross product to find orientation of ordered triplet (p1, p2, p3)
// > 0 -> Counterclockwise (left turn)
// < 0 -> Clockwise (right turn)
// = 0 -> Collinear
double cross_product(Point p1, Point p2, Point p3) {
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_points> <seed>\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if(num_points < 3) {
        fprintf(stderr, "Error: num_points must be at least 3.\n");
        exit(1);
    }

    mt_seed(seed);

    points = (Point*)malloc(num_points * sizeof(Point));
    hull = (Point*)malloc(num_points * sizeof(Point)); // Max possible hull size is num_points

    if (!points || !hull) {
        fprintf(stderr, "Failed to allocate memory\n");
        exit(1);
    }

    // Generate random points in a square
    for (int i = 0; i < num_points; i++) {
        points[i].x = rand_double() * 10000.0;
        points[i].y = rand_double() * 10000.0;
    }
}

void run_computation() {
    // Andrew's monotone chain algorithm to find the convex hull.

    // Step 1: Sort points lexicographically.
    qsort(points, num_points, sizeof(Point), compare_points);

    int k = 0; // Index for hull points array

    // Step 2: Build the lower hull.
    for (int i = 0; i < num_points; ++i) {
        // While the last two points on the hull and the current point make a non-left turn,
        // pop the last point from the hull.
        while (k >= 2 && cross_product(hull[k - 2], hull[k - 1], points[i]) <= 0) {
            k--;
        }
        hull[k++] = points[i];
    }

    // Step 3: Build the upper hull.
    // t stores the size of the lower hull, to not disturb it.
    for (int i = num_points - 2, t = k + 1; i >= 0; i--) {
        // While the last two points on the hull and the current point make a non-left turn,
        // pop the last point from the hull.
        while (k >= t && cross_product(hull[k - 2], hull[k - 1], points[i]) <= 0) {
            k--;
        }
        hull[k++] = points[i];
    }

    // The last point added is a duplicate of the first point, so we exclude it.
    result_hull_size = k - 1;
}

void cleanup() {
    if (points) free(points);
    if (hull) free(hull);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print the number of points in the convex hull to stdout
    printf("%d\n", result_hull_size);

    // Print the computation time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
