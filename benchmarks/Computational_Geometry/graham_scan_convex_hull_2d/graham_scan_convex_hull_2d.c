#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND GLOBALS ---
typedef struct {
    double x, y;
} Point;

int num_points;
Point *points;
Point p0; // The pivot point with the lowest y-coordinate

Point *hull; // Array to store the resulting convex hull
int hull_size; // Number of points in the hull
double final_result; // Accumulated result to prevent dead code elimination


// --- GRAHAM SCAN HELPER FUNCTIONS ---
// Swap two points
void swap(Point *a, Point *b) {
    Point temp = *a;
    *a = *b;
    *b = temp;
}

// Returns the square of the Euclidean distance between two points
double distSq(Point p1, Point p2) {
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

// Finds the orientation of the ordered triplet (p, q, r).
// Returns: 0 if collinear, 1 if clockwise, 2 if counter-clockwise
int orientation(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (fabs(val) < 1e-9) return 0; // Use tolerance for floating point collinearity
    return (val > 0) ? 2 : 1;      // Counter-clockwise or Clockwise
}

// Comparison function for qsort. Sorts points based on polar angle with p0.
int compare(const void *vp1, const void *vp2) {
    Point *p1 = (Point *)vp1;
    Point *p2 = (Point *)vp2;

    int o = orientation(p0, *p1, *p2);
    if (o == 0) { // Collinear points
        // The point closer to p0 comes first
        return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;
    }
    // Sort counter-clockwise
    return (o == 2) ? -1 : 1;
}

// --- BENCHMARK CORE FUNCTIONS ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_points> <seed>\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_points <= 0) {
        fprintf(stderr, "FATAL: num_points must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    points = (Point *)malloc(num_points * sizeof(Point));
    hull = (Point *)malloc(num_points * sizeof(Point));
    if (!points || !hull) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_points; i++) {
        points[i].x = (double)mt_rand() / (double)UINT32_MAX * 10000.0;
        points[i].y = (double)mt_rand() / (double)UINT32_MAX * 10000.0;
    }
}

void run_computation() {
    if (num_points < 3) {
        // The hull is just the points themselves
        for (int i = 0; i < num_points; i++) {
            hull[i] = points[i];
        }
        hull_size = num_points;
    } else {
        // 1. Find the bottom-most point (the pivot, p0)
        int min_idx = 0;
        for (int i = 1; i < num_points; i++) {
            if (points[i].y < points[min_idx].y || (points[i].y == points[min_idx].y && points[i].x < points[min_idx].x)) {
                min_idx = i;
            }
        }

        // Place the pivot at the beginning of the array
        swap(&points[0], &points[min_idx]);
        p0 = points[0];

        // 2. Sort the n-1 remaining points based on polar angle with p0
        qsort(&points[1], num_points - 1, sizeof(Point), compare);

        // 3. Create the hull using a stack-like approach
        // The 'hull' array will act as the stack
        hull[0] = points[0];
        hull[1] = points[1];
        hull[2] = points[2];
        hull_size = 3;

        for (int i = 3; i < num_points; i++) {
            // Keep popping from the stack while the turn is not counter-clockwise (left)
            while (hull_size > 1 && orientation(hull[hull_size - 2], hull[hull_size - 1], points[i]) != 2) {
                hull_size--;
            }
            hull[hull_size++] = points[i];
        }
    }
    
    // Calculate a final result to prevent optimization
    final_result = 0.0;
    for (int i = 0; i < hull_size; i++) {
        final_result += hull[i].x + hull[i].y;
    }
}

void cleanup() {
    free(points);
    free(hull);
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
    printf("%.2f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
