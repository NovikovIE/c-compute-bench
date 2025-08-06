#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Data Structures ---
struct Point {
    double x;
    double y;
};

// Global variables to hold benchmark data
int num_polygon_vertices;
int num_test_points;
struct Point *polygon;
struct Point *test_points;
int final_result;

// Helper to generate a random double in [0, 1)
double mt_rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_polygon_vertices> <num_test_points> <seed>\n", argv[0]);
        exit(1);
    }

    num_polygon_vertices = atoi(argv[1]);
    num_test_points = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_polygon_vertices < 3) {
        fprintf(stderr, "Error: A polygon must have at least 3 vertices.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory on the heap
    polygon = (struct Point *)malloc(num_polygon_vertices * sizeof(struct Point));
    test_points = (struct Point *)malloc(num_test_points * sizeof(struct Point));

    if (!polygon || !test_points) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a simple, star-shaped polygon to avoid self-intersections.
    // Points are generated around a circle with some radius variation.
    #ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif
    for (int i = 0; i < num_polygon_vertices; ++i) {
        double angle = 2.0 * M_PI * i / num_polygon_vertices;
        double radius = 1.0 + 0.5 * (mt_rand_double() - 0.5); // Radius variation
        polygon[i].x = radius * cos(angle);
        polygon[i].y = radius * sin(angle);
    }

    // Generate random test points in a square region [-2, 2] x [-2, 2]
    // This region fully encloses the polygon.
    for (int i = 0; i < num_test_points; ++i) {
        test_points[i].x = 4.0 * mt_rand_double() - 2.0;
        test_points[i].y = 4.0 * mt_rand_double() - 2.0;
    }
}

void run_computation() {
    int points_inside = 0;
    for (int i = 0; i < num_test_points; ++i) {
        const struct Point p = test_points[i];
        int intersections = 0;
        
        for (int j = 0; j < num_polygon_vertices; ++j) {
            const struct Point p1 = polygon[j];
            const struct Point p2 = polygon[(j + 1) % num_polygon_vertices];

            // Ray-casting algorithm
            // Check if the point's y-coordinate is between the edge's y-coordinates
            if (((p1.y > p.y) != (p2.y > p.y)) &&
                 // and the point is to the left of the line segment
                 (p.x < (p2.x - p1.x) * (p.y - p1.y) / (p2.y - p1.y) + p1.x)) {
                intersections++;
            }
        }

        // If the number of intersections is odd, the point is inside
        if (intersections % 2 != 0) {
            points_inside++;
        }
    }
    final_result = points_inside;
}

void cleanup() {
    free(polygon);
    free(test_points);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
