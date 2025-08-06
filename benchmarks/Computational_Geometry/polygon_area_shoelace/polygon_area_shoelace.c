#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (Do Not Modify - Include This Verbatim) ---
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

// Benchmark data structures and parameters
typedef struct {
    double x;
    double y;
} Point;

int num_polygon_vertices;
Point *polygon;
double final_area; // To hold the result

// Generates a random double between 0.0 and max_val
double rand_double(double max_val) {
    return ((double)mt_rand() / (double)UINT32_MAX) * max_val;
}

// Setup function: parses arguments, allocates memory, generates data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_polygon_vertices> <seed>\n", argv[0]);
        exit(1);
    }

    num_polygon_vertices = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_polygon_vertices <= 2) {
        fprintf(stderr, "Error: Number of vertices must be at least 3.\n");
        exit(1);
    }

    mt_seed(seed);

    polygon = (Point *)malloc(num_polygon_vertices * sizeof(Point));
    if (polygon == NULL) {
        fprintf(stderr, "Failed to allocate memory for the polygon.\n");
        exit(1);
    }

    // Generate random vertices for the polygon.
    // Coordinates will be in the range [0.0, 1000.0]
    for (int i = 0; i < num_polygon_vertices; i++) {
        polygon[i].x = rand_double(1000.0);
        polygon[i].y = rand_double(1000.0);
    }
}

// Computation function: calculates the polygon area using the Shoelace formula
void run_computation() {
    double area = 0.0;
    int n = num_polygon_vertices;

    // Shoelace formula: 0.5 * |sum(x_i * y_{i+1} - x_{i+1} * y_i)|
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n; // wrap around for the last vertex
        area += (polygon[i].x * polygon[j].y) - (polygon[j].x * polygon[i].y);
    }

    final_area = 0.5 * fabs(area);
}

// Cleanup function: frees allocated memory
void cleanup() {
    free(polygon);
}

// Main driver
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    // Calculate time taken in seconds
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%.6f\n", final_area);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
