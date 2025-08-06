#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Start of Mersenne Twister (Do Not Modify) ---
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

// Benchmark-specific data structures and globals
typedef struct {
    double x;
    double y;
} Point;

int num_polygon_vertices;
Point* polygon_vertices;
int* triangles; // Stores indices: (v1, v2, v3), (v1, v2, v3), ...
int* active_indices; // Helper array for tracking vertices
long long final_result;

// Function to handle allocation errors
void handle_alloc_error(void* ptr, const char* name) {
    if (!ptr) {
        fprintf(stderr, "FATAL: Memory allocation failed for %s\n", name);
        exit(1);
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_polygon_vertices> <seed>\n", argv[0]);
        exit(1);
    }

    num_polygon_vertices = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_polygon_vertices < 3) {
        fprintf(stderr, "FATAL: Number of vertices must be at least 3.\n");
        exit(1);
    }

    mt_seed(seed);

    polygon_vertices = (Point*)malloc(num_polygon_vertices * sizeof(Point));
    handle_alloc_error(polygon_vertices, "polygon_vertices");

    // We will generate N-2 triangles
    triangles = (int*)malloc((num_polygon_vertices - 2) * 3 * sizeof(int));
    handle_alloc_error(triangles, "triangles");
    
    active_indices = (int*)malloc(num_polygon_vertices * sizeof(int));
    handle_alloc_error(active_indices, "active_indices");

    // Generate a simple, convex polygon by generating points on a circle with radius noise
    // This ensures the ear-clipping algorithm can be simplified for the benchmark.
    double center_x = 5000.0;
    double center_y = 5000.0;
    double radius_base = 4500.0;

    for (int i = 0; i < num_polygon_vertices; ++i) {
        double angle = 2.0 * M_PI * i / num_polygon_vertices;
        double radius_noise = (double)(mt_rand() % 1000) / 100.0 - 5.0;
        double current_radius = radius_base + radius_noise;
        polygon_vertices[i].x = center_x + current_radius * cos(angle);
        polygon_vertices[i].y = center_y + current_radius * sin(angle);
    }

    final_result = 0;
}

void run_computation() {
    if (num_polygon_vertices < 3) return;

    // Initialize the list of active vertex indices
    for (int i = 0; i < num_polygon_vertices; ++i) {
        active_indices[i] = i;
    }

    int num_active = num_polygon_vertices;
    int triangle_count = 0;

    // This is a simplified O(n^2) ear-clipping algorithm for CONVEX polygons.
    // We repeatedly clip the ear formed by the first three active vertices.
    while (num_active > 3) {
        int v0_original_idx = active_indices[0];
        int v1_original_idx = active_indices[1];
        int v2_original_idx = active_indices[2];

        // Store the new triangle using original vertex indices
        triangles[triangle_count * 3 + 0] = v0_original_idx;
        triangles[triangle_count * 3 + 1] = v1_original_idx;
        triangles[triangle_count * 3 + 2] = v2_original_idx;
        triangle_count++;

        // Accumulate a result to prevent dead code elimination.
        // We use the coordinates of the first vertex of each new triangle.
        final_result += (long long)(polygon_vertices[v0_original_idx].x);

        // Remove the ear tip (the second vertex in our current list) by shifting the array.
        // This is the O(n) part of the loop, making the whole algorithm O(n^2).
        for (int i = 1; i < num_active - 1; ++i) {
            active_indices[i] = active_indices[i + 1];
        }
        num_active--;
    }

    // Add the final remaining triangle
    triangles[triangle_count * 3 + 0] = active_indices[0];
    triangles[triangle_count * 3 + 1] = active_indices[1];
    triangles[triangle_count * 3 + 2] = active_indices[2];
    final_result += (long long)(polygon_vertices[active_indices[0]].x);
}

void cleanup() {
    free(polygon_vertices);
    free(triangles);
    free(active_indices);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%lld\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
