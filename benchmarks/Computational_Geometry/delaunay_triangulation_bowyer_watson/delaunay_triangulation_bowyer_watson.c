#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

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

#define BOX_SIZE 10000.0

// --- Data Structures ---
typedef struct {
    double x, y;
} Point;

typedef struct {
    int p1_idx, p2_idx, p3_idx;
} Triangle;

typedef struct {
    int p1_idx, p2_idx;
} Edge;

// --- Global Data ---
struct BenchmarkData {
    int num_points;
    int total_points; // num_points + 3 super triangle vertices
    Point* points;
    Triangle* triangles;
    int num_triangles;
    int triangles_capacity;
    long long final_checksum;
} g_data;

// --- Helper Functions for Computation ---

// Comparison for sorting edges
int compare_edges(const void* a, const void* b) {
    const Edge* e1 = (const Edge*)a;
    const Edge* e2 = (const Edge*)b;
    if (e1->p1_idx < e2->p1_idx) return -1;
    if (e1->p1_idx > e2->p1_idx) return 1;
    if (e1->p2_idx < e2->p2_idx) return -1;
    if (e1->p2_idx > e2->p2_idx) return 1;
    return 0;
}

// Checks if a point is inside the circumcircle of a triangle.
// Uses a robust determinant method.
int is_in_circumcircle(const Point p, const Triangle t) {
    Point p1 = g_data.points[t.p1_idx];
    Point p2 = g_data.points[t.p2_idx];
    Point p3 = g_data.points[t.p3_idx];

    double dx1 = p1.x - p.x, dy1 = p1.y - p.y;
    double dx2 = p2.x - p.x, dy2 = p2.y - p.y;
    double dx3 = p3.x - p.x, dy3 = p3.y - p.y;

    double det = (dx1 * dx1 + dy1 * dy1) * (dx2 * dy3 - dx3 * dy2) -
                 (dx2 * dx2 + dy2 * dy2) * (dx1 * dy3 - dx3 * dy1) +
                 (dx3 * dx3 + dy3 * dy3) * (dx1 * dy2 - dx2 * dy1);
    
    // The super-triangle is defined in CCW order. For a CCW triangle, a positive determinant means the point is inside.
    return det > 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_points> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_points = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (g_data.num_points <= 0) {
        fprintf(stderr, "Error: num_points must be positive.\n");
        exit(1);
    }
    
    mt_seed(seed);

    g_data.total_points = g_data.num_points + 3;
    g_data.points = (Point*)malloc(g_data.total_points * sizeof(Point));
    if (!g_data.points) {
        fprintf(stderr, "Failed to allocate memory for points.\n");
        exit(1);
    }

    // Generate random points within a box
    for (int i = 0; i < g_data.num_points; i++) {
        g_data.points[i].x = ((double)mt_rand() / (double)UINT32_MAX) * BOX_SIZE;
        g_data.points[i].y = ((double)mt_rand() / (double)UINT32_MAX) * BOX_SIZE;
    }
    
    // Create a super-triangle that encloses all points
    double margin = BOX_SIZE * 3.0;
    g_data.points[g_data.num_points + 0] = (Point){-margin + BOX_SIZE/2.0, -margin/2.0};
    g_data.points[g_data.num_points + 1] = (Point){BOX_SIZE + margin, -margin/2.0};
    g_data.points[g_data.num_points + 2] = (Point){BOX_SIZE / 2.0, BOX_SIZE + margin};

    // Initialize triangulation with the super-triangle
    g_data.triangles_capacity = 2 * g_data.num_points + 1;
    g_data.triangles = (Triangle*)malloc(g_data.triangles_capacity * sizeof(Triangle));
    if (!g_data.triangles) {
        fprintf(stderr, "Failed to allocate memory for triangles.\n");
        exit(1);
    }
    
    g_data.triangles[0] = (Triangle){g_data.num_points, g_data.num_points + 1, g_data.num_points + 2};
    g_data.num_triangles = 1;
    g_data.final_checksum = 0;
}

void run_computation() {
    int max_bad_triangles = 50; // Initial guess, will grow
    int* bad_triangles_idx = (int*)malloc(max_bad_triangles * sizeof(int));
    if (!bad_triangles_idx) exit(1);

    int max_polygon_edges = max_bad_triangles * 3;
    Edge* polygon_edges = (Edge*)malloc(max_polygon_edges * sizeof(Edge));
    if (!polygon_edges) { free(bad_triangles_idx); exit(1); }

    for (int i = 0; i < g_data.num_points; ++i) {
        Point current_point = g_data.points[i];
        int num_bad_triangles = 0;

        // Find bad triangles
        for (int j = 0; j < g_data.num_triangles; ++j) {
            if (is_in_circumcircle(current_point, g_data.triangles[j])) {
                if (num_bad_triangles >= max_bad_triangles) {
                    max_bad_triangles *= 2;
                    bad_triangles_idx = (int*)realloc(bad_triangles_idx, max_bad_triangles * sizeof(int));
                    if (!bad_triangles_idx) exit(1);
                }
                bad_triangles_idx[num_bad_triangles++] = j;
            }
        }
        
        // Collect edges of bad triangles
        int num_polygon_edges = 0;
        if (num_bad_triangles * 3 > max_polygon_edges) {
            max_polygon_edges = num_bad_triangles * 3;
            polygon_edges = (Edge*)realloc(polygon_edges, max_polygon_edges * sizeof(Edge));
            if (!polygon_edges) exit(1);
        }

        for (int j = 0; j < num_bad_triangles; ++j) {
            Triangle t = g_data.triangles[bad_triangles_idx[j]];
            Edge e1 = {t.p1_idx, t.p2_idx}, e2 = {t.p2_idx, t.p3_idx}, e3 = {t.p3_idx, t.p1_idx};
            // Normalize edges for consistent comparison
            if (e1.p1_idx > e1.p2_idx) { int temp = e1.p1_idx; e1.p1_idx = e1.p2_idx; e1.p2_idx = temp; }
            if (e2.p1_idx > e2.p2_idx) { int temp = e2.p1_idx; e2.p1_idx = e2.p2_idx; e2.p2_idx = temp; }
            if (e3.p1_idx > e3.p2_idx) { int temp = e3.p1_idx; e3.p1_idx = e3.p2_idx; e3.p2_idx = temp; }
            polygon_edges[num_polygon_edges++] = e1;
            polygon_edges[num_polygon_edges++] = e2;
            polygon_edges[num_polygon_edges++] = e3;
        }

        // Remove bad triangles by replacing them with triangles from the end of the list
        // Sort indices descending to avoid invalidation
        for(int m=0; m<num_bad_triangles; m++) {
            for(int n=m+1; n<num_bad_triangles; n++) {
                if(bad_triangles_idx[m] < bad_triangles_idx[n]) {
                    int temp = bad_triangles_idx[m];
                    bad_triangles_idx[m] = bad_triangles_idx[n];
                    bad_triangles_idx[n] = temp;
                }
            }
        }
        for (int j = 0; j < num_bad_triangles; j++) {
            g_data.triangles[bad_triangles_idx[j]] = g_data.triangles[--g_data.num_triangles];
        }

        // Find unique edges which form the polygonal hole
        qsort(polygon_edges, num_polygon_edges, sizeof(Edge), compare_edges);
        
        for (int j = 0; j < num_polygon_edges; ) {
            if (j + 1 < num_polygon_edges && polygon_edges[j].p1_idx == polygon_edges[j+1].p1_idx && polygon_edges[j].p2_idx == polygon_edges[j+1].p2_idx) {
                j += 2; // This edge is shared by two bad triangles, so it's internal. Skip.
            } else {
                // This is a boundary edge. Form a new triangle with the current point.
                if (g_data.num_triangles >= g_data.triangles_capacity) {
                    g_data.triangles_capacity *= 2;
                    g_data.triangles = (Triangle*)realloc(g_data.triangles, g_data.triangles_capacity * sizeof(Triangle));
                    if (!g_data.triangles) exit(1);
                }
                g_data.triangles[g_data.num_triangles++] = (Triangle){polygon_edges[j].p1_idx, polygon_edges[j].p2_idx, i};
                j++;
            }
        }
    }
    
    free(bad_triangles_idx);
    free(polygon_edges);

    // Final pass: remove any triangles that include vertices from the super-triangle
    int write_idx = 0;
    for (int i = 0; i < g_data.num_triangles; ++i) {
        Triangle t = g_data.triangles[i];
        if (t.p1_idx < g_data.num_points && t.p2_idx < g_data.num_points && t.p3_idx < g_data.num_points) {
            if (i != write_idx) {
                g_data.triangles[write_idx] = t;
            }
            write_idx++;
        }
    }
    g_data.num_triangles = write_idx;

    // Calculate a checksum to prevent dead code elimination and provide a result
    for (int i = 0; i < g_data.num_triangles; ++i) {
        g_data.final_checksum += g_data.triangles[i].p1_idx + g_data.triangles[i].p2_idx + g_data.triangles[i].p3_idx;
    }
}

void cleanup() {
    free(g_data.points);
    free(g_data.triangles);
}

// --- Main Function ---
int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print result to stdout
    printf("%lld\n", g_data.final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
