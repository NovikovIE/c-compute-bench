/*
 * Benchmark: visibility_graph_construction
 * Theme:     Computational Geometry
 * Author:    C-Expert
 * Description: Computes the visibility graph for a set of non-overlapping simple polygonal obstacles.
 *              The visibility graph connects any two vertices from the set of all obstacle vertices
 *              if the line segment between them does not intersect the interior of any obstacle.
 *              The benchmark's complexity is dominated by O(V^3) segment intersection tests,
 *              where V is the total number of vertices.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) generator (DO NOT MODIFY) ---
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

// --- Benchmark Data Structures ---
typedef struct {
    double x, y;
} Point;

typedef struct {
    int start_vertex_index;
    int num_vertices;
} Polygon;

// --- Global Variables ---
static int num_obstacles_g;
static int total_vertices_g;

static Point* vertices_g;
static Polygon* obstacles_g;
static long long final_result_g;


// --- Helper functions for geometry ---

// To find orientation of ordered triplet (p, q, r).
// Returns: 0 (collinear), 1 (clockwise), 2 (counter-clockwise)
static int orientation(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (fabs(val) < 1e-10) return 0; // Use a small epsilon for float comparison
    return (val > 0) ? 1 : 2;
}

// Checks if point q lies on line segment 'pr', given that p, q, r are collinear.
static int on_segment(Point p, Point q, Point r) {
    if (q.x <= fmax(p.x, r.x) && q.x >= fmin(p.x, r.x) &&
        q.y <= fmax(p.y, r.y) && q.y >= fmin(p.y, r.y)) {
        return 1;
    }
    return 0;
}

// Returns 1 if line segment 'p1q1' and 'p2q2' intersect.
static int segments_intersect(Point p1, Point q1, Point p2, Point q2) {
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    // General case: segments cross each other
    if (o1 != o2 && o3 != o4) {
        return 1;
    }

    // Special Cases for collinear points
    if (o1 == 0 && on_segment(p1, p2, q1)) return 1;
    if (o2 == 0 && on_segment(p1, q2, q1)) return 1;
    if (o3 == 0 && on_segment(p2, p1, q2)) return 1;
    if (o4 == 0 && on_segment(p2, q1, q2)) return 1;

    return 0; // No intersection
}


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_obstacles> <total_vertices> <seed>\n", argv[0]);
        exit(1);
    }

    num_obstacles_g = atoi(argv[1]);
    total_vertices_g = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    
    if (num_obstacles_g <= 0 || total_vertices_g < num_obstacles_g * 3) {
        fprintf(stderr, "Error: num_obstacles must be > 0 and total_vertices must be at least 3 * num_obstacles\n");
        exit(1);
    }

    mt_seed(seed);

    vertices_g = (Point*)malloc(total_vertices_g * sizeof(Point));
    obstacles_g = (Polygon*)malloc(num_obstacles_g * sizeof(Polygon));
    if (!vertices_g || !obstacles_g) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    int vertices_per_obstacle = total_vertices_g / num_obstacles_g;
    int remaining_vertices = total_vertices_g % num_obstacles_g;
    int current_vertex_index = 0;

    double world_size = 10000.0;
    int grid_dim = (int)ceil(sqrt((double)num_obstacles_g));
    double cell_size = world_size / grid_dim;

    for (int i = 0; i < num_obstacles_g; ++i) {
        obstacles_g[i].start_vertex_index = current_vertex_index;
        obstacles_g[i].num_vertices = vertices_per_obstacle + (i < remaining_vertices ? 1 : 0);

        int grid_x = i % grid_dim;
        int grid_y = i / grid_dim;

        double offset_x = grid_x * cell_size;
        double offset_y = grid_y * cell_size;
        
        Point centroid = {
            .x = offset_x + (mt_rand() / (double)UINT32_MAX) * cell_size * 0.8 + cell_size * 0.1,
            .y = offset_y + (mt_rand() / (double)UINT32_MAX) * cell_size * 0.8 + cell_size * 0.1
        };

        double base_radius = cell_size / 5.0;
        double radius_flutter = cell_size / 10.0;

        for (int j = 0; j < obstacles_g[i].num_vertices; ++j) {
            double angle = 2.0 * M_PI * j / obstacles_g[i].num_vertices;
            double radius = base_radius + (mt_rand() / (double)UINT32_MAX) * radius_flutter;
            vertices_g[current_vertex_index].x = centroid.x + radius * cos(angle);
            vertices_g[current_vertex_index].y = centroid.y + radius * sin(angle);
            current_vertex_index++;
        }
    }
}

void run_computation() {
    long long visible_edge_count = 0;

    for (int i = 0; i < total_vertices_g; i++) {
        for (int j = i + 1; j < total_vertices_g; j++) {
            Point p1 = vertices_g[i];
            Point p2 = vertices_g[j];
            int is_blocked = 0;

            for (int k = 0; k < num_obstacles_g; k++) {
                Polygon obs = obstacles_g[k];
                for (int v_offset = 0; v_offset < obs.num_vertices; v_offset++) {
                    int v_idx1 = obs.start_vertex_index + v_offset;
                    int v_idx2 = obs.start_vertex_index + (v_offset + 1) % obs.num_vertices;
                    
                    // A segment is not blocked by non-proper intersections at its endpoints.
                    if (i == v_idx1 || i == v_idx2 || j == v_idx1 || j == v_idx2) {
                        continue;
                    }

                    Point a = vertices_g[v_idx1];
                    Point b = vertices_g[v_idx2];
                    
                    if (segments_intersect(p1, p2, a, b)) {
                        is_blocked = 1;
                        break;
                    }
                }
                if (is_blocked) {
                    break;
                }
            }

            if (!is_blocked) {
                visible_edge_count++;
            }
        }
    }
    final_result_g = visible_edge_count;
}

void cleanup() {
    free(vertices_g);
    free(obstacles_g);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    
    setup_benchmark(argc, argv);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%lld\n", final_result_g);
    fprintf(stderr, "%.6f", time_taken);
    
    return 0;
}
