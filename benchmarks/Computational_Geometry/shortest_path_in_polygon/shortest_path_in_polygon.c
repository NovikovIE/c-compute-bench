#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (MT19937) --- (DO NOT MODIFY)
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

double mt_rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}
// --- End of Mersenne Twister ---

// --- Data Structures ---
typedef struct {
    double x, y;
} Point;

// --- Global Benchmark State ---
static int N;                 // Number of polygon vertices
static int start_vertex;
static int end_vertex;
static Point *polygon;
static double *visibility_graph; // Adjacency matrix
static double path_length;

// --- Helper Functions for Geometry ---
double distance_sq(Point p1, Point p2) {
    return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

// Given three collinear points p, q, r, the function checks if point q lies on line segment 'pr'
static bool on_segment(Point p, Point q, Point r) {
    return (q.x <= fmax(p.x, r.x) && q.x >= fmin(p.x, r.x) &&
            q.y <= fmax(p.y, r.y) && q.y >= fmin(p.y, r.y));
}

// To find orientation of ordered triplet (p, q, r).
// 0 --> p, q and r are collinear
// 1 --> Clockwise
// 2 --> Counterclockwise
static int orientation(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (fabs(val) < 1e-9) return 0; // collinear
    return (val > 0) ? 1 : 2;      // clock or counterclock wise
}

// The main function that returns true if line segment 'p1q1' and 'p2q2' intersect.
static bool segments_intersect(Point p1, Point q1, Point p2, Point q2) {
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    if (o1 != o2 && o3 != o4) return true;

    // Special Cases for collinear points
    if (o1 == 0 && on_segment(p1, p2, q1)) return true;
    if (o2 == 0 && on_segment(p1, q2, q1)) return true;
    if (o3 == 0 && on_segment(p2, p1, q2)) return true;
    if (o4 == 0 && on_segment(p2, q1, q2)) return true;

    return false;
}

// Returns true if the point p lies inside the polygon[] of size n.
static bool is_inside_polygon(Point p) {
    if (N < 3) return false;
    Point extreme = {DBL_MAX, p.y};
    int count = 0, i = 0;
    do {
        int next = (i + 1) % N;
        if (segments_intersect(polygon[i], polygon[next], p, extreme)) {
            if (orientation(polygon[i], p, polygon[next]) == 0) {
                return on_segment(polygon[i], p, polygon[next]);
            }
            count++;
        }
        i = next;
    } while (i != 0);
    return count & 1; // Same as count % 2 == 1
}


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_polygon_vertices start_vertex_index end_vertex_index seed\n", argv[0]);
        exit(1);
    }

    N = atoi(argv[1]);
    start_vertex = atoi(argv[2]);
    end_vertex = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);
    mt_seed(seed);

    if (N <= 2 || start_vertex >= N || end_vertex >= N || start_vertex < 0 || end_vertex < 0) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }

    polygon = (Point *)malloc(N * sizeof(Point));
    visibility_graph = (double *)malloc((size_t)N * N * sizeof(double));
    if (!polygon || !visibility_graph) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate a simple, star-shaped polygon
    double center_x = 0.0, center_y = 0.0;
    double base_radius = 200.0, random_range = 100.0;
    for (int i = 0; i < N; ++i) {
        double angle = 2.0 * M_PI * i / N;
        double radius = base_radius + random_range * mt_rand_double();
        polygon[i].x = center_x + radius * cos(angle);
        polygon[i].y = center_y + radius * sin(angle);
    }

    path_length = 0.0;
}

void run_computation() {
    // 1. Build visibility graph
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
             visibility_graph[i * N + j] = DBL_MAX;
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            bool is_visible = true;
            // An edge is visible if it doesn't properly intersect any other edge
            // and its midpoint is inside the polygon.
            for (int k = 0; k < N; ++k) {
                int k_next = (k + 1) % N;
                if (k == i || k == j || k_next == i || k_next == j) continue;
                if (segments_intersect(polygon[i], polygon[j], polygon[k], polygon[k_next])) {
                    is_visible = false;
                    break;
                }
            }

            if (is_visible) {
                // Midpoint check for segments that don't cross edges but could pass outside
                Point mid_point = {(polygon[i].x + polygon[j].x) / 2.0, (polygon[i].y + polygon[j].y) / 2.0};
                if (!is_inside_polygon(mid_point)) {
                    is_visible = false;
                }
            }
            
            // Polygon edges are always visible
            if (abs(i-j) == 1 || abs(i-j) == N-1) {
                is_visible = true;
            }

            if (is_visible) {
                double dist = sqrt(distance_sq(polygon[i], polygon[j]));
                visibility_graph[i * N + j] = dist;
                visibility_graph[j * N + i] = dist;
            }
        }
    }

    // 2. Run Dijkstra's algorithm
    double *dist = (double *)malloc(N * sizeof(double));
    bool *visited = (bool *)malloc(N * sizeof(bool));
    if (!dist || !visited) { exit(1); } // Should not happen

    for (int i = 0; i < N; ++i) {
        dist[i] = DBL_MAX;
        visited[i] = false;
    }

    dist[start_vertex] = 0.0;

    for (int count = 0; count < N - 1; ++count) {
        double min_dist = DBL_MAX;
        int u = -1;

        for (int v = 0; v < N; ++v) {
            if (!visited[v] && dist[v] <= min_dist) {
                min_dist = dist[v];
                u = v;
            }
        }

        if (u == -1 || u == end_vertex) break;

        visited[u] = true;

        for (int v = 0; v < N; ++v) {
            if (!visited[v] && visibility_graph[u * N + v] != DBL_MAX &&
                dist[u] != DBL_MAX && dist[u] + visibility_graph[u * N + v] < dist[v]) {
                dist[v] = dist[u] + visibility_graph[u * N + v];
            }
        }
    }

    path_length = dist[end_vertex];

    free(dist);
    free(visited);
}

void cleanup() {
    free(polygon);
    free(visibility_graph);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", path_length);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
