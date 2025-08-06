#define _GNU_SOURCE
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
typedef struct {
    double x;
    double y;
} Point;

// Global data structure
struct {
    int num_points;
    Point *points;
    int *hull; // Stores indices of points on the convex hull
    int hull_size;
    int final_result;
} g_data;

// --- Benchmark Functions ---

// Helper to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <num_hull_vertices> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_points = atoi(argv[1]);
    int num_hull_vertices_target = atoi(argv[2]);
    uint32_t seed = (uint32_t)atol(argv[3]);

    if (g_data.num_points <= 3 || num_hull_vertices_target < 3 || num_hull_vertices_target > g_data.num_points) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.points = (Point *)malloc(g_data.num_points * sizeof(Point));
    g_data.hull = (int *)malloc(g_data.num_points * sizeof(int)); // Max possible hull size is num_points
    if (!g_data.points || !g_data.hull) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }
    
    double hull_radius = 10000.0;
    double inner_radius = hull_radius * 0.95;
    
    // Generate num_hull_vertices_target points on a circle to ensure they form the hull
    for (int i = 0; i < num_hull_vertices_target; ++i) {
        double angle = 2.0 * M_PI * i / num_hull_vertices_target;
        double r = hull_radius * (0.99 + 0.02 * rand_double()); // Add jitter
        g_data.points[i].x = r * cos(angle);
        g_data.points[i].y = r * sin(angle);
    }

    // Generate remaining points inside the inner circle
    for (int i = num_hull_vertices_target; i < g_data.num_points; ++i) {
        double r = inner_radius * sqrt(rand_double()); // sqrt for uniform distribution over area
        double angle = 2.0 * M_PI * rand_double();
        g_data.points[i].x = r * cos(angle);
        g_data.points[i].y = r * sin(angle);
    }

    // Shuffle all points to mix hull and inner points
    for (int i = g_data.num_points - 1; i > 0; i--) {
        int j = mt_rand() % (i + 1);
        Point temp = g_data.points[i];
        g_data.points[i] = g_data.points[j];
        g_data.points[j] = temp;
    }
    
    g_data.final_result = 0;
    g_data.hull_size = 0;
}

void cleanup() {
    free(g_data.points);
    free(g_data.hull);
}

// Helper function to find orientation of ordered triplet (p, q, r).
// Returns: 0 (collinear), 1 (clockwise), 2 (counter-clockwise)
static inline int orientation(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val < 1e-9 && val > -1e-9) return 0; // Collinear
    return (val > 0) ? 1 : 2; // 1 for Clockwise, 2 for Counter-clockwise
}

// Helper to calculate squared distance between two points
static inline double dist_sq(Point p1, Point p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return dx * dx + dy * dy;
}

void run_computation() {
    int n = g_data.num_points;
    if (n < 3) return;

    // 1. Find the leftmost point (guaranteed to be on the hull)
    int leftmost_idx = 0;
    for (int i = 1; i < n; i++) {
        if (g_data.points[i].x < g_data.points[leftmost_idx].x) {
            leftmost_idx = i;
        } else if (fabs(g_data.points[i].x - g_data.points[leftmost_idx].x) < 1e-9 && g_data.points[i].y < g_data.points[leftmost_idx].y) {
            leftmost_idx = i;
        }
    }

    // 2. Start from the leftmost point, keep moving counterclockwise
    //    until we reach the start point again. This is the Jarvis March.
    int p_idx = leftmost_idx;
    int q_idx;
    int count = 0;
    do {
        g_data.hull[count++] = p_idx;

        // Find the most counter-clockwise point 'q' with respect to 'p'.
        q_idx = (p_idx + 1) % n;
        for (int i = 0; i < n; i++) {
            // If point i is more counter-clockwise than current q, update q_idx
            int o = orientation(g_data.points[p_idx], g_data.points[q_idx], g_data.points[i]);
            if (o == 2) { // 2 is Counter-clockwise
                q_idx = i;
            }
            // If p, q, i are collinear, pick the point furthest from p.
            else if (o == 0 && dist_sq(g_data.points[p_idx], g_data.points[i]) > dist_sq(g_data.points[p_idx], g_data.points[q_idx])) {
                q_idx = i;
            }
        }

        p_idx = q_idx;

    } while (p_idx != leftmost_idx);

    g_data.hull_size = count;

    // Calculate a final result (checksum of hull point indices) to prevent dead code elimination.
    int sum = 0;
    for (int i = 0; i < g_data.hull_size; i++) {
        sum += g_data.hull[i];
    }
    g_data.final_result = sum;
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

    printf("%d\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
