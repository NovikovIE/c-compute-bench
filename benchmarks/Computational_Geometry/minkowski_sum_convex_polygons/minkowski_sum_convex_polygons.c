#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// START of Mersenne Twister (Do Not Modify)
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
// END of Mersenne Twister

// Benchmark-specific data structures and globals
typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    int num_vertices_A;
    int num_vertices_B;
    Point *polygon_A;
    Point *polygon_B;
    Point *minkowski_sum;
    int num_vertices_sum;
    int final_result; // To prevent dead code elimination (sum of vertex indices)
} BenchmarkData;

static BenchmarkData g_data;
static Point g_centroid; // Global for qsort comparison function

// Helper to generate a random double
double rand_double(double min, double max) {
    return min + (mt_rand() / (double)UINT32_MAX) * (max - min);
}

// Comparison function for qsort to sort points by angle around a centroid
int compare_points(const void* a, const void* b) {
    Point* p1 = (Point*)a;
    Point* p2 = (Point*)b;
    double angle1 = atan2(p1->y - g_centroid.y, p1->x - g_centroid.x);
    double angle2 = atan2(p2->y - g_centroid.y, p2->x - g_centroid.x);
    if (angle1 < angle2) return -1;
    if (angle1 > angle2) return 1;
    return 0;
}

// Generates a simple polygon by sorting points around their centroid
void generate_polygon(Point* polygon, int num_vertices) {
    if (num_vertices <= 0) return;

    g_centroid.x = 0;
    g_centroid.y = 0;

    for (int i = 0; i < num_vertices; ++i) {
        polygon[i].x = rand_double(-1000.0, 1000.0);
        polygon[i].y = rand_double(-1000.0, 1000.0);
        g_centroid.x += polygon[i].x;
        g_centroid.y += polygon[i].y;
    }
    g_centroid.x /= num_vertices;
    g_centroid.y /= num_vertices;

    qsort(polygon, num_vertices, sizeof(Point), compare_points);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_vertices_A num_vertices_B seed\n", argv[0]);
        exit(1);
    }

    g_data.num_vertices_A = atoi(argv[1]);
    g_data.num_vertices_B = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_vertices_A <= 2 || g_data.num_vertices_B <= 2) {
        fprintf(stderr, "Error: Number of vertices for each polygon must be greater than 2.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.polygon_A = (Point*)malloc(g_data.num_vertices_A * sizeof(Point));
    g_data.polygon_B = (Point*)malloc(g_data.num_vertices_B * sizeof(Point));
    g_data.minkowski_sum = (Point*)malloc((g_data.num_vertices_A + g_data.num_vertices_B) * sizeof(Point));

    if (!g_data.polygon_A || !g_data.polygon_B || !g_data.minkowski_sum) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    generate_polygon(g_data.polygon_A, g_data.num_vertices_A);
    generate_polygon(g_data.polygon_B, g_data.num_vertices_B);

    g_data.final_result = 0;
}

void run_computation() {
    int nA = g_data.num_vertices_A;
    int nB = g_data.num_vertices_B;
    Point* pA = g_data.polygon_A;
    Point* pB = g_data.polygon_B;
    Point* sum = g_data.minkowski_sum;

    sum[0].x = pA[0].x + pB[0].x;
    sum[0].y = pA[0].y + pB[0].y;

    int i = 0, j = 0, k = 0;
    int max_edges = nA + nB - 1;

    while (i < nA && j < nB) {
        if (k >= max_edges) break;
        Point pA_curr = pA[i];
        Point pA_next = pA[(i + 1) % nA];
        Point pB_curr = pB[j];
        Point pB_next = pB[(j + 1) % nB];

        double edgeA_x = pA_next.x - pA_curr.x;
        double edgeA_y = pA_next.y - pA_curr.y;
        double edgeB_x = pB_next.x - pB_curr.x;
        double edgeB_y = pB_next.y - pB_curr.y;

        double angleA = atan2(edgeA_y, edgeA_x);
        double angleB = atan2(edgeB_y, edgeB_x);

        if (angleA < angleB) {
            sum[k+1].x = sum[k].x + edgeA_x;
            sum[k+1].y = sum[k].y + edgeA_y;
            i++;
        } else {
            sum[k+1].x = sum[k].x + edgeB_x;
            sum[k+1].y = sum[k].y + edgeB_y;
            j++;
        }
        k++;
    }

    while (i < nA) {
        if (k >= max_edges) break;
        Point pA_curr = pA[i];
        Point pA_next = pA[(i + 1) % nA];
        sum[k+1].x = sum[k].x + (pA_next.x - pA_curr.x);
        sum[k+1].y = sum[k].y + (pA_next.y - pA_curr.y);
        i++;
        k++;
    }

    while (j < nB) {
        if (k >= max_edges) break;
        Point pB_curr = pB[j];
        Point pB_next = pB[(j + 1) % nB];
        sum[k+1].x = sum[k].x + (pB_next.x - pB_curr.x);
        sum[k+1].y = sum[k].y + (pB_next.y - pB_curr.y);
        j++;
        k++;
    }

    g_data.num_vertices_sum = k + 1;

    // Calculate area of the resulting polygon to prevent dead code elimination
    double area = 0.0;
    for (int v = 0; v < g_data.num_vertices_sum; ++v) {
        Point p1 = sum[v];
        Point p2 = sum[(v + 1) % g_data.num_vertices_sum];
        area += (p1.x * p2.y - p2.x * p1.y);
    }
    g_data.final_result = (int)(fabs(area * 0.5));
}

void cleanup() {
    free(g_data.polygon_A);
    free(g_data.polygon_B);
    free(g_data.minkowski_sum);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    printf("%d\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
