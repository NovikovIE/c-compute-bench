#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- BEGIN: Mersenne Twister (Do Not Modify) ---
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
// --- END: Mersenne Twister ---

// --- Data Structures ---
typedef struct {
    double x, y;
} Point;

typedef struct {
    Point* vertices;
    int num_vertices;
} Polygon;

// --- Global Data ---
static int num_polygons_A, num_polygons_B;
static Polygon* polygons_A;
static Polygon* polygons_B;
static long long final_result = 0;

// --- Polygon Generation Helpers ---
static Point g_centroid_for_sort;

int compare_angles(const void* a, const void* b) {
    Point* p1 = (Point*)a;
    Point* p2 = (Point*)b;
    double angle1 = atan2(p1->y - g_centroid_for_sort.y, p1->x - g_centroid_for_sort.x);
    double angle2 = atan2(p2->y - g_centroid_for_sort.y, p2->x - g_centroid_for_sort.x);
    if (angle1 < angle2) return -1;
    if (angle1 > angle2) return 1;
    return 0;
}

void generate_convex_polygon(Polygon* poly, int avg_vertices) {
    int num_vertices = avg_vertices + (mt_rand() % 7) - 3;
    if (num_vertices < 3) num_vertices = 3;
    poly->num_vertices = num_vertices;
    poly->vertices = (Point*)malloc(num_vertices * sizeof(Point));

    double center_x = (double)mt_rand() / UINT32_MAX * 1000.0;
    double center_y = (double)mt_rand() / UINT32_MAX * 1000.0;
    double radius = (double)mt_rand() / UINT32_MAX * 50.0 + 20.0;

    g_centroid_for_sort = (Point){.x = 0, .y = 0};
    for (int i = 0; i < num_vertices; ++i) {
        poly->vertices[i].x = center_x + ((double)mt_rand() / UINT32_MAX - 0.5) * 2 * radius;
        poly->vertices[i].y = center_y + ((double)mt_rand() / UINT32_MAX - 0.5) * 2 * radius;
        g_centroid_for_sort.x += poly->vertices[i].x;
        g_centroid_for_sort.y += poly->vertices[i].y;
    }
    g_centroid_for_sort.x /= num_vertices;
    g_centroid_for_sort.y /= num_vertices;

    qsort(poly->vertices, num_vertices, sizeof(Point), compare_angles);
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_polygons_A avg_vertices_A num_polygons_B avg_vertices_B seed\n", argv[0]);
        exit(1);
    }

    num_polygons_A = atoi(argv[1]);
    int avg_vertices_A = atoi(argv[2]);
    num_polygons_B = atoi(argv[3]);
    int avg_vertices_B = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    polygons_A = (Polygon*)malloc(num_polygons_A * sizeof(Polygon));
    for (int i = 0; i < num_polygons_A; ++i) {
        generate_convex_polygon(&polygons_A[i], avg_vertices_A);
    }

    polygons_B = (Polygon*)malloc(num_polygons_B * sizeof(Polygon));
    for (int i = 0; i < num_polygons_B; ++i) {
        generate_convex_polygon(&polygons_B[i], avg_vertices_B);
    }
}

void cleanup() {
    for (int i = 0; i < num_polygons_A; ++i) {
        free(polygons_A[i].vertices);
    }
    free(polygons_A);

    for (int i = 0; i < num_polygons_B; ++i) {
        free(polygons_B[i].vertices);
    }
    free(polygons_B);
}

// --- Core Computation: Sutherland-Hodgman Polygon Clipping ---
int is_inside(Point p, Point edge_start, Point edge_end) {
    return (edge_end.x - edge_start.x) * (p.y - edge_start.y) > (edge_end.y - edge_start.y) * (p.x - edge_start.x);
}

Point get_intersection_point(Point s1, Point e1, Point s2, Point e2) {
    double dx1 = e1.x - s1.x;
    double dy1 = e1.y - s1.y;
    double dx2 = e2.x - s2.x;
    double dy2 = e2.y - s2.y;
    double denominator = dx1 * dy2 - dy1 * dx2;
    // Assumes non-parallel lines, as convex polygon edges won't be parallel in a problematic way.
    double t = ((s2.x - s1.x) * dy2 - (s2.y - s1.y) * dx2) / denominator;
    return (Point){s1.x + t * dx1, s1.y + t * dy1};
}

Polygon clip_polygon(Polygon subject, Point clip_p1, Point clip_p2) {
    // FIX: The original allocation size was (subject.num_vertices + 1), which was too small
    // and could lead to a buffer overflow, causing heap corruption and a SIGABRT.
    // A safe upper bound for the output is 2 * num_vertices. We also ensure a minimum allocation size.
    int max_output_vertices = subject.num_vertices * 2;
    if (max_output_vertices < 4) {
        max_output_vertices = 4;
    }
    Point* output_list = (Point*)malloc(max_output_vertices * sizeof(Point));
    int output_count = 0;

    // FIX: The original code would access subject.vertices[-1] if subject.num_vertices was 0,
    // leading to a crash. This check prevents the out-of-bounds read.
    if (subject.num_vertices > 0) {
        Point s = subject.vertices[subject.num_vertices - 1];
        for (int j = 0; j < subject.num_vertices; ++j) {
            Point e = subject.vertices[j];
            int s_inside = is_inside(s, clip_p1, clip_p2);
            int e_inside = is_inside(e, clip_p1, clip_p2);

            if (s_inside && e_inside) {
                output_list[output_count++] = e;
            } else if (s_inside && !e_inside) {
                output_list[output_count++] = get_intersection_point(s, e, clip_p1, clip_p2);
            } else if (!s_inside && e_inside) {
                output_list[output_count++] = get_intersection_point(s, e, clip_p1, clip_p2);
                output_list[output_count++] = e;
            }
            s = e;
        }
    }

    free(subject.vertices);
    return (Polygon){output_list, output_count};
}

Polygon intersect(const Polygon* subject, const Polygon* clipper) {
    // Create a working copy of the subject polygon
    Polygon temp_poly;
    temp_poly.num_vertices = subject->num_vertices;
    temp_poly.vertices = (Point*)malloc(temp_poly.num_vertices * sizeof(Point));
    for(int i = 0; i < temp_poly.num_vertices; ++i) {
        temp_poly.vertices[i] = subject->vertices[i];
    }

    for (int i = 0; i < clipper->num_vertices; ++i) {
        Point p1 = clipper->vertices[i];
        Point p2 = clipper->vertices[(i + 1) % clipper->num_vertices];
        temp_poly = clip_polygon(temp_poly, p1, p2);
    }

    return temp_poly;
}

void run_computation() {
    long long total_intersection_vertices = 0;
    for (int i = 0; i < num_polygons_A; ++i) {
        for (int j = 0; j < num_polygons_B; ++j) {
            Polygon result = intersect(&polygons_A[i], &polygons_B[j]);
            total_intersection_vertices += result.num_vertices;
            if (result.vertices) {
                free(result.vertices);
            }
        }
    }
    final_result = total_intersection_vertices;
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

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
