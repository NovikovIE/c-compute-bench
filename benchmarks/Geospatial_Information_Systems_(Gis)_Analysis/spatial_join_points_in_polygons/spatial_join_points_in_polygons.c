#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Define M_PI if not already defined by math.h
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- START MERSENNE TWISTER (Verbatim) ---
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
// --- END MERSENNE TWISTER ---

// --- Data Structures and Globals ---
typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    int start_index;
    int num_vertices;
} Polygon;

typedef struct {
    double angle;
    double radius;
} PolarCoord;

// Benchmark Parameters
int NUM_POINTS;
int NUM_POLYGONS;
int AVG_VERTICES_PER_POLYGON;

// Data
Point *points;
Polygon *polygons;
Point *all_vertices;
long long total_num_vertices = 0;

// Result
int final_result;

// --- Helper Functions ---
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

int compare_angles(const void *a, const void *b) {
    PolarCoord *p1 = (PolarCoord *)a;
    PolarCoord *p2 = (PolarCoord *)b;
    if (p1->angle < p2->angle) return -1;
    if (p1->angle > p2->angle) return 1;
    return 0;
}

// Ray-casting point-in-polygon test
int is_inside(const Point* point, const Polygon* polygon) {
    int i, j, c = 0;
    int nvert = polygon->num_vertices;
    Point* vertices = &all_vertices[polygon->start_index];
    double testx = point->x;
    double testy = point->y;

    for (i = 0, j = nvert - 1; i < nvert; j = i++) {
        if (((vertices[i].y > testy) != (vertices[j].y > testy)) &&
            (testx < (vertices[j].x - vertices[i].x) * (testy - vertices[i].y) / (vertices[j].y - vertices[i].y) + vertices[i].x)) {
            c = !c;
        }
    }
    return c;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_points num_polygons avg_vertices_per_polygon seed\n", argv[0]);
        exit(1);
    }

    NUM_POINTS = atoi(argv[1]);
    NUM_POLYGONS = atoi(argv[2]);
    AVG_VERTICES_PER_POLYGON = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    mt_seed(seed);

    // Allocate points array
    points = (Point *)malloc(NUM_POINTS * sizeof(Point));
    if (!points) { perror("Failed to allocate points"); exit(1); }

    // Allocate polygons array
    polygons = (Polygon *)malloc(NUM_POLYGONS * sizeof(Polygon));
    if (!polygons) { perror("Failed to allocate polygons"); exit(1); }

    // First pass: determine vertex counts and polygon start indices
    total_num_vertices = 0;
    for (int i = 0; i < NUM_POLYGONS; ++i) {
        int deviation = AVG_VERTICES_PER_POLYGON > 4 ? AVG_VERTICES_PER_POLYGON / 2 : 2;
        int vert_count = (AVG_VERTICES_PER_POLYGON - deviation) + (mt_rand() % (2 * deviation + 1));
        if (vert_count < 3) vert_count = 3;

        polygons[i].start_index = total_num_vertices;
        polygons[i].num_vertices = vert_count;
        total_num_vertices += vert_count;
    }

    // Allocate all vertices at once
    all_vertices = (Point *)malloc(total_num_vertices * sizeof(Point));
    if (!all_vertices) { perror("Failed to allocate all_vertices"); exit(1); }
    
    // Second pass: generate polygon vertices
    PolarCoord* polar_coords_temp = (PolarCoord*)malloc(AVG_VERTICES_PER_POLYGON * 2 * sizeof(PolarCoord));
    if(!polar_coords_temp) { perror("Failed to allocate polar coords temp"); exit(1); }

    for (int i = 0; i < NUM_POLYGONS; ++i) {
        double center_x = rand_double() * 1000.0;
        double center_y = rand_double() * 1000.0;
        double avg_radius = rand_double() * 40.0 + 10.0; 
        int num_v = polygons[i].num_vertices;

        for (int j = 0; j < num_v; ++j) {
            polar_coords_temp[j].angle = rand_double() * 2.0 * M_PI;
            polar_coords_temp[j].radius = avg_radius * (0.75 + rand_double() * 0.5);
        }

        qsort(polar_coords_temp, num_v, sizeof(PolarCoord), compare_angles);

        Point* current_vertices = &all_vertices[polygons[i].start_index];
        for (int j = 0; j < num_v; ++j) {
            current_vertices[j].x = center_x + polar_coords_temp[j].radius * cos(polar_coords_temp[j].angle);
            current_vertices[j].y = center_y + polar_coords_temp[j].radius * sin(polar_coords_temp[j].angle);
        }
    }
    free(polar_coords_temp);

    // Generate points
    for (int i = 0; i < NUM_POINTS; ++i) {
        points[i].x = rand_double() * 1000.0;
        points[i].y = rand_double() * 1000.0;
    }

    final_result = 0;
}

void run_computation() {
    int count = 0;
    for (int i = 0; i < NUM_POINTS; ++i) {
        for (int j = 0; j < NUM_POLYGONS; ++j) {
            if (is_inside(&points[i], &polygons[j])) {
                count++;
                break; // A point can be in only one polygon for this problem
            }
        }
    }
    final_result = count;
}

void cleanup() {
    free(points);
    free(polygons);
    free(all_vertices);
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

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
