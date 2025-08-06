#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- Data Structures ---
typedef struct {
    float x, y;
} Point;

typedef struct {
    Point* vertices;
    int num_vertices;
} Polygon;

// --- Global Benchmark Data ---
static int num_polygons;
static int avg_vertices_per_polygon;
static int output_grid_width;
static int output_grid_height;

static Polygon* polygons;
static int* output_grid;
static int final_result;


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_polygons avg_vertices_per_polygon output_grid_width output_grid_height seed\n", argv[0]);
        exit(1);
    }

    num_polygons = atoi(argv[1]);
    avg_vertices_per_polygon = atoi(argv[2]);
    output_grid_width = atoi(argv[3]);
    output_grid_height = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    polygons = (Polygon*)malloc(num_polygons * sizeof(Polygon));
    if (!polygons) {
        perror("Failed to allocate polygons array");
        exit(1);
    }

    for (int i = 0; i < num_polygons; ++i) {
        int n_verts = (avg_vertices_per_polygon / 2) + (mt_rand() % avg_vertices_per_polygon);
        if (n_verts < 3) n_verts = 3;

        polygons[i].num_vertices = n_verts;
        polygons[i].vertices = (Point*)malloc(n_verts * sizeof(Point));
        if (!polygons[i].vertices) {
            perror("Failed to allocate vertex array");
            exit(1);
        }

        float cx = (mt_rand() / (float)UINT32_MAX) * (output_grid_width * 0.8f) + (output_grid_width * 0.1f);
        float cy = (mt_rand() / (float)UINT32_MAX) * (output_grid_height * 0.8f) + (output_grid_height * 0.1f);
        float radius = (mt_rand() / (float)UINT32_MAX) * 50.0f + 20.0f;

        for (int j = 0; j < n_verts; ++j) {
            float angle = 2.0f * M_PI * j / n_verts;
            float r = radius * (0.8f + (mt_rand() / (float)UINT32_MAX) * 0.4f);
            polygons[i].vertices[j].x = cx + r * cosf(angle);
            polygons[i].vertices[j].y = cy + r * sinf(angle);
        }
    }

    size_t grid_size_bytes = (size_t)output_grid_width * output_grid_height * sizeof(int);
    output_grid = (int*)malloc(grid_size_bytes);
    if (!output_grid) {
        perror("Failed to allocate output grid");
        exit(1);
    }
    memset(output_grid, 0, grid_size_bytes);
}

static int compare_floats(const void* a, const void* b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

void run_computation() {
    float* intersections = (float*)malloc(avg_vertices_per_polygon * 2 * sizeof(float));
    if (!intersections) {
        perror("Failed to allocate intersection buffer");
        exit(1);
    }

    for (int i = 0; i < num_polygons; ++i) {
        Polygon* p = &polygons[i];
        int n = p->num_vertices;
        Point* v = p->vertices;

        float min_y = v[0].y, max_y = v[0].y;
        for (int j = 1; j < n; ++j) {
            if (v[j].y < min_y) min_y = v[j].y;
            if (v[j].y > max_y) max_y = v[j].y;
        }

        int start_y = (int)ceilf(min_y);
        int end_y = (int)floorf(max_y);
        if (start_y < 0) start_y = 0;
        if (end_y >= output_grid_height) end_y = output_grid_height - 1;

        for (int y = start_y; y <= end_y; ++y) {
            int num_intersections = 0;
            for (int j = 0; j < n; ++j) {
                Point p1 = v[j];
                Point p2 = v[(j + 1) % n];

                if ((p1.y <= y && p2.y > y) || (p2.y <= y && p1.y > y)) {
                    if (num_intersections < avg_vertices_per_polygon * 2) {
                         intersections[num_intersections++] = (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y) + p1.x;
                    }
                }
            }

            qsort(intersections, num_intersections, sizeof(float), compare_floats);

            for (int k = 0; k < num_intersections; k += 2) {
                if (k + 1 < num_intersections) {
                    int x_start = (int)ceilf(intersections[k]);
                    int x_end = (int)floorf(intersections[k + 1]);

                    if (x_start < 0) x_start = 0;
                    if (x_end >= output_grid_width) x_end = output_grid_width - 1;

                    for (int x = x_start; x <= x_end; ++x) {
                         if (y >= 0 && y < output_grid_height && x >= 0 && x < output_grid_width) {
                             output_grid[y * output_grid_width + x]++;
                         }
                    }
                }
            }
        }
    }

    free(intersections);

    int checksum = 0;
    size_t grid_size = (size_t)output_grid_width * output_grid_height;
    for (size_t i = 0; i < grid_size; ++i) {
        checksum = (checksum + output_grid[i] * (int)((i % 128) + 1)) & 0x7FFFFFFF;
    }
    final_result = checksum;
}

void cleanup() {
    if (polygons) {
        for (int i = 0; i < num_polygons; ++i) {
            free(polygons[i].vertices);
        }
        free(polygons);
    }
    free(output_grid);
}

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
