/*
 * BENCHMARK: Geospatial Information Systems (GIS) Analysis
 * PROGRAM: dissolve_polygons_by_attribute
 * DESCRIPTION: This benchmark simulates a common GIS operation called "dissolve".
 *              It merges adjacent polygons that share a common attribute value.
 *              The simulation involves iterating through pairs of polygons, checking
 *              for a matching attribute, and then performing a spatial check (bounding
 *              box overlap) to determine adjacency. If conditions are met, one
 *              polygon is marked as "dissolved". The final result is the total number
 *              of vertices in the remaining, non-dissolved polygons.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// MERSENNE TWISTER (DO NOT MODIFY) 
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
// END MERSENNE TWISTER

// --- BENCHMARK DATA STRUCTURES ---
typedef struct {
    double x;
    double y;
} Vertex;

typedef struct {
    double min_x, max_x, min_y, max_y;
} BoundingBox;

typedef struct {
    Vertex* vertices;
    int num_vertices;
    int attribute_id;
    int is_dissolved;
    BoundingBox bbox; // Pre-calculated for performance
} Polygon;

// --- GLOBAL STATE ---
static Polygon* g_polygons;
static int g_num_input_polygons;
static long long g_final_result; // Use long long for large vertex counts

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_polygons> <avg_vertices> <num_attributes> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_input_polygons = atoi(argv[1]);
    int avg_vertices_per_polygon = atoi(argv[2]);
    int num_attributes = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (g_num_input_polygons <= 0 || avg_vertices_per_polygon <= 0 || num_attributes <= 0) {
        fprintf(stderr, "FATAL: All numerical arguments must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    g_polygons = (Polygon*)malloc(g_num_input_polygons * sizeof(Polygon));
    if (!g_polygons) {
        fprintf(stderr, "FATAL: Failed to allocate memory for polygons\n");
        exit(1);
    }

    const double world_size = 10000.0;

    for (int i = 0; i < g_num_input_polygons; i++) {
        // Randomize number of vertices around the average, ensuring at least 3.
        int num_verts = (avg_vertices_per_polygon / 2) + (mt_rand() % avg_vertices_per_polygon);
        if (num_verts < 3) num_verts = 3;

        g_polygons[i].num_vertices = num_verts;
        g_polygons[i].vertices = (Vertex*)malloc(num_verts * sizeof(Vertex));
        if (!g_polygons[i].vertices) {
            fprintf(stderr, "FATAL: Failed to allocate memory for vertices\n");
            exit(1);
        }

        // Generate a somewhat realistic 'convex-like' random polygon
        double center_x = (double)mt_rand() / UINT32_MAX * world_size;
        double center_y = (double)mt_rand() / UINT32_MAX * world_size;
        double radius = ((double)mt_rand() / UINT32_MAX) * 100.0 + 50.0;

        // Initialize bounding box for calculation
        BoundingBox* bbox = &g_polygons[i].bbox;
        bbox->min_x = world_size * 2; bbox->max_x = -world_size * 2;
        bbox->min_y = world_size * 2; bbox->max_y = -world_size * 2;

        for (int j = 0; j < num_verts; j++) {
            double angle = 2.0 * M_PI * j / num_verts;
            double jitter = ((double)mt_rand() / UINT32_MAX - 0.5) * radius * 0.5;
            double current_radius = radius + jitter;

            double x = center_x + current_radius * cos(angle);
            double y = center_y + current_radius * sin(angle);

            g_polygons[i].vertices[j].x = x;
            g_polygons[i].vertices[j].y = y;

            // Update bounding box on-the-fly
            if (x < bbox->min_x) bbox->min_x = x;
            if (x > bbox->max_x) bbox->max_x = x;
            if (y < bbox->min_y) bbox->min_y = y;
            if (y > bbox->max_y) bbox->max_y = y;
        }

        g_polygons[i].attribute_id = mt_rand() % num_attributes;
        g_polygons[i].is_dissolved = 0;
    }
}

void cleanup() {
    if (g_polygons) {
        for (int i = 0; i < g_num_input_polygons; i++) {
            free(g_polygons[i].vertices);
        }
        free(g_polygons);
    }
}

// Fast check for bounding box overlap (spatial adjacency proxy)
static inline int check_bbox_overlap(const BoundingBox* b1, const BoundingBox* b2) {
    if (b1->max_x < b2->min_x || b1->min_x > b2->max_x ||
        b1->max_y < b2->min_y || b1->min_y > b2->max_y) {
        return 0; // No overlap
    }
    return 1; // Overlap
}

void run_computation() {
    for (int i = 0; i < g_num_input_polygons; ++i) {
        if (g_polygons[i].is_dissolved) {
            continue;
        }

        for (int j = i + 1; j < g_num_input_polygons; ++j) {
            if (g_polygons[j].is_dissolved) {
                continue;
            }
            
            // Main condition for dissolving: same attribute and spatial adjacency
            if (g_polygons[i].attribute_id == g_polygons[j].attribute_id) {
                if (check_bbox_overlap(&g_polygons[i].bbox, &g_polygons[j].bbox)) {
                    // A real GIS would perform a geometric union here. We simulate
                    // this by simply marking one polygon as dissolved.
                    g_polygons[j].is_dissolved = 1;
                }
            }
        }
    }

    // Calculate a final result to prevent dead code elimination.
    long long total_vertices = 0;
    for (int i = 0; i < g_num_input_polygons; ++i) {
        if (!g_polygons[i].is_dissolved) {
            total_vertices += g_polygons[i].num_vertices;
        }
    }
    g_final_result = total_vertices;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    double time_taken;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
