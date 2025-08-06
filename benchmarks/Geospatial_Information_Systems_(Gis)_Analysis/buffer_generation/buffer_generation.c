#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// M_PI is not in C standard, but defined in POSIX and common
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

// --- Data Structures ---
typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    Point* vertices;
    int num_vertices;
} Feature;

// --- Global Benchmark State ---
static int g_num_features;
static int g_buffer_resolution_segments;
static double g_buffer_distance;

static Feature* g_input_features;
static Feature* g_output_buffers;

static double g_final_result = 0.0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_features avg_vertices_per_feature buffer_resolution_segments buffer_distance seed\n", argv[0]);
        exit(1);
    }

    g_num_features = atoi(argv[1]);
    int avg_vertices_per_feature = atoi(argv[2]);
    g_buffer_resolution_segments = atoi(argv[3]);
    g_buffer_distance = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    g_input_features = (Feature*)malloc(g_num_features * sizeof(Feature));
    g_output_buffers = (Feature*)malloc(g_num_features * sizeof(Feature));
    if (!g_input_features || !g_output_buffers) {
        fprintf(stderr, "Failed to allocate memory for features.\n");
        exit(1);
    }

    double max_coord = 10000.0;

    for (int i = 0; i < g_num_features; ++i) {
        // Generate a random number of vertices for this feature, ensuring at least a line segment (2 vertices)
        int num_vertices = 2;
        if (avg_vertices_per_feature > 2) {
             num_vertices += (mt_rand() % (avg_vertices_per_feature * 2 - 4));
        }

        g_input_features[i].num_vertices = num_vertices;
        g_input_features[i].vertices = (Point*)malloc(num_vertices * sizeof(Point));
        if (!g_input_features[i].vertices) {
            fprintf(stderr, "Failed to allocate memory for vertices.\n");
            exit(1);
        }

        for (int j = 0; j < num_vertices; ++j) {
            g_input_features[i].vertices[j].x = ((double)mt_rand() / (double)UINT32_MAX) * max_coord;
            g_input_features[i].vertices[j].y = ((double)mt_rand() / (double)UINT32_MAX) * max_coord;
        }

        // Allocate memory for the output buffer polygon
        // The simplified model generates a buffer for each vertex
        int buffer_vertex_count = num_vertices * g_buffer_resolution_segments;
        g_output_buffers[i].num_vertices = buffer_vertex_count;
        g_output_buffers[i].vertices = (Point*)malloc(buffer_vertex_count * sizeof(Point));
        if (!g_output_buffers[i].vertices) {
            fprintf(stderr, "Failed to allocate memory for buffer vertices.\n");
            exit(1);
        }
    }
}

void run_computation() {
    double total_x_sum = 0.0;
    const double two_pi = 2.0 * M_PI;

    for (int i = 0; i < g_num_features; ++i) {
        Feature* in_feat = &g_input_features[i];
        Feature* out_buff = &g_output_buffers[i];
        int buffer_vertex_idx = 0;

        for (int j = 0; j < in_feat->num_vertices; ++j) {
            Point center = in_feat->vertices[j];

            for (int k = 0; k < g_buffer_resolution_segments; ++k) {
                double angle = two_pi * (double)k / (double)g_buffer_resolution_segments;
                double buffer_x = center.x + g_buffer_distance * cos(angle);
                double buffer_y = center.y + g_buffer_distance * sin(angle);
                
                out_buff->vertices[buffer_vertex_idx].x = buffer_x;
                out_buff->vertices[buffer_vertex_idx].y = buffer_y;
                buffer_vertex_idx++;

                total_x_sum += buffer_x; // Accumulate to prevent dead-code elimination
            }
        }
    }
    g_final_result = total_x_sum;
}

void cleanup() {
    for (int i = 0; i < g_num_features; ++i) {
        free(g_input_features[i].vertices);
        free(g_output_buffers[i].vertices);
    }
    free(g_input_features);
    free(g_output_buffers);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
