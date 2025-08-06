#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- DO NOT MODIFY ---
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
// --- END DO NOT MODIFY ---

// Benchmark parameters
int num_features;
int avg_vertices_per_feature;

// Data structures
typedef struct {
    double x;
    double y;
} Vertex;

typedef struct {
    Vertex* vertices;
    int num_vertices;
} Feature;

// Global data pointers
Feature* source_layer = NULL;
Feature* reprojected_layer = NULL;
double final_checksum = 0.0;

// Helper function to generate a random double
double rand_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_features> <avg_vertices_per_feature> <seed>\n", argv[0]);
        exit(1);
    }

    num_features = atoi(argv[1]);
    avg_vertices_per_feature = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (num_features <= 0 || avg_vertices_per_feature <= 0) {
        fprintf(stderr, "Parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    source_layer = (Feature*)malloc(num_features * sizeof(Feature));
    reprojected_layer = (Feature*)malloc(num_features * sizeof(Feature));

    if (!source_layer || !reprojected_layer) {
        fprintf(stderr, "Failed to allocate memory for layers.\n");
        exit(1);
    }
    
    for (int i = 0; i < num_features; ++i) {
        int num_vertices = (avg_vertices_per_feature / 2) + (mt_rand() % avg_vertices_per_feature);
        if (num_vertices < 3) num_vertices = 3;

        source_layer[i].num_vertices = num_vertices;
        reprojected_layer[i].num_vertices = num_vertices;

        source_layer[i].vertices = (Vertex*)malloc(num_vertices * sizeof(Vertex));
        reprojected_layer[i].vertices = (Vertex*)malloc(num_vertices * sizeof(Vertex));

        if (!source_layer[i].vertices || !reprojected_layer[i].vertices) {
            fprintf(stderr, "Failed to allocate memory for vertices.\n");
            exit(1);
        }

        for (int j = 0; j < num_vertices; ++j) {
            // Generate coordinates roughly simulating longitude/latitude
            source_layer[i].vertices[j].x = rand_double(-180.0, 180.0);
            source_layer[i].vertices[j].y = rand_double(-90.0, 90.0);
        }
    }
}

void run_computation() {
    // A simple affine transformation matrix for reprojection
    const double T[2][3] = {
        {1.0001, 0.0002, 500.0},  // {a, b, tx}
        {-0.0002, 0.9999, -250.0}  // {c, d, ty}
    };

    double checksum = 0.0;

    for (int i = 0; i < num_features; ++i) {
        for (int j = 0; j < source_layer[i].num_vertices; ++j) {
            double old_x = source_layer[i].vertices[j].x;
            double old_y = source_layer[i].vertices[j].y;

            double new_x = T[0][0] * old_x + T[0][1] * old_y + T[0][2];
            double new_y = T[1][0] * old_x + T[1][1] * old_y + T[1][2];

            reprojected_layer[i].vertices[j].x = new_x;
            reprojected_layer[i].vertices[j].y = new_y;

            checksum += new_x + new_y;
        }
    }
    final_checksum = checksum;
}

void cleanup() {
    if (source_layer) {
        for (int i = 0; i < num_features; ++i) {
            free(source_layer[i].vertices);
        }
        free(source_layer);
    }
    if (reprojected_layer) {
        for (int i = 0; i < num_features; ++i) {
            free(reprojected_layer[i].vertices);
        }
        free(reprojected_layer);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%.4f\n", final_checksum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
