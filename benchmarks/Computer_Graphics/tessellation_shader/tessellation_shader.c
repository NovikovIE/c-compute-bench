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
// --- End Mersenne Twister ---

// --- Benchmark Data Structures and Globals ---
typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    Vec3 cp[4]; // Control points for a quad patch (P00, P10, P01, P11)
} Patch;

// Parameters
static int num_patches;
static int tessellation_level;

// Data buffers
static Patch *patches;
static Vec3 *output_vertices;

// Final result
double final_result;

// --- Benchmark Functions ---

// Helper to generate a random float between -1.0 and 1.0
float random_float() {
    return ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_patches> <tessellation_level> <seed>\n", argv[0]);
        exit(1);
    }

    num_patches = atoi(argv[1]);
    tessellation_level = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    patches = (Patch *)malloc(num_patches * sizeof(Patch));
    if (!patches) {
        fprintf(stderr, "Failed to allocate memory for patches.\n");
        exit(1);
    }

    for (int i = 0; i < num_patches; ++i) {
        for (int j = 0; j < 4; ++j) {
            patches[i].cp[j] = (Vec3){random_float(), random_float(), random_float()};
        }
    }

    size_t verts_per_patch = (size_t)(tessellation_level + 1) * (tessellation_level + 1);
    size_t total_verts = (size_t)num_patches * verts_per_patch;
    output_vertices = (Vec3 *)malloc(total_verts * sizeof(Vec3));
    if (!output_vertices) {
        fprintf(stderr, "Failed to allocate memory for output vertices.\n");
        free(patches);
        exit(1);
    }

    final_result = 0.0;
}

// Linear interpolation for vectors
Vec3 lerp(Vec3 a, Vec3 b, float t) {
    Vec3 result;
    result.x = a.x + t * (b.x - a.x);
    result.y = a.y + t * (b.y - a.y);
    result.z = a.z + t * (b.z - a.z);
    return result;
}

void run_computation() {
    double accumulator = 0.0;
    size_t verts_per_patch = (size_t)(tessellation_level + 1) * (tessellation_level + 1);
    float tess_level_f = (float)tessellation_level;

    for (int i = 0; i < num_patches; ++i) {
        const Vec3 p00 = patches[i].cp[0];
        const Vec3 p10 = patches[i].cp[1];
        const Vec3 p01 = patches[i].cp[2];
        const Vec3 p11 = patches[i].cp[3];

        for (int v_step = 0; v_step <= tessellation_level; ++v_step) {
            for (int u_step = 0; u_step <= tessellation_level; ++u_step) {
                float u = (float)u_step / tess_level_f;
                float v = (float)v_step / tess_level_f;

                // Bilinear interpolation to find the base vertex position
                Vec3 temp1 = lerp(p00, p10, u);
                Vec3 temp2 = lerp(p01, p11, u);
                Vec3 new_vertex = lerp(temp1, temp2, v);

                // Apply a procedural displacement (simulating a displacement shader)
                new_vertex.z += 0.1f * sinf(new_vertex.x * 10.0f) * cosf(new_vertex.y * 10.0f);

                // Store the result and accumulate a value to prevent dead code elimination
                size_t index = (size_t)i * verts_per_patch + (size_t)v_step * (tessellation_level + 1) + u_step;
                output_vertices[index] = new_vertex;
                accumulator += new_vertex.z;
            }
        }
    }
    final_result = accumulator;
}

void cleanup() {
    free(patches);
    free(output_vertices);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
