/*
 * Copyright (c) 2024, The Benchkit Project
 * All rights reserved. 
 * 
 * The Vertex Shader Skinning benchmark simulates a common GPU task on the CPU.
 * It calculates the final position of mesh vertices based on the influence of a
 * set of animated "bones" (a skeletal animation technique).
 * Each vertex is transformed by multiple bone matrices, and the results are
 * blended together using weighting factors.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>


// ==========================================================================
// Mersenne Twister (MT19937) PRNG (DO NOT MODIFY)
// ==========================================================================
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

// ==========================================================================
// Benchmark-specific data structures
// ==========================================================================

#define BONES_PER_VERTEX 4

typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    float m[16]; // Column-major order: m[0], m[1], m[2], m[3] is the first column
} Mat4;

// Global parameters
static int num_vertices;
static int num_bones;

// Global data arrays
static Vec3 *base_vertices = NULL;
static Vec3 *skinned_vertices = NULL;
static Mat4 *bone_matrices = NULL;
static float *vertex_weights = NULL;
static int *bone_indices = NULL;

// Global variable to hold the final result
static volatile float final_result;

// ==========================================================================
// Helper functions
// ==========================================================================

// Generates a random float in [0.0, 1.0)
float rand_float() {
    return (float)mt_rand() / ((float)UINT32_MAX + 1.0f);
}

// ==========================================================================
// Benchmark Setup, Computation, and Cleanup
// ==========================================================================

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_bones> <seed>\n", argv[0]);
        exit(1);
    }

    num_vertices = atoi(argv[1]);
    num_bones = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_vertices <= 0 || num_bones <= 0) {
        fprintf(stderr, "Error: num_vertices and num_bones must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory from the heap
    base_vertices = (Vec3 *)malloc(num_vertices * sizeof(Vec3));
    skinned_vertices = (Vec3 *)malloc(num_vertices * sizeof(Vec3));
    bone_matrices = (Mat4 *)malloc(num_bones * sizeof(Mat4));
    vertex_weights = (float *)malloc(num_vertices * BONES_PER_VERTEX * sizeof(float));
    bone_indices = (int *)malloc(num_vertices * BONES_PER_VERTEX * sizeof(int));

    if (!base_vertices || !skinned_vertices || !bone_matrices || !vertex_weights || !bone_indices) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize base vertices with random positions
    for (int i = 0; i < num_vertices; i++) {
        base_vertices[i] = (Vec3){rand_float() * 2.0f - 1.0f, rand_float() * 2.0f - 1.0f, rand_float() * 2.0f - 1.0f};
    }

    // Initialize vertex weights and bone indices
    for (int i = 0; i < num_vertices; i++) {
        float total_weight = 0.0f;
        int influence_offset = i * BONES_PER_VERTEX;
        for (int j = 0; j < BONES_PER_VERTEX; j++) {
            float weight = rand_float();
            vertex_weights[influence_offset + j] = weight;
            total_weight += weight;
            bone_indices[influence_offset + j] = mt_rand() % num_bones;
        }

        // Normalize weights so they sum to 1.0
        for (int j = 0; j < BONES_PER_VERTEX; j++) {
            vertex_weights[influence_offset + j] /= total_weight;
        }
    }

    // Initialize bone matrices (simple random translation)
    for (int i = 0; i < num_bones; i++) {
        memset(bone_matrices[i].m, 0, sizeof(Mat4));
        bone_matrices[i].m[0] = 1.0f;  // [0][0]
        bone_matrices[i].m[5] = 1.0f;  // [1][1]
        bone_matrices[i].m[10] = 1.0f; // [2][2]
        bone_matrices[i].m[15] = 1.0f; // [3][3]

        // Translation components in column 3
        bone_matrices[i].m[12] = rand_float() * 10.0f - 5.0f; // tx
        bone_matrices[i].m[13] = rand_float() * 10.0f - 5.0f; // ty
        bone_matrices[i].m[14] = rand_float() * 10.0f - 5.0f; // tz
    }
}

void run_computation() {
    float accumulator = 0.0f;
    
    for (int i = 0; i < num_vertices; ++i) {
        const Vec3 *src_pos = &base_vertices[i];
        
        float final_x = 0.0f;
        float final_y = 0.0f;
        float final_z = 0.0f;

        const int influence_offset = i * BONES_PER_VERTEX;
        for (int j = 0; j < BONES_PER_VERTEX; ++j) {
            const float weight = vertex_weights[influence_offset + j];
            const int bone_idx = bone_indices[influence_offset + j];
            const Mat4 *bone_matrix = &bone_matrices[bone_idx];

            // Transform vertex by bone matrix (assuming w=1 for affine transform)
            // Using column-major matrix layout
            const float tx = bone_matrix->m[0] * src_pos->x + bone_matrix->m[4] * src_pos->y + bone_matrix->m[8]  * src_pos->z + bone_matrix->m[12];
            const float ty = bone_matrix->m[1] * src_pos->x + bone_matrix->m[5] * src_pos->y + bone_matrix->m[9]  * src_pos->z + bone_matrix->m[13];
            const float tz = bone_matrix->m[2] * src_pos->x + bone_matrix->m[6] * src_pos->y + bone_matrix->m[10] * src_pos->z + bone_matrix->m[14];
            
            // Accumulate weighted transformed position
            final_x += tx * weight;
            final_y += ty * weight;
            final_z += tz * weight;
        }

        skinned_vertices[i] = (Vec3){final_x, final_y, final_z};
    }

    // Sum a component of the results to prevent dead code elimination.
    // This mimics a subsequent pass that would use the skinned data.
    for (int i = 0; i < num_vertices; ++i) {
        accumulator += skinned_vertices[i].x;
    }
    final_result = accumulator;
}

void cleanup() {
    free(base_vertices);
    free(skinned_vertices);
    free(bone_matrices);
    free(vertex_weights);
    free(bone_indices);
    base_vertices = NULL;
    skinned_vertices = NULL;
    bone_matrices = NULL;
    vertex_weights = NULL;
    bone_indices = NULL;
}

// ==========================================================================
// Main function
// ==========================================================================

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final computational result to stdout for verification
    printf("%f\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
