#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

// --- Mersenne Twister (MT19937) Generator ---
// Do Not Modify
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
// --- End of Mersenne Twister ---

// --- Benchmark Configuration ---

// Benchmark-specific data structures
typedef struct { float x, y, z; } Point3D;

// Parameters from command line
static int NUM_POINTS_CLOUD_A;
static int NUM_POINTS_CLOUD_B;
static int MAX_ITERATIONS;

// Global data pointers
static Point3D *cloud_a; // Source cloud
static Point3D *cloud_b; // Target cloud
static Point3D *transformed_a; // Source cloud being transformed during iterations

// Final result to prevent dead code elimination
float final_result_sum;

// Helper to generate a random float in [0, 1)
float rand_float_01() {
    return (float)mt_rand() / 4294967296.0f; // 2^32
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    // 1. Argument parsing
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_points_a> <num_points_b> <max_iterations> <seed>\n", argv[0]);
        exit(1);
    }
    NUM_POINTS_CLOUD_A = atoi(argv[1]);
    NUM_POINTS_CLOUD_B = atoi(argv[2]);
    MAX_ITERATIONS = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    // 2. Seed the random number generator
    mt_seed(seed);

    // 3. Allocate memory
    cloud_a = (Point3D*) malloc(NUM_POINTS_CLOUD_A * sizeof(Point3D));
    cloud_b = (Point3D*) malloc(NUM_POINTS_CLOUD_B * sizeof(Point3D));
    transformed_a = (Point3D*) malloc(NUM_POINTS_CLOUD_A * sizeof(Point3D));

    if (!cloud_a || !cloud_b || !transformed_a) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // 4. Generate data
    // Generate cloud B as a random set of points
    for (int i = 0; i < NUM_POINTS_CLOUD_B; i++) {
        cloud_b[i].x = rand_float_01() * 10.0f;
        cloud_b[i].y = rand_float_01() * 10.0f;
        cloud_b[i].z = rand_float_01() * 10.0f;
    }

    // Generate cloud A as another random set of points, shifted slightly
    // A real ICP test would use a transformed copy of B, but for a benchmark
    // of computational load, two distinct clouds are sufficient.
    for (int i = 0; i < NUM_POINTS_CLOUD_A; i++) {
        cloud_a[i].x = rand_float_01() * 10.0f + 0.5f; // Small offset
        cloud_a[i].y = rand_float_01() * 10.0f - 0.3f;
        cloud_a[i].z = rand_float_01() * 10.0f + 0.1f;
        transformed_a[i] = cloud_a[i]; // Initialize transformed cloud
    }
}

void run_computation() {
    float total_dist_sq_sum = 0.0f;

    // A small, fixed rotation and translation to be applied each iteration
    // to simulate the transformation update step of ICP. This is not a
    // true ICP update but provides a consistent computational load.
    const float theta = 0.01f; // ~0.57 degrees
    const float cos_t = cosf(theta);
    const float sin_t = sinf(theta);
    
    // Rotation around Z axis and a small translation
    float rot[3][3] = {
        {cos_t, -sin_t, 0.0f},
        {sin_t,  cos_t, 0.0f},
        {0.0f,   0.0f,  1.0f}
    };
    float trans[3] = {0.01f, -0.01f, 0.005f};

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        // Step 1: Find closest point correspondences (the most expensive part)
        // This is a naive O(N*M) search.
        for (int i = 0; i < NUM_POINTS_CLOUD_A; ++i) {
            float min_dist_sq = FLT_MAX;
            Point3D p_a = transformed_a[i];
            
            for (int j = 0; j < NUM_POINTS_CLOUD_B; ++j) {
                Point3D p_b = cloud_b[j];
                float dist_sq = (p_a.x - p_b.x) * (p_a.x - p_b.x) +
                                (p_a.y - p_b.y) * (p_a.y - p_b.y) +
                                (p_a.z - p_b.z) * (p_a.z - p_b.z);
                if (dist_sq < min_dist_sq) {
                    min_dist_sq = dist_sq;
                }
            }
            // We don't need the correspondences for this benchmark, just the computation.
        }

        // Step 2: Apply the fixed transformation to the source cloud
        for (int i = 0; i < NUM_POINTS_CLOUD_A; ++i) {
            Point3D p = transformed_a[i];
            transformed_a[i].x = rot[0][0] * p.x + rot[0][1] * p.y + rot[0][2] * p.z + trans[0];
            transformed_a[i].y = rot[1][0] * p.x + rot[1][1] * p.y + rot[1][2] * p.z + trans[1];
            transformed_a[i].z = rot[2][0] * p.x + rot[2][1] * p.y + rot[2][2] * p.z + trans[2];
        }
    }

    // Final calculation to generate a result based on the final transformed cloud
    for (int i = 0; i < NUM_POINTS_CLOUD_A; ++i) {
        float min_dist_sq = FLT_MAX;
        Point3D p_a = transformed_a[i];
        
        for (int j = 0; j < NUM_POINTS_CLOUD_B; ++j) {
            Point3D p_b = cloud_b[j];
            float dist_sq = (p_a.x - p_b.x) * (p_a.x - p_b.x) +
                            (p_a.y - p_b.y) * (p_a.y - p_b.y) +
                            (p_a.z - p_b.z) * (p_a.z - p_b.z);
            if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
            }
        }
        total_dist_sq_sum += min_dist_sq;
    }
    
    final_result_sum = total_dist_sq_sum;
}

void cleanup() {
    free(cloud_a);
    free(cloud_b);
    free(transformed_a);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print final result to stdout
    printf("%f\n", final_result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
