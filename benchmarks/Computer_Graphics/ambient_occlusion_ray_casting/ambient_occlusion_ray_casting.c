#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Mersenne Twister (MT19937) generator - DO NOT MODIFY
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
// End of Mersenne Twister

// --- Utility Functions ---

// Generate a random float between 0.0 and 1.0
static inline float rand_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

// --- Data Structures ---

typedef struct {
    float x, y, z;
} vec3;

typedef struct {
    vec3 center;
    float radius;
    float radius_sq; // Precomputed for efficiency
} Sphere;

typedef struct {
    vec3 origin;
    vec3 dir;
} Ray;

// --- Vector Math ---

static inline vec3 vec3_sub(vec3 a, vec3 b) {
    return (vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline float vec3_dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// --- Benchmark Globals ---

// Parameters
static int image_width;
static int image_height;
static int num_samples_per_point;
static int num_objects;

// Scene data
static Sphere *objects;
static float *image_buffer;

// Final result
static float final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <image_width> <image_height> <num_samples_per_point> <num_objects> <seed>\n", argv[0]);
        exit(1);
    }

    image_width = atoi(argv[1]);
    image_height = atoi(argv[2]);
    num_samples_per_point = atoi(argv[3]);
    num_objects = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    objects = (Sphere *)malloc(num_objects * sizeof(Sphere));
    if (!objects) {
        fprintf(stderr, "Failed to allocate memory for objects.\n");
        exit(1);
    }
    
    image_buffer = (float *)malloc(image_width * image_height * sizeof(float));
    if (!image_buffer) {
        fprintf(stderr, "Failed to allocate memory for image buffer.\n");
        free(objects);
        exit(1);
    }
    
    // Generate random spheres within a unit cube [0,1] on x,y and [0.1, 1.1] on z
    for (int i = 0; i < num_objects; i++) {
        objects[i].center.x = rand_float();
        objects[i].center.y = rand_float();
        objects[i].center.z = rand_float() + 0.1f; // Place objects slightly above the plane
        objects[i].radius = 0.05f + rand_float() * 0.1f; // Small-ish spheres
        objects[i].radius_sq = objects[i].radius * objects[i].radius;
    }
}

// Check for ray-sphere intersection. Returns 1 if hit, 0 otherwise.
static inline int intersect_sphere(const Ray *ray, const Sphere *sphere) {
    vec3 oc = vec3_sub(ray->origin, sphere->center);
    float b = vec3_dot(oc, ray->dir);
    float c = vec3_dot(oc, oc) - sphere->radius_sq;
    
    // Optimization: if ray origin is outside sphere (c > 0) and 
    // ray is pointing away from sphere (b > 0), no hit is possible.
    if (c > 0.0f && b > 0.0f) return 0;

    // The discriminant of the quadratic equation t^2 + 2bt + c = 0.
    // We only need to know if there's a real root, so we check b^2 - c >= 0.
    float discriminant = b * b - c;
    
    return discriminant >= 0.0f;
}

void run_computation() {
    final_result = 0.0f;

    for (int y = 0; y < image_height; ++y) {
        for (int x = 0; x < image_width; ++x) {
            
            // The point on the xy-plane we're calculating occlusion for
            vec3 point_on_plane = {
                (float)x / (float)image_width,
                (float)y / (float)image_height,
                0.0f
            };

            int occluded_samples = 0;
            for (int s = 0; s < num_samples_per_point; ++s) {
                // Generate a random direction in the upper hemisphere (cosine-weighted)
                float r1 = rand_float(); // for z component (cos_theta)
                float r2 = rand_float(); // for angle phi
                
                float sin_theta = sqrtf(1.0f - r1 * r1);
                float phi = 2.0f * M_PI * r2;
                
                vec3 rand_dir = {
                    cosf(phi) * sin_theta,
                    sinf(phi) * sin_theta,
                    r1 // z component
                };

                Ray occlusion_ray = {point_on_plane, rand_dir};

                // Check for intersection with any object
                for (int i = 0; i < num_objects; ++i) {
                    if (intersect_sphere(&occlusion_ray, &objects[i])) {
                        occluded_samples++;
                        break; // One hit is enough for this sample ray
                    }
                }
            }

            // AO factor is the proportion of non-occluded rays
            float ao_factor = 1.0f - ((float)occluded_samples / (float)num_samples_per_point);
            image_buffer[y * image_width + x] = ao_factor;
        }
    }

    // Accumulate the result to prevent dead code elimination
    for (int i = 0; i < image_width * image_height; ++i) {
        final_result += image_buffer[i];
    }
}

void cleanup() {
    free(objects);
    free(image_buffer);
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
