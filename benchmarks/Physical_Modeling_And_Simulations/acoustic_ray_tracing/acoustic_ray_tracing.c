#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

// --- Utility Functions and Structs ---

#define RAND_DOUBLE() ((double)mt_rand() / (double)UINT32_MAX)
const double EPSILON = 1e-6;

typedef struct {
    double x, y, z;
} Vector3D;

typedef struct {
    Vector3D origin;
    Vector3D direction;
    double energy;
} Ray;

typedef struct {
    Vector3D normal;
    double d; // Distance from origin for plane equation normal.P + d = 0
    double absorption_coeff;
} Plane;

// --- Global Data Structure ---

typedef struct {
    int num_rays;
    int num_reflections;
    int num_planes; // from room_complexity
    Ray* rays;
    Plane* room_planes;
    double total_absorbed_energy;
} BenchmarkData;

static BenchmarkData g_data;

// --- Vector Math --- 

static inline double vec_dot(Vector3D a, Vector3D b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vector3D vec_scale(Vector3D v, double s) {
    return (Vector3D){v.x * s, v.y * s, v.z * s};
}

static inline Vector3D vec_sub(Vector3D a, Vector3D b) {
    return (Vector3D){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vector3D vec_add(Vector3D a, Vector3D b) {
    return (Vector3D){a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vector3D vec_normalize(Vector3D v) {
    double mag = sqrt(vec_dot(v, v));
    if (mag > EPSILON) {
        return vec_scale(v, 1.0 / mag);
    }
    return (Vector3D){0, 0, 0};
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_rays num_reflections room_complexity seed\n", argv[0]);
        exit(1);
    }

    g_data.num_rays = atoi(argv[1]);
    g_data.num_reflections = atoi(argv[2]);
    g_data.num_planes = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    mt_seed(seed);

    g_data.rays = (Ray*)malloc(g_data.num_rays * sizeof(Ray));
    g_data.room_planes = (Plane*)malloc(g_data.num_planes * sizeof(Plane));
    if (!g_data.rays || !g_data.room_planes) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize room planes
    for (int i = 0; i < g_data.num_planes; i++) {
        Vector3D normal;
        do {
            normal.x = 2.0 * RAND_DOUBLE() - 1.0;
            normal.y = 2.0 * RAND_DOUBLE() - 1.0;
            normal.z = 2.0 * RAND_DOUBLE() - 1.0;
        } while (vec_dot(normal, normal) < EPSILON); // Ensure non-zero vector
        g_data.room_planes[i].normal = vec_normalize(normal);
        // Place planes within a bounding sphere of radius ~10
        g_data.room_planes[i].d = (2.0 * RAND_DOUBLE() - 1.0) * 10.0; 
        g_data.room_planes[i].absorption_coeff = 0.05 + 0.2 * RAND_DOUBLE(); // 5% to 25% absorption
    }

    // Initialize rays
    Vector3D initial_origin = {0.0, 0.0, 0.0};
    for (int i = 0; i < g_data.num_rays; i++) {
        g_data.rays[i].origin = initial_origin;
        Vector3D dir;
        double mag_sq;
        do {
          dir.x = 2.0 * RAND_DOUBLE() - 1.0;
          dir.y = 2.0 * RAND_DOUBLE() - 1.0;
          dir.z = 2.0 * RAND_DOUBLE() - 1.0;
          mag_sq = vec_dot(dir, dir);
        } while (mag_sq > 1.0 || mag_sq < EPSILON); // Rejection sampling for uniform sphere distribution
        g_data.rays[i].direction = vec_normalize(dir);
        g_data.rays[i].energy = 1.0;
    }
    
    g_data.total_absorbed_energy = 0.0;
}

void run_computation() {
    double total_absorbed = 0.0;

    for (int i = 0; i < g_data.num_rays; i++) {
        Ray current_ray = g_data.rays[i];

        for (int j = 0; j < g_data.num_reflections; j++) {
            double min_t = 1e30; // effectively infinity
            int hit_plane_idx = -1;

            // Find the closest intersection
            for (int k = 0; k < g_data.num_planes; k++) {
                double denom = vec_dot(g_data.room_planes[k].normal, current_ray.direction);
                
                // Only consider intersections with planes we're facing
                if (fabs(denom) > EPSILON) {
                    double numer = -(vec_dot(g_data.room_planes[k].normal, current_ray.origin) + g_data.room_planes[k].d);
                    double t = numer / denom;
                    
                    // Intersection must be in front of the ray and closer than previous hits
                    if (t > EPSILON && t < min_t) {
                        min_t = t;
                        hit_plane_idx = k;
                    }
                }
            }

            if (hit_plane_idx == -1) {
                // Ray escaped the scene, no more reflections
                break; 
            }

            // --- Update Ray State ---
            Plane hit_plane = g_data.room_planes[hit_plane_idx];

            // 1. Move ray to intersection point
            current_ray.origin = vec_add(current_ray.origin, vec_scale(current_ray.direction, min_t));

            // 2. Absorb energy
            double absorbed = current_ray.energy * hit_plane.absorption_coeff;
            total_absorbed += absorbed;
            current_ray.energy -= absorbed;
            
            if (current_ray.energy < 1e-4) { // Energy threshold to terminate ray
                break;
            }

            // 3. Reflect direction
            Vector3D reflection = vec_sub(current_ray.direction, vec_scale(hit_plane.normal, 2.0 * vec_dot(current_ray.direction, hit_plane.normal)));
            current_ray.direction = vec_normalize(reflection); 

            // 4. Nudge ray off surface to avoid self-intersection
            current_ray.origin = vec_add(current_ray.origin, vec_scale(current_ray.direction, EPSILON * 10));
        }
    }

    g_data.total_absorbed_energy = total_absorbed;
}

void cleanup() {
    free(g_data.rays);
    free(g_data.room_planes);
    g_data.rays = NULL;
    g_data.room_planes = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%.4f\n", g_data.total_absorbed_energy);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
