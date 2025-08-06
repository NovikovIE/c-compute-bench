#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

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

// Benchmark-specific code

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Data Structures
typedef struct {
    double x, y, z;
} Vec3;

typedef struct {
    Vec3 origin, dir;
} Ray;

typedef struct {
    Vec3 center;
    double radius;
    Vec3 color;
    Vec3 emission;
} Sphere;

typedef struct {
    Vec3 pos;
    Vec3 power;
} Photon;

typedef struct {
    double dist_sq;
    int photon_idx;
} Neighbor;

// Global variables for benchmark data
static int image_width;
static int image_height;
static int num_photons_to_emit;
static int k_nearest_neighbors;

static Sphere *scene = NULL;
static const int num_spheres = 4;
static Photon *photon_map = NULL;
static int photons_stored = 0;
static Vec3 *image_buffer = NULL;
static unsigned long long final_result = 0;

// Helper functions
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

Vec3 vec_add(Vec3 a, Vec3 b) { return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 vec_sub(Vec3 a, Vec3 b) { return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 vec_scale(Vec3 v, double s) { return (Vec3){v.x * s, v.y * s, v.z * s}; }
double vec_dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double vec_length_sq(Vec3 v) { return vec_dot(v, v); }
Vec3 vec_normalize(Vec3 v) {
    double len = sqrt(vec_length_sq(v));
    return (len > 1e-9) ? vec_scale(v, 1.0 / len) : v;
}

int intersect_sphere(const Ray *r, const Sphere *s, double *t) {
    Vec3 oc = vec_sub(r->origin, s->center);
    double b = vec_dot(oc, r->dir);
    double c = vec_dot(oc, oc) - s->radius * s->radius;
    double disc = b * b - c;
    if (disc < 0) return 0;
    disc = sqrt(disc);
    double t0 = -b - disc;
    double t1 = -b + disc;
    if (t0 > 1e-4) {
        *t = t0;
        return 1;
    } else if (t1 > 1e-4) {
        *t = t1;
        return 1;
    }
    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s image_width image_height num_photons k_neighbors seed\n", argv[0]);
        exit(1);
    }

    image_width = atoi(argv[1]);
    image_height = atoi(argv[2]);
    num_photons_to_emit = atoi(argv[3]);
    k_nearest_neighbors = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    scene = (Sphere*)malloc(sizeof(Sphere) * num_spheres);
    photon_map = (Photon*)malloc(sizeof(Photon) * num_photons_to_emit);
    image_buffer = (Vec3*)malloc(sizeof(Vec3) * image_width * image_height);

    if (!scene || !photon_map || !image_buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Scene: Left wall, Right wall, Back wall, Light source sphere
    scene[0] = (Sphere){{1e5 + 1, 40.8, 81.6}, 1e5, {0.75, 0.25, 0.25}, {0,0,0}}; // Left
    scene[1] = (Sphere){{-1e5 + 99, 40.8, 81.6}, 1e5, {0.25, 0.25, 0.75}, {0,0,0}}; // Right
    scene[2] = (Sphere){{50, 40.8, 1e5}, 1e5, {0.75, 0.75, 0.75}, {0,0,0}}; // Back
    scene[3] = (Sphere){{50, 681.6 - .27, 81.6}, 600, {0,0,0}, {12,12,12}}; // Ground/Floor
}

void run_computation() {
    // Phase 1: Photon Emission
    photons_stored = 0;
    Vec3 light_pos = {50, 60.0, 81.6};
    Vec3 light_power = {2500, 2500, 2500};

    for (int i = 0; i < num_photons_to_emit; ++i) {
        double phi = 2.0 * M_PI * rand_double();
        double cos_theta = 1.0 - rand_double();
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        Vec3 dir = {cos(phi) * sin_theta, cos_theta, sin(phi) * sin_theta};
        Ray photon_ray = {light_pos, dir};
        
        double t = DBL_MAX;
        int hit_idx = -1;
        for (int s = 0; s < num_spheres; ++s) {
            double d;
            if (intersect_sphere(&photon_ray, &scene[s], &d) && d < t) {
                t = d;
                hit_idx = s;
            }
        }

        if (hit_idx != -1) {
            Vec3 hit_pos = vec_add(photon_ray.origin, vec_scale(photon_ray.dir, t));
             photon_map[photons_stored].pos = hit_pos;
             photon_map[photons_stored].power = vec_scale(light_power, 1.0 / num_photons_to_emit);
             photons_stored++;
        }
    }

    // Phase 2: Rendering
    Ray camera = {{50, 52, 295.6}, vec_normalize((Vec3){0, -0.042612, -1})};
    double fov = 0.5;
    Vec3 cx = {image_width * fov / image_height, 0, 0};
    Vec3 cy = vec_scale(vec_normalize((Vec3){cx.y*camera.dir.z - cx.z*camera.dir.y, cx.z*camera.dir.x - cx.x*camera.dir.z, cx.x*camera.dir.y - cx.y*camera.dir.x}), fov);

    for (int y = 0; y < image_height; ++y) {
        for (int x = 0; x < image_width; ++x) {
            int idx = y * image_width + x;
            double dx = (double)x / image_width - 0.5;
            double dy = (double)y / image_height - 0.5;
            Vec3 dir = vec_add(camera.dir, vec_add(vec_scale(cx, dx), vec_scale(cy, dy)));
            Ray r = {camera.origin, vec_normalize(dir)};
            
            double t = DBL_MAX;
            int hit_idx = -1;
            for (int s = 0; s < num_spheres; ++s) {
                double d;
                if (intersect_sphere(&r, &scene[s], &d) && d < t) {
                    t = d;
                    hit_idx = s;
                }
            }

            if (hit_idx == -1) {
                image_buffer[idx] = (Vec3){0, 0, 0};
                continue;
            }
            
            const Sphere *obj = &scene[hit_idx];
            Vec3 hit_pos = vec_add(r.origin, vec_scale(r.dir, t));

            // KNN Photon Search
            Neighbor* neighbors = (Neighbor*)malloc(k_nearest_neighbors * sizeof(Neighbor));
            for(int k=0; k<k_nearest_neighbors; ++k) neighbors[k].dist_sq = DBL_MAX;

            for (int p_idx = 0; p_idx < photons_stored; ++p_idx) {
                double dist2 = vec_length_sq(vec_sub(hit_pos, photon_map[p_idx].pos));
                if (dist2 < neighbors[k_nearest_neighbors - 1].dist_sq) {
                    neighbors[k_nearest_neighbors-1].dist_sq = dist2;
                    neighbors[k_nearest_neighbors-1].photon_idx = p_idx;

                    // Insertion sort to keep neighbors sorted by distance
                    for (int k = k_nearest_neighbors-2; k >= 0; --k) {
                        if (neighbors[k+1].dist_sq < neighbors[k].dist_sq) {
                            Neighbor temp = neighbors[k];
                            neighbors[k] = neighbors[k+1];
                            neighbors[k+1] = temp;
                        } else {
                            break;
                        }
                    }
                }
            }

            Vec3 accumulated_power = {0,0,0};
            double max_dist_sq = neighbors[k_nearest_neighbors-1].dist_sq;
            
            for(int k=0; k<k_nearest_neighbors; ++k) {
                if(neighbors[k].dist_sq < DBL_MAX) {
                    accumulated_power = vec_add(accumulated_power, photon_map[neighbors[k].photon_idx].power);
                }
            }
            free(neighbors);

            double area = M_PI * max_dist_sq;
            if (area > 1e-6) {
                 Vec3 radiance = vec_scale(accumulated_power, 1.0 / area);
                 image_buffer[idx] = (Vec3){radiance.x * obj->color.x, radiance.y * obj->color.y, radiance.z * obj->color.z};
            } else {
                image_buffer[idx] = (Vec3){0,0,0};
            }
        }
    }

    // Accumulate final result to prevent dead code elimination
    for (int i = 0; i < image_width * image_height; ++i) {
        double r = fmin(1.0, fmax(0.0, image_buffer[i].x));
        double g = fmin(1.0, fmax(0.0, image_buffer[i].y));
        double b = fmin(1.0, fmax(0.0, image_buffer[i].z));
        final_result += (unsigned long long)(r * 255.0);
        final_result += (unsigned long long)(g * 255.0);
        final_result += (unsigned long long)(b * 255.0);
    }
}

void cleanup() {
    if (scene) free(scene);
    if (photon_map) free(photon_map);
    if (image_buffer) free(image_buffer);
    scene = NULL;
    photon_map = NULL;
    image_buffer = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%llu\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
