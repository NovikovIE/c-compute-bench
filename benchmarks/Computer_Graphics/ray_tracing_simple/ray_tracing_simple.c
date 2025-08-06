#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// ---- Mersenne Twister (DO NOT MODIFY) ----
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
// ---- END Mersenne Twister ----

// ---- Benchmark Data Structures ----
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
    double x, y, z;
} Vec3;

typedef struct {
    Vec3 origin;
    Vec3 direction;
} Ray;

typedef struct {
    Vec3 center;
    double radius;
    Vec3 color;
    double specular; // Shininess factor
    double reflective; // 0 to 1
} Sphere;

typedef struct {
    Vec3 position;
    Vec3 color;
    double intensity;
} Light;

// Global parameters
static int g_width;
static int g_height;
static int g_num_objects;
static int g_max_depth;
static const int g_num_lights = 2;

// Global scene data
static Sphere* g_spheres = NULL;
static Light* g_lights = NULL;
static unsigned char* g_image_buffer = NULL;
static unsigned long long g_checksum = 0;

// ---- Vector Math ----
Vec3 vec3_add(Vec3 a, Vec3 b) { return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z}; }
Vec3 vec3_sub(Vec3 a, Vec3 b) { return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
Vec3 vec3_scale(Vec3 v, double s) { return (Vec3){v.x * s, v.y * s, v.z * s}; }
double vec3_dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double vec3_length_sq(Vec3 v) { return vec3_dot(v, v); }
double vec3_length(Vec3 v) { return sqrt(vec3_length_sq(v)); }
Vec3 vec3_normalize(Vec3 v) { double len = vec3_length(v); return len > 1e-9 ? vec3_scale(v, 1.0 / len) : v; }

// ---- Random Number Generation ----
double rand_double(double min, double max) {
    return min + (max - min) * (double)mt_rand() / (double)UINT32_MAX;
}

// ---- Ray Tracing Core ----
int intersect_sphere(const Ray* ray, const Sphere* sphere, double* t) {
    Vec3 oc = vec3_sub(ray->origin, sphere->center);
    double a = vec3_dot(ray->direction, ray->direction);
    double b = 2.0 * vec3_dot(oc, ray->direction);
    double c = vec3_dot(oc, oc) - sphere->radius * sphere->radius;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return 0; // No intersection
    }

    double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
    if (t1 > 1e-4) { *t = t1; return 1; }
    
    double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
    if (t2 > 1e-4) { *t = t2; return 1; }
    
    return 0;
}

int find_closest_object(const Ray* ray, double* closest_t, int* hit_object_idx) {
    *closest_t = 1e20;
    *hit_object_idx = -1;
    for (int i = 0; i < g_num_objects; ++i) {
        double t;
        if (intersect_sphere(ray, &g_spheres[i], &t)) {
            if (t < *closest_t) {
                *closest_t = t;
                *hit_object_idx = i;
            }
        }
    }
    return (*hit_object_idx != -1);
}

Vec3 trace(const Ray* ray, int depth) {
    double closest_t;
    int hit_object_idx;

    if (depth >= g_max_depth || !find_closest_object(ray, &closest_t, &hit_object_idx)) {
        return (Vec3){0.0, 0.0, 0.0};
    }
    
    const Sphere* hit_sphere = &g_spheres[hit_object_idx];
    Vec3 hit_point = vec3_add(ray->origin, vec3_scale(ray->direction, closest_t));
    Vec3 normal = vec3_normalize(vec3_sub(hit_point, hit_sphere->center));

    Vec3 total_color = {0.0, 0.0, 0.0};
    for (int i = 0; i < g_num_lights; ++i) {
        Vec3 light_dir = vec3_normalize(vec3_sub(g_lights[i].position, hit_point));
        Ray shadow_ray = {vec3_add(hit_point, vec3_scale(normal, 1e-4)), light_dir};
        double shadow_t;
        int shadow_hit_idx;
        if (find_closest_object(&shadow_ray, &shadow_t, &shadow_hit_idx)) {
            double light_dist_sq = vec3_length_sq(vec3_sub(g_lights[i].position, hit_point));
            if (shadow_t * shadow_t < light_dist_sq) continue;
        }

        double diff = fmax(0.0, vec3_dot(normal, light_dir));
        Vec3 diffuse_color = vec3_scale(hit_sphere->color, diff * g_lights[i].intensity);
        
        Vec3 view_dir = vec3_normalize(vec3_scale(ray->direction, -1.0));
        Vec3 half_vector = vec3_normalize(vec3_add(light_dir, view_dir));
        double spec_angle = fmax(0.0, vec3_dot(normal, half_vector));
        double specular = pow(spec_angle, hit_sphere->specular);
        Vec3 specular_color = vec3_scale(g_lights[i].color, specular * g_lights[i].intensity);

        total_color = vec3_add(total_color, diffuse_color);
        total_color = vec3_add(total_color, specular_color);
    }
    
    Vec3 reflected_color = {0.0, 0.0, 0.0};
    if (hit_sphere->reflective > 0) {
        Vec3 reflection_dir = vec3_sub(ray->direction, vec3_scale(normal, 2.0 * vec3_dot(ray->direction, normal)));
        Ray reflected_ray = {vec3_add(hit_point, vec3_scale(normal, 1e-4)), reflection_dir};
        reflected_color = trace(&reflected_ray, depth + 1);
    }

    Vec3 local_color = vec3_add((Vec3){0.1, 0.1, 0.1}, total_color);
    Vec3 final_color = vec3_add(vec3_scale(local_color, 1.0 - hit_sphere->reflective), 
                               vec3_scale(reflected_color, hit_sphere->reflective));

    return final_color;
}

// ---- Benchmark Setup, Run, Cleanup ----
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s image_width image_height num_objects max_recursion_depth seed\n", argv[0]);
        exit(1);
    }
    g_width = atoi(argv[1]);
    g_height = atoi(argv[2]);
    g_num_objects = atoi(argv[3]);
    g_max_depth = atoi(argv[4]);
    uint32_t seed = (uint32_t)atol(argv[5]);

    mt_seed(seed);

    g_spheres = (Sphere*)malloc(g_num_objects * sizeof(Sphere));
    g_lights = (Light*)malloc(g_num_lights * sizeof(Light));
    g_image_buffer = (unsigned char*)malloc((size_t)g_width * g_height * 3 * sizeof(unsigned char));

    if (!g_spheres || !g_lights || !g_image_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < g_num_objects; ++i) {
        g_spheres[i].center = (Vec3){
            rand_double(-10.0, 10.0),
            rand_double(-10.0, 10.0),
            rand_double(5.0, 25.0)
        };
        g_spheres[i].radius = rand_double(0.5, 2.0);
        g_spheres[i].color = (Vec3){rand_double(0.1, 1.0), rand_double(0.1, 1.0), rand_double(0.1, 1.0)};
        g_spheres[i].specular = rand_double(10.0, 250.0);
        g_spheres[i].reflective = rand_double(0.0, 0.75);
    }
    
    g_lights[0] = (Light){(Vec3){-20, 20, -10}, (Vec3){1,1,1}, 1.0};
    g_lights[1] = (Light){(Vec3){20, 20, -10}, (Vec3){1,1,1}, 0.8};
}

void run_computation() {
    g_checksum = 0;
    Vec3 eye = {0, 0, -5};
    double fov = M_PI / 3.0;
    double aspect_ratio = (double)g_width / (double)g_height;
    double tan_half_fov = tan(fov / 2.0);

    for (int j = 0; j < g_height; ++j) {
        for (int i = 0; i < g_width; ++i) {
            double px = (2 * ((i + 0.5) / g_width) - 1) * tan_half_fov * aspect_ratio;
            double py = (1 - 2 * ((j + 0.5) / g_height)) * tan_half_fov;
            
            Vec3 dir = vec3_normalize((Vec3){px, py, 1.0});
            Ray ray = {eye, dir};

            Vec3 color = trace(&ray, 0);

            unsigned char r = (unsigned char)(fmin(1.0, fmax(0.0, color.x)) * 255.0);
            unsigned char g = (unsigned char)(fmin(1.0, fmax(0.0, color.y)) * 255.0);
            unsigned char b = (unsigned char)(fmin(1.0, fmax(0.0, color.z)) * 255.0);
            
            size_t idx = ((size_t)j * g_width + i) * 3;
            g_image_buffer[idx] = r;
            g_image_buffer[idx+1] = g;
            g_image_buffer[idx+2] = b;

            g_checksum += r + g + b;
        }
    }
}

void cleanup() {
    free(g_spheres);
    g_spheres = NULL;
    free(g_lights);
    g_lights = NULL;
    free(g_image_buffer);
    g_image_buffer = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%llu\n", g_checksum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
