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
// --- End Mersenne Twister ---

// --- Benchmark Globals ---
typedef struct {
    // Parameters
    int volume_width;
    int volume_height;
    int volume_depth;
    int samples_per_ray;

    // Data
    float* volume_data; // 3D volume grid
    float* image_data;  // 2D output image

    // Result
    double accumulated_result;
} BenchmarkData;

static BenchmarkData g_data;

// --- Helper Functions for Computation ---

// Clamp a value to a range
static inline float clamp(float val, float low, float high) {
    return fmaxf(low, fminf(val, high));
}

// Linear interpolation
static inline float lerp(float a, float b, float t) {
    return a + t * (b - a);
}

// Sample the volume using trilinear interpolation
static float sample_volume(float x, float y, float z) {
    // Clamp coordinates to be within the volume bounds
    x = clamp(x, 0.0f, g_data.volume_width - 1.001f);
    y = clamp(y, 0.0f, g_data.volume_height - 1.001f);
    z = clamp(z, 0.0f, g_data.volume_depth - 1.001f);

    int x0 = (int)x;
    int y0 = (int)y;
    int z0 = (int)z;

    int x1 = x0 + 1;
    int y1 = y0 + 1;
    int z1 = z0 + 1;

    float tx = x - x0;
    float ty = y - y0;
    float tz = z - z0;
    
    // Indices for the 8 corner voxels
    const int w = g_data.volume_width;
    const int h = g_data.volume_height;
    size_t idx000 = (size_t)z0 * w * h + (size_t)y0 * w + x0;
    size_t idx100 = (size_t)z0 * w * h + (size_t)y0 * w + x1;
    size_t idx010 = (size_t)z0 * w * h + (size_t)y1 * w + x0;
    size_t idx110 = (size_t)z0 * w * h + (size_t)y1 * w + x1;
    size_t idx001 = (size_t)z1 * w * h + (size_t)y0 * w + x0;
    size_t idx101 = (size_t)z1 * w * h + (size_t)y0 * w + x1;
    size_t idx011 = (size_t)z1 * w * h + (size_t)y1 * w + x0;
    size_t idx111 = (size_t)z1 * w * h + (size_t)y1 * w + x1;

    // Fetch the 8 voxel values
    float v000 = g_data.volume_data[idx000];
    float v100 = g_data.volume_data[idx100];
    float v010 = g_data.volume_data[idx010];
    float v110 = g_data.volume_data[idx110];
    float v001 = g_data.volume_data[idx001];
    float v101 = g_data.volume_data[idx101];
    float v011 = g_data.volume_data[idx011];
    float v111 = g_data.volume_data[idx111];

    // Interpolate along x
    float c00 = lerp(v000, v100, tx);
    float c10 = lerp(v010, v110, tx);
    float c01 = lerp(v001, v101, tx);
    float c11 = lerp(v011, v111, tx);
    
    // Interpolate along y
    float c0 = lerp(c00, c10, ty);
    float c1 = lerp(c01, c11, ty);

    // Interpolate along z
    return lerp(c0, c1, tz);
}


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <volume_width> <volume_height> <volume_depth> <samples_per_ray> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.volume_width = atoi(argv[1]);
    g_data.volume_height = atoi(argv[2]);
    g_data.volume_depth = atoi(argv[3]);
    g_data.samples_per_ray = atoi(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);
    
    mt_seed(seed);

    size_t volume_size = (size_t)g_data.volume_width * g_data.volume_height * g_data.volume_depth;
    size_t image_size = (size_t)g_data.volume_width * g_data.volume_height;

    g_data.volume_data = (float*)calloc(volume_size, sizeof(float));
    g_data.image_data = (float*)malloc(image_size * sizeof(float));
    
    if (!g_data.volume_data || !g_data.image_data) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    
    g_data.accumulated_result = 0.0;
    
    // Procedurally generate a density volume made of several spheres
    const int num_spheres = 25;
    for (int i = 0; i < num_spheres; ++i) {
        float r = ((mt_rand() % 1000) / 1000.0f) * (fminf(g_data.volume_width, fminf(g_data.volume_height, g_data.volume_depth)) * 0.15f) + 5.0f;
        float r2 = r * r;
        float cx = ((mt_rand() % 1000) / 1000.0f) * g_data.volume_width;
        float cy = ((mt_rand() % 1000) / 1000.0f) * g_data.volume_height;
        float cz = ((mt_rand() % 1000) / 1000.0f) * g_data.volume_depth;
        float density = ((mt_rand() % 1000) / 1000.0f) * 0.1f;
        
        int x_min = (int)fmaxf(0, cx - r);
        int x_max = (int)fminf(g_data.volume_width - 1, cx + r);
        int y_min = (int)fmaxf(0, cy - r);
        int y_max = (int)fminf(g_data.volume_height - 1, cy + r);
        int z_min = (int)fmaxf(0, cz - r);
        int z_max = (int)fminf(g_data.volume_depth - 1, cz + r);

        for (int z = z_min; z <= z_max; ++z) {
            for (int y = y_min; y <= y_max; ++y) {
                for (int x = x_min; x <= x_max; ++x) {
                    float dx = x - cx;
                    float dy = y - cy;
                    float dz = z - cz;
                    if (dx*dx + dy*dy + dz*dz < r2) {
                        size_t index = (size_t)z * g_data.volume_width * g_data.volume_height + (size_t)y * g_data.volume_width + x;
                        g_data.volume_data[index] += density;
                    }
                }
            }
        }
    }
}

void run_computation() {
    const int width = g_data.volume_width;
    const int height = g_data.volume_height;
    const int depth = g_data.volume_depth;
    const int samples = g_data.samples_per_ray;
    
    const float step = (float)depth / (float)samples;

    // Orthographic view along the Z-axis
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            float total_density = 0.0f;
            float ray_x = (float)x + 0.5f;
            float ray_y = (float)y + 0.5f;

            // March ray through the volume
            for (int s = 0; s < samples; ++s) {
                float ray_z = s * step;
                total_density += sample_volume(ray_x, ray_y, ray_z);
            }
            g_data.image_data[y * width + x] = total_density;
        }
    }

    // Accumulate final result to prevent dead code elimination
    double sum = 0.0;
    size_t image_size = (size_t)width * height;
    for (size_t i = 0; i < image_size; ++i) {
        sum += g_data.image_data[i];
    }
    g_data.accumulated_result = sum;
}

void cleanup() {
    free(g_data.volume_data);
    free(g_data.image_data);
    g_data.volume_data = NULL;
    g_data.image_data = NULL;
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print result to stdout
    printf("%f\n", g_data.accumulated_result);
    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
