#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
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

// Benchmark parameters and data
static int TEXTURE_WIDTH;
static int TEXTURE_HEIGHT;
static int TEXTURE_DEPTH;
static int NUM_OCTAVES;
static float *texture_data;
static int *p; // Permutation table

// Result accumulator
static double final_hash = 0.0;

// --- Perlin Noise Implementation ---

static double fade(double t) {
    return t * t * t * (t * (t * 6 - 15) + 10);
}

static double lerp(double t, double a, double b) {
    return a + t * (b - a);
}

static double grad(int hash, double x, double y, double z) {
    int h = hash & 15;      
    double u = h < 8 ? x : y;
    double v = h < 4 ? y : (h == 12 || h == 14 ? x : z);
    return ((h & 1) == 0 ? u : -u) + ((h & 2) == 0 ? v : -v);
}

static double perlin_noise_3d(double x, double y, double z) {
    int X = (int)floor(x) & 255;
    int Y = (int)floor(y) & 255;
    int Z = (int)floor(z) & 255;

    x -= floor(x);
    y -= floor(y);
    z -= floor(z);

    double u = fade(x);
    double v = fade(y);
    double w = fade(z);

    int A = p[X] + Y;
    int AA = p[A] + Z;
    int AB = p[A + 1] + Z;
    int B = p[X + 1] + Y;
    int BA = p[B] + Z;
    int BB = p[B + 1] + Z;

    double res = lerp(w, lerp(v, lerp(u, grad(p[AA], x, y, z),
                                       grad(p[BA], x - 1, y, z)),
                               lerp(u, grad(p[AB], x, y - 1, z),
                                       grad(p[BB], x - 1, y - 1, z))),
                      lerp(v, lerp(u, grad(p[AA + 1], x, y, z - 1),
                                       grad(p[BA + 1], x - 1, y, z - 1)),
                               lerp(u, grad(p[AB + 1], x, y - 1, z - 1),
                                       grad(p[BB + 1], x - 1, y - 1, z - 1))));
    return (res + 1.0) / 2.0; // Normalize to [0, 1]
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <texture_width> <texture_height> <texture_depth> <num_octaves> <seed>\n", argv[0]);
        exit(1);
    }
    TEXTURE_WIDTH = atoi(argv[1]);
    TEXTURE_HEIGHT = atoi(argv[2]);
    TEXTURE_DEPTH = atoi(argv[3]);
    NUM_OCTAVES = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);
    
    size_t total_voxels = (size_t)TEXTURE_WIDTH * TEXTURE_HEIGHT * TEXTURE_DEPTH;
    texture_data = (float*)malloc(total_voxels * sizeof(float));
    if (!texture_data) {
        fprintf(stderr, "Failed to allocate texture_data\n");
        exit(1);
    }
    
    p = (int*)malloc(512 * sizeof(int));
    if (!p) {
        fprintf(stderr, "Failed to allocate permutation table\n");
        exit(1);
    }
    
    int permutation[256];
    for (int i = 0; i < 256; i++) {
        permutation[i] = i;
    }
    
    for (int i = 255; i > 0; i--) {
        int j = mt_rand() % (i + 1);
        int temp = permutation[i];
        permutation[i] = permutation[j];
        permutation[j] = temp;
    }
    
    for (int i = 0; i < 256; i++) {
        p[i] = p[i + 256] = permutation[i];
    }
}

void run_computation() {
    double total_noise_sum = 0.0;
    double freq_scale_x = 1.0 / TEXTURE_WIDTH;
    double freq_scale_y = 1.0 / TEXTURE_HEIGHT;
    double freq_scale_z = 1.0 / TEXTURE_DEPTH;

    for (int z = 0; z < TEXTURE_DEPTH; ++z) {
        for (int y = 0; y < TEXTURE_HEIGHT; ++y) {
            for (int x = 0; x < TEXTURE_WIDTH; ++x) {
                double noise_val = 0.0;
                double freq = 1.0;
                double amp = 1.0;
                
                for (int i = 0; i < NUM_OCTAVES; ++i) {
                    double sample_x = x * freq_scale_x * freq;
                    double sample_y = y * freq_scale_y * freq;
                    double sample_z = z * freq_scale_z * freq;
                    
                    noise_val += perlin_noise_3d(sample_x, sample_y, sample_z) * amp;
                    
                    freq *= 2.0; // Lacunarity
                    amp *= 0.5;  // Persistence
                }
                
                size_t index = (size_t)z * TEXTURE_WIDTH * TEXTURE_HEIGHT + (size_t)y * TEXTURE_WIDTH + x;
                texture_data[index] = (float)noise_val;
                total_noise_sum += noise_val;
            }
        }
    }
    final_hash = total_noise_sum;
}

void cleanup() {
    free(texture_data);
    free(p);
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
    printf("%f\n", final_hash);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
