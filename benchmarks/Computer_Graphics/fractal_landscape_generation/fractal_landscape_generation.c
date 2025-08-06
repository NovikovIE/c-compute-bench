#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// stdio.h, stdlib.h, time.h should be included at the top of the file

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

/* --- Benchmark Globals --- */
typedef struct {
    int num_recursion_levels;
    float roughness;
    uint32_t seed;
    int grid_size;
    float *heightmap;
    double final_result;
} BenchmarkData;

BenchmarkData g_data;

// Helper to generate a random float between -1.0 and 1.0
float rand_float_neg1_to_1() {
    return ((mt_rand() / (float)UINT32_MAX) * 2.0f) - 1.0f;
}

/* --- Benchmark Functions --- */

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_recursion_levels> <roughness> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_recursion_levels = atoi(argv[1]);
    g_data.roughness = atof(argv[2]);
    g_data.seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_recursion_levels <= 0) {
        fprintf(stderr, "Number of recursion levels must be positive.\n");
        exit(1);
    }

    // Grid size for diamond-square is (2^n) + 1
    g_data.grid_size = (1 << g_data.num_recursion_levels) + 1;
    size_t map_bytes = (size_t)g_data.grid_size * g_data.grid_size * sizeof(float);
    
    g_data.heightmap = (float*)malloc(map_bytes);
    if (!g_data.heightmap) {
        fprintf(stderr, "Failed to allocate memory for heightmap.\n");
        exit(1);
    }

    mt_seed(g_data.seed);

    // Initialize the four corners with random values
    int gs = g_data.grid_size;
    g_data.heightmap[0] = rand_float_neg1_to_1(); // Top-left
    g_data.heightmap[gs - 1] = rand_float_neg1_to_1(); // Top-right
    g_data.heightmap[(gs - 1) * gs] = rand_float_neg1_to_1(); // Bottom-left
    g_data.heightmap[(gs - 1) * gs + (gs - 1)] = rand_float_neg1_to_1(); // Bottom-right
    
    g_data.final_result = 0.0;
}

void run_computation() {
    int gs = g_data.grid_size;
    float *hm = g_data.heightmap;
    float displacement_range = 1.0f; // Initial displacement range

    for (int step = gs - 1; step > 1; step /= 2) {
        int half_step = step / 2;
        
        // Diamond step
        for (int y = 0; y < gs - 1; y += step) {
            for (int x = 0; x < gs - 1; x += step) {
                float c1 = hm[y * gs + x];
                float c2 = hm[y * gs + (x + step)];
                float c3 = hm[(y + step) * gs + x];
                float c4 = hm[(y + step) * gs + (x + step)];
                
                float avg = (c1 + c2 + c3 + c4) / 4.0f;
                float displacement = rand_float_neg1_to_1() * displacement_range;
                
                hm[(y + half_step) * gs + (x + half_step)] = avg + displacement;
            }
        }
        
        // Square step
        for (int y = 0; y < gs; y += half_step) {
            for (int x = (y + half_step) % step; x < gs; x += step) {
                float sum = 0.0f;
                int count = 0;
                
                if (y >= half_step) { sum += hm[(y - half_step) * gs + x]; count++; } // North
                if (y < gs - half_step) { sum += hm[(y + half_step) * gs + x]; count++; } // South
                if (x >= half_step) { sum += hm[y * gs + (x - half_step)]; count++; } // West
                if (x < gs - half_step) { sum += hm[y * gs + (x + half_step)]; count++; } // East
                
                if (count > 0) {
                    float avg = sum / count;
                    float displacement = rand_float_neg1_to_1() * displacement_range;
                    hm[y * gs + x] = avg + displacement;
                }
            }
        }
        
        // Reduce the random displacement range according to roughness
        displacement_range *= powf(2.0f, -g_data.roughness);
    }

    // Calculate a final result to prevent dead code elimination
    double sum = 0.0;
    for (int i = 0; i < gs * gs; ++i) {
        sum += hm[i];
    }
    g_data.final_result = sum;
}

void cleanup() {
    if (g_data.heightmap) {
        free(g_data.heightmap);
        g_data.heightmap = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
