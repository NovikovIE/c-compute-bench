#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Data and Globals ---
typedef struct {
    int grid_width;
    int grid_height;
    int num_layers;
    
    // 3D array for raster layers: layers[layer][height][width]
    float*** layers; 
    
    // 1D array for layer weights
    float* weights;
    
    // 2D grid for the final output
    float** output_grid;
    
    // Final result to prevent dead code elimination
    double final_sum; 
} BenchmarkData;

static BenchmarkData g_data;

// --- Utility Functions ---
float rand_float(float min, float max) {
    return min + ((float)mt_rand() / (float)UINT32_MAX) * (max - min);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <grid_width> <grid_height> <num_layers> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.grid_width = atoi(argv[1]);
    g_data.grid_height = atoi(argv[2]);
    g_data.num_layers = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    // Allocate layers (3D array)
    g_data.layers = (float***)malloc(g_data.num_layers * sizeof(float**));
    for (int k = 0; k < g_data.num_layers; k++) {
        g_data.layers[k] = (float**)malloc(g_data.grid_height * sizeof(float*));
        for (int i = 0; i < g_data.grid_height; i++) {
            g_data.layers[k][i] = (float*)malloc(g_data.grid_width * sizeof(float));
        }
    }

    // Allocate weights (1D array)
    g_data.weights = (float*)malloc(g_data.num_layers * sizeof(float));

    // Allocate output grid (2D array)
    g_data.output_grid = (float**)malloc(g_data.grid_height * sizeof(float*));
    for (int i = 0; i < g_data.grid_height; i++) {
        g_data.output_grid[i] = (float*)malloc(g_data.grid_width * sizeof(float));
    }

    // Initialize data with random values
    // Raster values typically represent something like elevation or temperature
    for (int k = 0; k < g_data.num_layers; k++) {
        for (int i = 0; i < g_data.grid_height; i++) {
            for (int j = 0; j < g_data.grid_width; j++) {
                g_data.layers[k][i][j] = rand_float(0.0f, 255.0f); // e.g., 8-bit sensor range
            }
        }
    }
    
    // Weights for combining layers
    for (int k = 0; k < g_data.num_layers; k++) {
        g_data.weights[k] = rand_float(0.1f, 1.0f); // Assign some significance to each layer
    }
}

void run_computation() {
    // Perform the weighted sum overlay
    for (int i = 0; i < g_data.grid_height; i++) {
        for (int j = 0; j < g_data.grid_width; j++) {
            float sum = 0.0f;
            for (int k = 0; k < g_data.num_layers; k++) {
                sum += g_data.layers[k][i][j] * g_data.weights[k];
            }
            g_data.output_grid[i][j] = sum;
        }
    }

    // Calculate a final sum to prevent the compiler from optimizing away the work
    double total_sum = 0.0;
    for (int i = 0; i < g_data.grid_height; i++) {
        for (int j = 0; j < g_data.grid_width; j++) {
            total_sum += g_data.output_grid[i][j];
        }
    }
    g_data.final_sum = total_sum;
}

void cleanup() {
    // Free layers
    for (int k = 0; k < g_data.num_layers; k++) {
        for (int i = 0; i < g_data.grid_height; i++) {
            free(g_data.layers[k][i]);
        }
        free(g_data.layers[k]);
    }
    free(g_data.layers);

    // Free weights
    free(g_data.weights);

    // Free output grid
    for (int i = 0; i < g_data.grid_height; i++) {
        free(g_data.output_grid[i]);
    }
    free(g_data.output_grid);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    // Print result to stdout
    printf("%f\n", g_data.final_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
