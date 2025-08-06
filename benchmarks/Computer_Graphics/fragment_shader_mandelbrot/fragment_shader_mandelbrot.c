#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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
// --- End of Mersenne Twister ---

// Global structure to hold benchmark data
typedef struct {
    int image_width;
    int image_height;
    int max_iterations;
    int* image_buffer; // Stores iteration count for each pixel
    long long total_iterations_sum; // Accumulated result
} BenchmarkData;

BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <image_width> <image_height> <max_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.image_width = atoi(argv[1]);
    g_data.image_height = atoi(argv[2]);
    g_data.max_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    g_data.image_buffer = (int*)malloc(g_data.image_width * g_data.image_height * sizeof(int));
    if (g_data.image_buffer == NULL) {
        fprintf(stderr, "Failed to allocate memory for image buffer\n");
        exit(1);
    }

    g_data.total_iterations_sum = 0;
}

void run_computation() {
    long long sum = 0;
    
    for (int py = 0; py < g_data.image_height; py++) {
        for (int px = 0; px < g_data.image_width; px++) {
            // Map pixel coordinate to complex plane
            // x0 in [-2.0, 1.0], y0 in [-1.5, 1.5]
            double x0 = (px / (double)(g_data.image_width)) * 3.0 - 2.0;
            double y0 = (py / (double)(g_data.image_height)) * 3.0 - 1.5;

            double x = 0.0;
            double y = 0.0;
            int iteration = 0;

            while (x*x + y*y <= 4.0 && iteration < g_data.max_iterations) {
                double xtemp = x*x - y*y + x0;
                y = 2*x*y + y0;
                x = xtemp;
                iteration++;
            }

            g_data.image_buffer[py * g_data.image_width + px] = iteration;
            sum += iteration;
        }
    }
    g_data.total_iterations_sum = sum;
}

void cleanup() {
    if (g_data.image_buffer) {
        free(g_data.image_buffer);
        g_data.image_buffer = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data.total_iterations_sum);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
