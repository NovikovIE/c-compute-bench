#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h> // For memset

// --- MERSENNE TWISTER (Verbatim as required) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND PARAMETERS ---

typedef struct {
    int width;
    int height;
    int blur_passes;
    float brightness_threshold; // Converted from an integer parameter

    float *original_image;    // The initial random image
    float *bright_pass_image; // Image containing only pixels above the brightness threshold
    float *temp_blur_buffer;  // Intermediate buffer for the separable blur
    float *final_image;       // The resulting image after adding the bloom effect

    double result_accumulator; // Accumulator for the final result to prevent dead code elimination
} BenchmarkData;

static BenchmarkData g_data;
const int BLUR_RADIUS = 5; // Fixed blur radius for simplicity

// --- HELPER FUNCTIONS ---

// Helper to clamp a value within a range. Used for handling image boundaries during blur.
static inline int clamp(int val, int min, int max) {
    if (val < min) return min;
    if (val > max) return max;
    return val;
}

// --- BENCHMARK CORE FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <image_width> <image_height> <blur_passes> <brightness_threshold> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.width = atoi(argv[1]);
    g_data.height = atoi(argv[2]);
    g_data.blur_passes = atoi(argv[3]);
    // The threshold is passed as an integer (0-255) and converted to a float (0.0-1.0)
    g_data.brightness_threshold = (float)atoi(argv[4]) / 255.0f;
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    if (g_data.width <= 0 || g_data.height <= 0 || g_data.blur_passes < 0) {
        fprintf(stderr, "FATAL: Invalid parameters.\n");
        exit(1);
    }

    size_t num_pixels = (size_t)g_data.width * g_data.height;
    size_t image_bytes = num_pixels * 3 * sizeof(float); // 3 channels (R, G, B)

    g_data.original_image = (float*)malloc(image_bytes);
    g_data.bright_pass_image = (float*)malloc(image_bytes);
    g_data.temp_blur_buffer = (float*)malloc(image_bytes);
    g_data.final_image = (float*)malloc(image_bytes);

    if (!g_data.original_image || !g_data.bright_pass_image || !g_data.temp_blur_buffer || !g_data.final_image) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a random initial image with float values between 0.0 and 1.0
    for (size_t i = 0; i < num_pixels * 3; ++i) {
        g_data.original_image[i] = (float)mt_rand() / (float)UINT32_MAX;
    }
}

void run_computation() {
    size_t width = g_data.width;
    size_t height = g_data.height;
    size_t num_pixels = width * height;

    // --- 1. Bright Pass Extraction ---
    // Creates an image containing only the pixels from the original that are brighter than the threshold.
    for (size_t i = 0; i < num_pixels; ++i) {
        size_t base_idx = i * 3;
        float r = g_data.original_image[base_idx + 0];
        float g = g_data.original_image[base_idx + 1];
        float b = g_data.original_image[base_idx + 2];

        // Calculate luminance (a measure of brightness)
        float luminance = 0.299f * r + 0.587f * g + 0.114f * b;

        if (luminance > g_data.brightness_threshold) {
            g_data.bright_pass_image[base_idx + 0] = r;
            g_data.bright_pass_image[base_idx + 1] = g;
            g_data.bright_pass_image[base_idx + 2] = b;
        } else {
            g_data.bright_pass_image[base_idx + 0] = 0.0f;
            g_data.bright_pass_image[base_idx + 1] = 0.0f;
            g_data.bright_pass_image[base_idx + 2] = 0.0f;
        }
    }

    // --- 2. Separable Box Blur ---
    // For each pass, we perform a horizontal blur into a temporary buffer,
    // and then a vertical blur from the temp buffer back into the source buffer.
    float *blur_buffer = g_data.bright_pass_image;
    float *temp_buffer = g_data.temp_blur_buffer;
    float blur_kernel_divisor = (float)(2 * BLUR_RADIUS + 1);

    for (int p = 0; p < g_data.blur_passes; ++p) {
        // Horizontal Pass (from blur_buffer to temp_buffer)
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                float sum_r = 0.0f, sum_g = 0.0f, sum_b = 0.0f;
                for (int kx = -BLUR_RADIUS; kx <= BLUR_RADIUS; ++kx) {
                    int sample_x = clamp(x + kx, 0, width - 1);
                    size_t sample_idx = ((size_t)y * width + sample_x) * 3;
                    sum_r += blur_buffer[sample_idx + 0];
                    sum_g += blur_buffer[sample_idx + 1];
                    sum_b += blur_buffer[sample_idx + 2];
                }
                size_t out_idx = ((size_t)y * width + x) * 3;
                temp_buffer[out_idx + 0] = sum_r / blur_kernel_divisor;
                temp_buffer[out_idx + 1] = sum_g / blur_kernel_divisor;
                temp_buffer[out_idx + 2] = sum_b / blur_kernel_divisor;
            }
        }
        
        // Vertical Pass (from temp_buffer to blur_buffer)
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                float sum_r = 0.0f, sum_g = 0.0f, sum_b = 0.0f;
                 for (int ky = -BLUR_RADIUS; ky <= BLUR_RADIUS; ++ky) {
                    int sample_y = clamp(y + ky, 0, height - 1);
                    size_t sample_idx = ((size_t)sample_y * width + x) * 3;
                    sum_r += temp_buffer[sample_idx + 0];
                    sum_g += temp_buffer[sample_idx + 1];
                    sum_b += temp_buffer[sample_idx + 2];
                }
                size_t out_idx = ((size_t)y * width + x) * 3;
                blur_buffer[out_idx + 0] = sum_r / blur_kernel_divisor;
                blur_buffer[out_idx + 1] = sum_g / blur_kernel_divisor;
                blur_buffer[out_idx + 2] = sum_b / blur_kernel_divisor;
            }
        }
    }
    // The final blurred result is now in g_data.bright_pass_image.

    // --- 3. Additive Blending ---
    // Add the blurred bright pass (bloom) back to the original image.
    for (size_t i = 0; i < num_pixels * 3; ++i) {
        float val = g_data.original_image[i] + g_data.bright_pass_image[i];
        // Clamp the final pixel values to the [0.0, 1.0] range
        g_data.final_image[i] = val > 1.0f ? 1.0f : val;
    }

    // --- 4. Final Accumulation ---
    // Sum up all pixel components of the final image to produce a single result.
    // This ensures the compiler doesn't optimize away the computation.
    g_data.result_accumulator = 0.0;
    for (size_t i = 0; i < num_pixels * 3; ++i) {
        g_data.result_accumulator += g_data.final_image[i];
    }
}

void cleanup() {
    free(g_data.original_image);
    free(g_data.bright_pass_image);
    free(g_data.temp_blur_buffer);
    free(g_data.final_image);

    g_data.original_image = NULL;
    g_data.bright_pass_image = NULL;
    g_data.temp_blur_buffer = NULL;
    g_data.final_image = NULL;
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", g_data.result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
