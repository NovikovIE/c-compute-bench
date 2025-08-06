#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// Benchmark parameters and data structures
long grid_width, grid_height, num_queries;
float *grid_values;
float *query_x, *query_y;
double accumulator = 0.0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_width grid_height num_queries seed\n", argv[0]);
        exit(1);
    }

    grid_width = atol(argv[1]);
    grid_height = atol(argv[2]);
    num_queries = atol(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (grid_width <= 1 || grid_height <= 1 || num_queries <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters. grid_width/height must be > 1 and num_queries > 0.\n");
        exit(1);
    }

    grid_values = (float *)malloc(grid_width * grid_height * sizeof(float));
    query_x = (float *)malloc(num_queries * sizeof(float));
    query_y = (float *)malloc(num_queries * sizeof(float));

    if (!grid_values || !query_x || !query_y) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Populate the grid with random values between 0.0 and 1.0
    for (long i = 0; i < grid_width * grid_height; ++i) {
        grid_values[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    // Generate random query points within the grid boundaries
    // The a small epsilon (1.0001f) ensures the point is never on the far boundary,
    // simplifying the interpolation logic (x1+1 is always a valid index).
    float max_x = (float)grid_width - 1.0001f;
    float max_y = (float)grid_height - 1.0001f;
    
    for (long i = 0; i < num_queries; ++i) {
        query_x[i] = ((float)mt_rand() / (float)UINT32_MAX) * max_x;
        query_y[i] = ((float)mt_rand() / (float)UINT32_MAX) * max_y;
    }
}

void run_computation() {
    double local_accumulator = 0.0;

    for (long i = 0; i < num_queries; ++i) {
        float qx = query_x[i];
        float qy = query_y[i];

        // Get the integer coordinates of the top-left corner of the cell
        long x1 = (long)qx;
        long y1 = (long)qy;
        long x2 = x1 + 1;
        long y2 = y1 + 1;

        // Get the values of the four corner points
        float Q11 = grid_values[y1 * grid_width + x1]; // Top-left
        float Q21 = grid_values[y1 * grid_width + x2]; // Top-right
        float Q12 = grid_values[y2 * grid_width + x1]; // Bottom-left
        float Q22 = grid_values[y2 * grid_width + x2]; // Bottom-right

        // Calculate fractional distances
        float frac_x = qx - (float)x1;
        float frac_y = qy - (float)y1;

        // Interpolate along the x-axis for the top and bottom edges
        float interp_top = Q11 * (1.0f - frac_x) + Q21 * frac_x;
        float interp_bottom = Q12 * (1.0f - frac_x) + Q22 * frac_x;

        // Interpolate along the y-axis between the two intermediate results
        float result = interp_top * (1.0f - frac_y) + interp_bottom * frac_y;

        local_accumulator += result;
    }
    accumulator = local_accumulator;
}

void cleanup() {
    free(grid_values);
    free(query_x);
    free(query_y);
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
    printf("%f\n", accumulator);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
