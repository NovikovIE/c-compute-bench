#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) --- Do Not Modify ---
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
// --- End of MT19937 ---

// --- Benchmark Data Structures ---
typedef struct {
    float x1, y1, x2, y2;
} Rectangle;

// --- Global Variables ---
static int NUM_RECTANGLES;
static Rectangle* rectangles;
static long long total_intersections; // Use long long to avoid overflow

// --- Utility for Data Generation ---
static inline float random_float(float max_val) {
    // Generates a float between 0.0 and max_val
    return ((float)mt_rand() / (float)UINT32_MAX) * max_val;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_rectangles> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_RECTANGLES = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (NUM_RECTANGLES <= 0) {
        fprintf(stderr, "Error: num_rectangles must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    rectangles = (Rectangle*)malloc(NUM_RECTANGLES * sizeof(Rectangle));
    if (!rectangles) {
        fprintf(stderr, "Error: Memory allocation failed for rectangles.\n");
        exit(1);
    }

    const float COORD_SPACE_SIZE = 10000.0f;
    for (int i = 0; i < NUM_RECTANGLES; ++i) {
        float p1x = random_float(COORD_SPACE_SIZE);
        float p1y = random_float(COORD_SPACE_SIZE);
        float p2x = random_float(COORD_SPACE_SIZE);
        float p2y = random_float(COORD_SPACE_SIZE);

        // Ensure x1 < x2 and y1 < y2
        rectangles[i].x1 = fminf(p1x, p2x);
        rectangles[i].y1 = fminf(p1y, p2y);
        rectangles[i].x2 = fmaxf(p1x, p2x);
        rectangles[i].y2 = fmaxf(p1y, p2y);
    }
}

void run_computation() {
    total_intersections = 0;
    for (int i = 0; i < NUM_RECTANGLES; ++i) {
        for (int j = i + 1; j < NUM_RECTANGLES; ++j) {
            Rectangle r1 = rectangles[i];
            Rectangle r2 = rectangles[j];

            // Check for intersection. Two rectangles intersect if their projections
            // on both the x and y axes overlap.
            if (r1.x1 < r2.x2 && r1.x2 > r2.x1 && r1.y1 < r2.y2 && r1.y2 > r2.y1) {
                total_intersections++;
            }
        }
    }
}

void cleanup() {
    free(rectangles);
    rectangles = NULL;
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

    // Print the final result to stdout
    printf("%lld\n", total_intersections);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
