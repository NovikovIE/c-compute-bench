#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// Data structure for an observer point
typedef struct {
    int x;
    int y;
    float height_offset; // Height of observer above the terrain at (x,y)
} Observer;

// --- Global Benchmark Data ---
int dem_width;
int dem_height;
int num_observer_points;

float* dem = NULL; // Digital Elevation Model (raster grid)
Observer* observers = NULL;
long long total_visible_cells = 0; // Final result accumulator

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s dem_width dem_height num_observer_points seed\n", argv[0]);
        exit(1);
    }

    dem_width = atoi(argv[1]);
    dem_height = atoi(argv[2]);
    num_observer_points = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory
    dem = (float*)malloc((long)dem_width * dem_height * sizeof(float));
    if (!dem) {
        fprintf(stderr, "Failed to allocate memory for DEM\n");
        exit(1);
    }

    observers = (Observer*)malloc(num_observer_points * sizeof(Observer));
    if (!observers) {
        fprintf(stderr, "Failed to allocate memory for observers\n");
        free(dem);
        exit(1);
    }

    // Generate random Digital Elevation Model data
    for (int i = 0; i < dem_width * dem_height; ++i) {
        dem[i] = (float)mt_rand() / (float)UINT32_MAX * 1000.0f; // Elevation 0-1000m
    }

    // Generate random observer points
    for (int i = 0; i < num_observer_points; ++i) {
        observers[i].x = mt_rand() % dem_width;
        observers[i].y = mt_rand() % dem_height;
        observers[i].height_offset = 1.5f + (float)mt_rand() / (float)UINT32_MAX * 20.0f; // Height 1.5m to 21.5m
    }
}

void run_computation() {
    total_visible_cells = 0;

    for (int i = 0; i < num_observer_points; ++i) {
        Observer obs = observers[i];
        long obs_idx = (long)obs.y * dem_width + obs.x;
        float obs_height = dem[obs_idx] + obs.height_offset;

        for (int ty = 0; ty < dem_height; ++ty) {
            for (int tx = 0; tx < dem_width; ++tx) {
                if (tx == obs.x && ty == obs.y) {
                    continue; // Skip the observer's own location
                }

                long target_idx = (long)ty * dem_width + tx;
                float target_height = dem[target_idx];

                // Line of Sight (LOS) check using a DDA-like algorithm
                float dx = tx - obs.x;
                float dy = ty - obs.y;

                int steps = (int)fmaxf(fabsf(dx), fabsf(dy));
                if (steps == 0) continue;

                float x_inc = dx / (float)steps;
                float y_inc = dy / (float)steps;

                float curr_x = obs.x;
                float curr_y = obs.y;

                int is_visible = 1;
                for (int k = 1; k < steps; ++k) {
                    curr_x += x_inc;
                    curr_y += y_inc;

                    int ix = (int)roundf(curr_x);
                    int iy = (int)roundf(curr_y);

                    float los_height = obs_height + (target_height - obs_height) * ((float)k / (float)steps);
                    float terrain_height = dem[(long)iy * dem_width + ix];

                    if (terrain_height > los_height) {
                        is_visible = 0;
                        break;
                    }
                }

                if (is_visible) {
                    total_visible_cells++;
                }
            }
        }
    }
}

void cleanup() {
    if (dem) free(dem);
    if (observers) free(observers);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", total_visible_cells);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
