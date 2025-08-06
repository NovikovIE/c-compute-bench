#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// Mersenne Twister (Do Not Modify - Include This Verbatim)
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
// End of Mersenne Twister

// --- Benchmark Data and Parameters ---

// Parameters
int WIDTH;
int HEIGHT;
int NUM_ZONES;
long long TOTAL_CELLS;

// Data grids
float *value_raster; // Represents a continuous value field (e.g., rainfall, elevation)
int *zone_raster;    // Defines the zone for each cell

// Zonal statistics arrays
double *zone_sums;     // Sum of values for each zone
long long *zone_counts; // Count of cells for each zone

// Final result accumulator to prevent dead-code elimination
double final_result_accumulator;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <width> <height> <num_zones> <seed>\n", argv[0]);
        exit(1);
    }

    WIDTH = atoi(argv[1]);
    HEIGHT = atoi(argv[2]);
    NUM_ZONES = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if(WIDTH <= 0 || HEIGHT <= 0 || NUM_ZONES <= 0) {
        fprintf(stderr, "Error: width, height, and num_zones must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    TOTAL_CELLS = (long long)WIDTH * HEIGHT;

    // Allocate memory for rasters
    value_raster = (float *)malloc(TOTAL_CELLS * sizeof(float));
    zone_raster = (int *)malloc(TOTAL_CELLS * sizeof(int));

    // Allocate memory for zonal statistics results
    zone_sums = (double *)malloc(NUM_ZONES * sizeof(double));
    zone_counts = (long long *)malloc(NUM_ZONES * sizeof(long long));

    if (!value_raster || !zone_raster || !zone_sums || !zone_counts) {
        fprintf(stderr, "Fatal: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize rasters with random data
    for (long long i = 0; i < TOTAL_CELLS; ++i) {
        // Value raster: random float between 0.0 and 1000.0
        value_raster[i] = (float)mt_rand() / (float)UINT32_MAX * 1000.0f;
        // Zone raster: random zone ID from 0 to NUM_ZONES-1
        zone_raster[i] = mt_rand() % NUM_ZONES;
    }

    // Initialize result arrays
    memset(zone_sums, 0, NUM_ZONES * sizeof(double));
    memset(zone_counts, 0, NUM_ZONES * sizeof(long long));
    
    final_result_accumulator = 0.0;
}

void run_computation() {
    // Phase 1: Iterate through all cells and accumulate sums and counts for each zone.
    // This loop's performance is dominated by memory access patterns, both sequential (rasters)
    // and random (zone_sums/counts), which is typical for raster analysis.
    for (long long i = 0; i < TOTAL_CELLS; ++i) {
        int zone_id = zone_raster[i];
        float value = value_raster[i];
        zone_sums[zone_id] += value;
        zone_counts[zone_id]++;
    }

    // Phase 2: Calculate the mean for each zone and sum the means.
    // This second pass ensures there's some computation after the main loop
    // and provides a single value for output, preventing dead-code elimination.
    double total_mean_sum = 0.0;
    for (int i = 0; i < NUM_ZONES; ++i) {
        if (zone_counts[i] > 0) {
            double mean = zone_sums[i] / (double)zone_counts[i];
            total_mean_sum += mean;
        }
    }
    final_result_accumulator = total_mean_sum;
}

void cleanup() {
    free(value_raster);
    free(zone_raster);
    free(zone_sums);
    free(zone_counts);
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
    printf("%f\n", final_result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
