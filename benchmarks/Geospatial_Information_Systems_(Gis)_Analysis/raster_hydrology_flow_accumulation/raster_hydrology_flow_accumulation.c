#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

/************************************************************************
 *                  MERSENNE TWISTER (from specification)               *
 ************************************************************************/

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

/************************************************************************
 *                        BENCHMARK SPECIFIC CODE                       *
 ************************************************************************/
 
// Benchmark: Raster Hydrology Flow Accumulation (D8 method)
// Description: This benchmark simulates a D8 flow accumulation algorithm. It first
//              creates a random Digital Elevation Model (DEM). It then determines a 
//              processing order by sorting cells from highest to lowest elevation. 
//              The core computation iterates through the sorted cells, finds the 
//              neighbor with the steepest descent (the D8 flow direction), and adds
//              the current cell's accumulated flow to that neighbor. This is a common
//              and computationally intensive task in hydrological analysis.

// Global data structure to hold all benchmark data
static struct {
    int dem_width;
    int dem_height;
    long long num_cells;
    float *dem;                 // Digital Elevation Model (raster)
    long long *flow_accumulation; // Accumulated flow for each cell
    int *processing_order;      // Cell indices sorted by elevation (descending)
    long long result_sum;       // Final result to prevent dead code elimination
} g_data;

// Temporary struct for sorting cells by elevation during setup
typedef struct {
    int index;
    float elevation;
} Cell;

// qsort comparison function for sorting cells by elevation in descending order
int compare_cells(const void* a, const void* b) {
    Cell* cell_a = (Cell*)a;
    Cell* cell_b = (Cell*)b;
    if (cell_a->elevation < cell_b->elevation) return 1;
    if (cell_a->elevation > cell_b->elevation) return -1;
    return 0;
}

void setup_benchmark(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <dem_width> <dem_height> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.dem_width = atoi(argv[1]);
    g_data.dem_height = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    mt_seed(seed);

    if (g_data.dem_width <= 0 || g_data.dem_height <= 0) {
        fprintf(stderr, "FATAL: DEM dimensions must be positive.\n");
        exit(1);
    }

    g_data.num_cells = (long long)g_data.dem_width * g_data.dem_height;

    // Allocate memory
    g_data.dem = (float*)malloc(g_data.num_cells * sizeof(float));
    g_data.flow_accumulation = (long long*)malloc(g_data.num_cells * sizeof(long long));
    g_data.processing_order = (int*)malloc(g_data.num_cells * sizeof(int));
    Cell* temp_cells = (Cell*)malloc(g_data.num_cells * sizeof(Cell));

    if (!g_data.dem || !g_data.flow_accumulation || !g_data.processing_order || !temp_cells) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize DEM with random elevations and prepare for sorting
    for (long long i = 0; i < g_data.num_cells; ++i) {
        g_data.dem[i] = (float)mt_rand() / (float)UINT32_MAX * 1000.0f;
        temp_cells[i].index = i;
        temp_cells[i].elevation = g_data.dem[i];
        // Each cell starts with one unit of flow (itself)
        g_data.flow_accumulation[i] = 1;
    }

    // Sort cells by elevation (descending) to get processing order
    qsort(temp_cells, g_data.num_cells, sizeof(Cell), compare_cells);

    // Store the sorted indices
    for (long long i = 0; i < g_data.num_cells; ++i) {
        g_data.processing_order[i] = temp_cells[i].index;
    }

    free(temp_cells); // Free temporary sorting structure
}

void run_computation() {
    // D8 neighbor directions: (dx, dy) pairs relative to current cell
    // E, NE, N, NW, W, SW, S, SE
    const int dx[] = {1,  1,  0, -1, -1, -1, 0, 1};
    const int dy[] = {1,  0, -1, -1,  0,  1, 1, 1}; // In image coordinates (Y points down)
    const float inv_dist[] = {1.0f / 1.41421356f, 1.0f, 1.0f / 1.41421356f, 1.0f, 1.0f / 1.41421356f, 1.0f, 1.0f / 1.41421356f, 1.0f};

    // Process cells from highest to lowest elevation
    for (long long i = 0; i < g_data.num_cells; ++i) {
        int current_idx = g_data.processing_order[i];
        int current_x = current_idx % g_data.dem_width;
        int current_y = current_idx / g_data.dem_width;
        float current_elev = g_data.dem[current_idx];

        float max_slope = -1.0f;
        int sink_idx = -1;

        // Find neighbor with the steepest descent
        for (int d = 0; d < 8; ++d) {
            int nx = current_x + dx[d];
            int ny = current_y + dy[d];

            // Check boundary conditions
            if (nx >= 0 && nx < g_data.dem_width && ny >= 0 && ny < g_data.dem_height) {
                int neighbor_idx = ny * g_data.dem_width + nx;
                float neighbor_elev = g_data.dem[neighbor_idx];

                if (neighbor_elev < current_elev) {
                    float slope = (current_elev - neighbor_elev) * inv_dist[d]; // slope = delta_h / distance
                    if (slope > max_slope) {
                        max_slope = slope;
                        sink_idx = neighbor_idx;
                    }
                }
            }
        }

        // If a sink was found, accumulate flow
        if (sink_idx != -1) {
            g_data.flow_accumulation[sink_idx] += g_data.flow_accumulation[current_idx];
        }
    }

    // Calculate final checksum to prevent dead code elimination
    long long sum = 0;
    for (long long i = 0; i < g_data.num_cells; ++i) {
        sum += g_data.flow_accumulation[i];
    }
    g_data.result_sum = sum;
}

void cleanup() {
    free(g_data.dem);
    free(g_data.flow_accumulation);
    free(g_data.processing_order);
    g_data.dem = NULL;
    g_data.flow_accumulation = NULL;
    g_data.processing_order = NULL;
}

int main(int argc, char **argv) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print checksum to stdout
    printf("%lld\n", g_data.result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}