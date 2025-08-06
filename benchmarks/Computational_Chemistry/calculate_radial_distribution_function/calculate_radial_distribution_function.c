#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Mersenne Twister (MT19937)
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

// Structure to hold all benchmark data
typedef struct {
    int num_frames;
    int num_particles;
    int num_bins;
    double box_size;
    
    // Flattened array for particle positions: layout is [frame][particle][coord]
    // Total size: num_frames * num_particles * 3
    double* particle_positions;
    
    // Histogram for the Radial Distribution Function
    double* rdf_histogram;
    
    // Final result to prevent dead code elimination
    double final_result;
} BenchmarkData;

// Global instance of the benchmark data
static BenchmarkData g_data;

// Function to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_frames> <num_particles> <num_bins> <box_size> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_frames = atoi(argv[1]);
    g_data.num_particles = atoi(argv[2]);
    g_data.num_bins = atoi(argv[3]);
    g_data.box_size = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    size_t pos_count = (size_t)g_data.num_frames * g_data.num_particles * 3;
    g_data.particle_positions = (double*)malloc(pos_count * sizeof(double));
    if (!g_data.particle_positions) {
        fprintf(stderr, "Failed to allocate memory for particle_positions.\n");
        exit(1);
    }

    g_data.rdf_histogram = (double*)malloc(g_data.num_bins * sizeof(double));
    if (!g_data.rdf_histogram) {
        fprintf(stderr, "Failed to allocate memory for rdf_histogram.\n");
        free(g_data.particle_positions);
        exit(1);
    }

    // Initialize particle positions randomly within the box
    for (size_t i = 0; i < pos_count; ++i) {
        g_data.particle_positions[i] = rand_double() * g_data.box_size;
    }

    // Initialize RDF histogram to zero
    for (int i = 0; i < g_data.num_bins; ++i) {
        g_data.rdf_histogram[i] = 0.0;
    }
    
    g_data.final_result = 0.0;
}

void run_computation() {
    double max_r = g_data.box_size / 2.0;
    double max_r_sq = max_r * max_r;
    double delta_r = max_r / g_data.num_bins;
    double total_count = 0.0;

    for (int f = 0; f < g_data.num_frames; ++f) {
        const double* frame_positions = g_data.particle_positions + (size_t)f * g_data.num_particles * 3;
        
        for (int i = 0; i < g_data.num_particles; ++i) {
            double ix = frame_positions[i * 3 + 0];
            double iy = frame_positions[i * 3 + 1];
            double iz = frame_positions[i * 3 + 2];

            for (int j = i + 1; j < g_data.num_particles; ++j) {
                double jx = frame_positions[j * 3 + 0];
                double jy = frame_positions[j * 3 + 1];
                double jz = frame_positions[j * 3 + 2];
                
                // Minimum Image Convention for periodic boundary conditions
                double dx = ix - jx;
                double dy = iy - jy;
                double dz = iz - jz;
                
                dx -= g_data.box_size * round(dx / g_data.box_size);
                dy -= g_data.box_size * round(dy / g_data.box_size);
                dz -= g_data.box_size * round(dz / g_data.box_size);

                double r_sq = dx * dx + dy * dy + dz * dz;

                if (r_sq < max_r_sq) {
                    double r = sqrt(r_sq);
                    int bin_index = (int)(r / delta_r);
                    if (bin_index < g_data.num_bins) {
                        // Each pair contributes to the count. Adding 2 is standard (i-j and j-i).
                        g_data.rdf_histogram[bin_index] += 2.0;
                    }
                }
            }
        }
    }

    // Calculate a final result to ensure computation is not optimized away
    for (int i = 0; i < g_data.num_bins; ++i) {
        total_count += g_data.rdf_histogram[i];
    }
    g_data.final_result = total_count;
}

void cleanup() {
    free(g_data.particle_positions);
    free(g_data.rdf_histogram);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%.2f\n", g_data.final_result);
    
    cleanup();

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
