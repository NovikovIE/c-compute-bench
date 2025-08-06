#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Specific benchmark parameters
static int grid_dim_x;
static int grid_dim_y;
static int num_time_steps;

// Data arrays
static double *Ez; // Electric field z-component
static double *Hx; // Magnetic field x-component
static double *Hy; // Magnetic field y-component

// Result accumulator
static double final_result = 0.0;


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

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <grid_dim_x> <grid_dim_y> <num_time_steps> <seed>\n", argv[0]);
        exit(1);
    }

    grid_dim_x = atoi(argv[1]);
    grid_dim_y = atoi(argv[2]);
    num_time_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    // Allocate memory for the 2D grid fields
    size_t grid_size = (size_t)grid_dim_x * grid_dim_y;
    Ez = (double *)malloc(grid_size * sizeof(double));
    Hx = (double *)malloc(grid_size * sizeof(double));
    Hy = (double *)malloc(grid_size * sizeof(double));

    if (!Ez || !Hx || !Hy) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize fields to zero. A real simulation would have more complex initial conditions.
    for (size_t i = 0; i < grid_size; i++) {
        Ez[i] = 0.0;
        Hx[i] = 0.0;
        Hy[i] = 0.0;
    }
}

void run_computation() {
    // Simulation coefficients (simplified for benchmark)
    const double Ce = 0.5;
    const double Ch = 0.5;

    for (int t = 0; t < num_time_steps; t++) {
        // Update magnetic field components (Hx, Hy)
        // Loops are bounded to avoid out-of-bounds access at edges.
        for (int i = 0; i < grid_dim_x; i++) {
            for (int j = 0; j < grid_dim_y - 1; j++) {
                Hx[i * grid_dim_y + j] -= Ch * (Ez[i * grid_dim_y + (j + 1)] - Ez[i * grid_dim_y + j]);
            }
        }

        for (int i = 0; i < grid_dim_x - 1; i++) {
            for (int j = 0; j < grid_dim_y; j++) {
                Hy[i * grid_dim_y + j] += Ch * (Ez[(i + 1) * grid_dim_y + j] - Ez[i * grid_dim_y + j]);
            }
        }


        // Inject a source pulse to generate a wave. A simple sinusoidal source.
        Ez[(grid_dim_x / 2) * grid_dim_y + (grid_dim_y / 2)] += sin(0.1 * t);


        // Update electric field component (Ez)
        // Loops start from 1 to implement a simple boundary condition.
        for (int i = 1; i < grid_dim_x - 1; i++) {
            for (int j = 1; j < grid_dim_y - 1; j++) {
                Ez[i * grid_dim_y + j] += Ce * 
                    ((Hy[i * grid_dim_y + j] - Hy[(i - 1) * grid_dim_y + j]) - 
                     (Hx[i * grid_dim_y + j] - Hx[i * grid_dim_y + (j - 1)]));
            }
        }
    }

    // Accumulate a result to prevent dead code elimination
    double sum = 0.0;
    for (size_t i = 0; i < (size_t)grid_dim_x * grid_dim_y; i++) {
        sum += Ez[i];
    }
    final_result = sum;
}

void cleanup() {
    free(Ez);
    free(Hx);
    free(Hy);
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
    printf("%f\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
