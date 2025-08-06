#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// Represents a single point vortex in 2D space
typedef struct {
    double x;      // x-coordinate
    double y;      // y-coordinate
    double gamma;  // circulation strength
} Vortex;

// Global struct to hold all benchmark data
static struct {
    int num_vortices;
    int num_time_steps;
    double dt;
    double delta;

    Vortex *vortices;
    double *vel_x; // Buffer for x-velocities
    double *vel_y; // Buffer for y-velocities

    double final_result;
} G;

// Generates a random double in [min, max]
double rand_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_vortices> <num_time_steps> <dt> <delta> <seed>\n", argv[0]);
        exit(1);
    }

    G.num_vortices = atoi(argv[1]);
    G.num_time_steps = atoi(argv[2]);
    G.dt = atof(argv[3]);
    G.delta = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    G.vortices = (Vortex *)malloc(G.num_vortices * sizeof(Vortex));
    G.vel_x = (double *)malloc(G.num_vortices * sizeof(double));
    G.vel_y = (double *)malloc(G.num_vortices * sizeof(double));

    if (!G.vortices || !G.vel_x || !G.vel_y) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < G.num_vortices; ++i) {
        G.vortices[i].x = rand_double(-1.0, 1.0);
        G.vortices[i].y = rand_double(-1.0, 1.0);
        G.vortices[i].gamma = rand_double(-0.1, 0.1);
    }
    G.final_result = 0.0;
}

void run_computation() {
    const double delta_sq = G.delta * G.delta;
    const double inv_2pi = 1.0 / (2.0 * M_PI);

    for (int t = 0; t < G.num_time_steps; ++t) {
        // Step 1: Calculate induced velocities for all vortices (N^2 interaction).
        for (int i = 0; i < G.num_vortices; ++i) {
            double total_ux = 0.0;
            double total_uy = 0.0;

            for (int j = 0; j < G.num_vortices; ++j) {
                if (i == j) continue;

                double dx = G.vortices[i].x - G.vortices[j].x;
                double dy = G.vortices[i].y - G.vortices[j].y;
                
                // Regularized squared distance to avoid singularity
                double r_sq_reg = dx * dx + dy * dy + delta_sq;
                double r_sq_inv = 1.0 / r_sq_reg;

                // Biot-Savart law for 2D point vortices: u = -gamma * dy / (2*pi*r^2), v = gamma * dx / (2*pi*r^2)
                double velocity_coeff = G.vortices[j].gamma * inv_2pi * r_sq_inv;
                total_ux -= velocity_coeff * dy;
                total_uy += velocity_coeff * dx;
            }
            G.vel_x[i] = total_ux;
            G.vel_y[i] = total_uy;
        }

        // Step 2: Update all vortex positions using the calculated velocities (Euler step).
        for (int i = 0; i < G.num_vortices; ++i) {
            G.vortices[i].x += G.vel_x[i] * G.dt;
            G.vortices[i].y += G.vel_y[i] * G.dt;
        }
    }

    // Accumulate a final result to prevent dead code elimination.
    double total_sum = 0.0;
    for (int i = 0; i < G.num_vortices; ++i) {
        total_sum += G.vortices[i].x + G.vortices[i].y;
    }
    G.final_result = total_sum;
}

void cleanup() {
    free(G.vortices);
    free(G.vel_x);
    free(G.vel_y);
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
    printf("%f\n", G.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
