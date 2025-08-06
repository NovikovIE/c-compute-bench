#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator --- DO NOT MODIFY ---
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

// Benchmark parameters
static int num_particles;
static int num_time_steps;
static double box_size;
static double cutoff_radius;
static double time_step;

// Particle data arrays (Global pointers)
static double *px, *py, *pz; // Positions
static double *vx, *vy, *vz; // Velocities
static double *fx, *fy, *fz; // Forces

// Final result accumulator
static double total_potential_energy;

// Helper to generate a random double in a range
double rand_double(double min, double max) {
    return min + (double)mt_rand() / (double)UINT32_MAX * (max - min);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_particles num_time_steps box_size cutoff_radius time_step seed\n", argv[0]);
        exit(1);
    }

    num_particles = atoi(argv[1]);
    num_time_steps = atoi(argv[2]);
    box_size = atof(argv[3]);
    cutoff_radius = atof(argv[4]);
    time_step = atof(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);

    mt_seed(seed);

    // Allocate memory on the heap
    size_t array_size = num_particles * sizeof(double);
    px = (double *)malloc(array_size);
    py = (double *)malloc(array_size);
    pz = (double *)malloc(array_size);
    vx = (double *)malloc(array_size);
    vy = (double *)malloc(array_size);
    vz = (double *)malloc(array_size);
    fx = (double *)malloc(array_size);
    fy = (double *)malloc(array_size);
    fz = (double *)malloc(array_size);

    if (!px || !py || !pz || !vx || !vy || !vz || !fx || !fy || !fz) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize particle positions randomly within the box
    for (int i = 0; i < num_particles; ++i) {
        px[i] = rand_double(0.0, box_size);
        py[i] = rand_double(0.0, box_size);
        pz[i] = rand_double(0.0, box_size);
    }

    // Initialize velocities and forces to zero
    memset(vx, 0, array_size);
    memset(vy, 0, array_size);
    memset(vz, 0, array_size);
}

void run_computation() {
    total_potential_energy = 0.0;
    double cutoff_sq = cutoff_radius * cutoff_radius;
    double box_half = box_size / 2.0;

    for (int step = 0; step < num_time_steps; ++step) {
        // Reset forces and step potential energy
        size_t array_size = num_particles * sizeof(double);
        memset(fx, 0, array_size);
        memset(fy, 0, array_size);
        memset(fz, 0, array_size);
        double step_potential_energy = 0.0;

        // 1. Calculate forces and potential energy
        for (int i = 0; i < num_particles; ++i) {
            for (int j = i + 1; j < num_particles; ++j) {
                double dx = px[i] - px[j];
                double dy = py[i] - py[j];
                double dz = pz[i] - pz[j];

                // Periodic Boundary Conditions (Minimum Image Convention)
                if (dx > box_half) dx -= box_size;
                else if (dx < -box_half) dx += box_size;
                if (dy > box_half) dy -= box_size;
                else if (dy < -box_half) dy += box_size;
                if (dz > box_half) dz -= box_size;
                else if (dz < -box_half) dz += box_size;

                double r2 = dx * dx + dy * dy + dz * dz;

                if (r2 < cutoff_sq && r2 > 1e-12) {
                    double r2_inv = 1.0 / r2;
                    double r6_inv = r2_inv * r2_inv * r2_inv;
                    double lj_potential_term = r6_inv * (r6_inv - 1.0);
                    double lj_force_term = r2_inv * lj_potential_term;

                    step_potential_energy += 4.0 * lj_potential_term;

                    double force_scalar = 48.0 * lj_force_term;
                    fx[i] += force_scalar * dx;
                    fy[i] += force_scalar * dy;
                    fz[i] += force_scalar * dz;

                    fx[j] -= force_scalar * dx;
                    fy[j] -= force_scalar * dy;
                    fz[j] -= force_scalar * dz;
                }
            }
        }

        // 2. Update particle positions and velocities (Velocity Verlet integration, mass=1)
        for (int i = 0; i < num_particles; ++i) {
            vx[i] += fx[i] * time_step;
            vy[i] += fy[i] * time_step;
            vz[i] += fz[i] * time_step;

            px[i] += vx[i] * time_step;
            py[i] += vy[i] * time_step;
            pz[i] += vz[i] * time_step;

            // Apply periodic boundary conditions (Wrap around)
            if (px[i] < 0.0) px[i] += box_size;
            if (px[i] >= box_size) px[i] -= box_size;
            if (py[i] < 0.0) py[i] += box_size;
            if (py[i] >= box_size) py[i] -= box_size;
            if (pz[i] < 0.0) pz[i] += box_size;
            if (pz[i] >= box_size) pz[i] -= box_size;
        }
        total_potential_energy += step_potential_energy;
    }
}

void cleanup() {
    free(px); free(py); free(pz);
    free(vx); free(vy); free(vz);
    free(fx); free(fy); free(fz);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated potential energy to stdout
    printf("%f\n", total_potential_energy);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
