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

// Benchmark parameters and data structures
int num_atoms;
int num_time_steps;
double non_bonded_cutoff_radius;
double box_size;
double *positions_x, *positions_y, *positions_z;
double *velocities_x, *velocities_y, *velocities_z;
double *forces_x, *forces_y, *forces_z;
double total_kinetic_energy; // The final result to print

// A psuedo-random number generator for doubles in [0, 1)
double random_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_atoms> <num_time_steps> <non_bonded_cutoff_radius> <seed>\n", argv[0]);
        exit(1);
    }

    num_atoms = atoi(argv[1]);
    num_time_steps = atoi(argv[2]);
    non_bonded_cutoff_radius = atof(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    mt_seed(seed);

    // Assume a constant reduced density (rho*) of 0.8
    double density = 0.8;
    double volume = (double)num_atoms / density;
    box_size = cbrt(volume);

    // Allocate memory
    positions_x = (double*)malloc(num_atoms * sizeof(double));
    positions_y = (double*)malloc(num_atoms * sizeof(double));
    positions_z = (double*)malloc(num_atoms * sizeof(double));
    velocities_x = (double*)malloc(num_atoms * sizeof(double));
    velocities_y = (double*)malloc(num_atoms * sizeof(double));
    velocities_z = (double*)malloc(num_atoms * sizeof(double));
    forces_x = (double*)malloc(num_atoms * sizeof(double));
    forces_y = (double*)malloc(num_atoms * sizeof(double));
    forces_z = (double*)malloc(num_atoms * sizeof(double));

    if (!positions_x || !positions_y || !positions_z ||
        !velocities_x || !velocities_y || !velocities_z ||
        !forces_x || !forces_y || !forces_z) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize atom positions and velocities
    for (int i = 0; i < num_atoms; ++i) {
        positions_x[i] = random_double() * box_size;
        positions_y[i] = random_double() * box_size;
        positions_z[i] = random_double() * box_size;
        // Initial velocities from -0.5 to 0.5
        velocities_x[i] = random_double() - 0.5;
        velocities_y[i] = random_double() - 0.5;
        velocities_z[i] = random_double() - 0.5;
    }

    total_kinetic_energy = 0.0;
}

void run_computation() {
    const double dt = 0.001; // Time step
    const double mass = 1.0; // Mass of each atom
    const double epsilon = 1.0; // Lennard-Jones potential depth
    const double sigma = 1.0; // Lennard-Jones finite distance where potential is zero
    const double cutoff_sq = non_bonded_cutoff_radius * non_bonded_cutoff_radius;

    for (int step = 0; step < num_time_steps; ++step) {
        // 1. Calculate forces (Lennard-Jones)
        memset(forces_x, 0, num_atoms * sizeof(double));
        memset(forces_y, 0, num_atoms * sizeof(double));
        memset(forces_z, 0, num_atoms * sizeof(double));

        for (int i = 0; i < num_atoms; ++i) {
            for (int j = i + 1; j < num_atoms; ++j) {
                double dx = positions_x[i] - positions_x[j];
                double dy = positions_y[i] - positions_y[j];
                double dz = positions_z[i] - positions_z[j];

                // Apply periodic boundary conditions (minimum image convention)
                dx -= box_size * rint(dx / box_size);
                dy -= box_size * rint(dy / box_size);
                dz -= box_size * rint(dz / box_size);

                double r_sq = dx * dx + dy * dy + dz * dz;

                if (r_sq < cutoff_sq) {
                    // Optimized LJ force calculation, avoids sqrt()
                    double r2i = 1.0 / r_sq;
                    double r6i = r2i * r2i * r2i;
                    double force_term = 48.0 * epsilon * r2i * r6i * (sigma*sigma*sigma*sigma*sigma*sigma * r6i - 0.5);
                    
                    forces_x[i] += force_term * dx;
                    forces_y[i] += force_term * dy;
                    forces_z[i] += force_term * dz;
                    forces_x[j] -= force_term * dx;
                    forces_y[j] -= force_term * dy;
                    forces_z[j] -= force_term * dz;
                }
            }
        }

        double current_ke = 0.0;
        // 2. Update velocities and positions (using Euler-Cromer integrator)
        for (int i = 0; i < num_atoms; ++i) {
            // Update velocities
            velocities_x[i] += (forces_x[i] / mass) * dt;
            velocities_y[i] += (forces_y[i] / mass) * dt;
            velocities_z[i] += (forces_z[i] / mass) * dt;

            // Update positions
            positions_x[i] += velocities_x[i] * dt;
            positions_y[i] += velocities_y[i] * dt;
            positions_z[i] += velocities_z[i] * dt;

            // Apply periodic boundary conditions to positions
            if (positions_x[i] < 0.0) positions_x[i] += box_size;
            if (positions_x[i] >= box_size) positions_x[i] -= box_size;
            if (positions_y[i] < 0.0) positions_y[i] += box_size;
            if (positions_y[i] >= box_size) positions_y[i] -= box_size;
            if (positions_z[i] < 0.0) positions_z[i] += box_size;
            if (positions_z[i] >= box_size) positions_z[i] -= box_size;

            // 3. Accumulate result (prevent dead code elimination)
            current_ke += 0.5 * mass * (velocities_x[i] * velocities_x[i] + 
                                        velocities_y[i] * velocities_y[i] + 
                                        velocities_z[i] * velocities_z[i]);
        }
        total_kinetic_energy += current_ke; 
    }
}

void cleanup() {
    free(positions_x);
    free(positions_y);
    free(positions_z);
    free(velocities_x);
    free(velocities_y);
    free(velocities_z);
    free(forces_x);
    free(forces_y);
    free(forces_z);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated kinetic energy to stdout
    printf("%f\n", total_kinetic_energy);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
