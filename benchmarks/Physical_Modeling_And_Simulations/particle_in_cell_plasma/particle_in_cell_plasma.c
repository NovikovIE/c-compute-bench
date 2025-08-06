#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <math.h>

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

// --- Benchmark Data and Globals ---
typedef struct {
    int num_particles;
    int grid_size;
    int num_time_steps;

    double* particle_pos;
    double* particle_vel;
    double* grid_charge_density;
    double* grid_field;

    double final_result;
} BenchmarkData;

static BenchmarkData gdata;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_particles grid_size num_time_steps seed\n", argv[0]);
        exit(1);
    }

    gdata.num_particles = atoi(argv[1]);
    gdata.grid_size = atoi(argv[2]);
    gdata.num_time_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    gdata.particle_pos = (double*)malloc(gdata.num_particles * sizeof(double));
    gdata.particle_vel = (double*)malloc(gdata.num_particles * sizeof(double));
    gdata.grid_charge_density = (double*)malloc(gdata.grid_size * sizeof(double));
    gdata.grid_field = (double*)malloc(gdata.grid_size * sizeof(double));

    if (!gdata.particle_pos || !gdata.particle_vel || !gdata.grid_charge_density || !gdata.grid_field) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    double grid_domain = (double)gdata.grid_size;
    for (int i = 0; i < gdata.num_particles; ++i) {
        gdata.particle_pos[i] = (mt_rand() / (double)UINT32_MAX) * grid_domain;
        // Initialize with small random velocities
        gdata.particle_vel[i] = (mt_rand() / (double)UINT32_MAX - 0.5) * 0.1;
    }
}

void run_computation() {
    const double dt = 0.1; // Time step size
    const double q = -1.0; // Charge of particles
    const double m = 1.0;  // Mass of particles
    const double q_div_m = q / m;
    const double grid_domain = (double)gdata.grid_size;

    for (int t = 0; t < gdata.num_time_steps; ++t) {
        // 1. Reset grid charge density
        memset(gdata.grid_charge_density, 0, gdata.grid_size * sizeof(double));

        // 2. Deposit charge from particles to grid (Cloud-in-Cell)
        for (int i = 0; i < gdata.num_particles; ++i) {
            double pos = gdata.particle_pos[i];
            int grid_idx = (int)pos;
            double frac = pos - grid_idx;

            // Apply charge to two nearest grid points, with periodic boundaries
            int next_idx = (grid_idx + 1) % gdata.grid_size;
            gdata.grid_charge_density[grid_idx] += (1.0 - frac) * q;
            gdata.grid_charge_density[next_idx] += frac * q;
        }

        // 3. Solve for electric field from charge density
        // Using a simple 1D solver: integrate charge density (dE/dx = rho)
        // Assume periodic system, so net field from all charges at a point
        // To simplify for benchmark, use a simple integration with a neutralizing background
        double total_charge = 0.0;
        for (int i = 0; i < gdata.grid_size; ++i) {
            total_charge += gdata.grid_charge_density[i];
        }
        double avg_charge = total_charge / gdata.grid_size;
        for (int i = 0; i < gdata.grid_size; ++i) {
            gdata.grid_charge_density[i] -= avg_charge;
        }
        
        gdata.grid_field[0] = 0.0;
        for (int i = 1; i < gdata.grid_size; ++i) {
            gdata.grid_field[i] = gdata.grid_field[i-1] + gdata.grid_charge_density[i-1]; // Assuming dx=1
        }

        // 4. Interpolate field to particles and push them
        for (int i = 0; i < gdata.num_particles; ++i) {
            double pos = gdata.particle_pos[i];
            int grid_idx = (int)pos;
            double frac = pos - grid_idx;

            // Interpolate E-field from grid to particle position
            int next_idx = (grid_idx + 1) % gdata.grid_size;
            double e_field = gdata.grid_field[grid_idx] * (1.0 - frac) + gdata.grid_field[next_idx] * frac;

            // Update velocity (v_new = v_old + (q*E/m) * dt)
            gdata.particle_vel[i] += q_div_m * e_field * dt;

            // Update position (x_new = x_old + v_new * dt)
            gdata.particle_pos[i] += gdata.particle_vel[i] * dt;

            // Apply periodic boundary conditions
            if (gdata.particle_pos[i] >= grid_domain) {
                gdata.particle_pos[i] -= grid_domain;
            } else if (gdata.particle_pos[i] < 0.0) {
                gdata.particle_pos[i] += grid_domain;
            }
        }
    }

    // Calculate final checksum to prevent dead code elimination
    double checksum = 0.0;
    for (int i = 0; i < gdata.num_particles; ++i) {
        checksum += gdata.particle_pos[i];
    }
    gdata.final_result = checksum;
}

void cleanup() {
    free(gdata.particle_pos);
    free(gdata.particle_vel);
    free(gdata.grid_charge_density);
    free(gdata.grid_field);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%.2f\n", gdata.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}