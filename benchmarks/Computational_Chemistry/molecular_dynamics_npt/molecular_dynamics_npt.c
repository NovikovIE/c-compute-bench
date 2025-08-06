#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- Mersenne Twister (verbatim) ---
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
// --- End Mersenne Twister ---


// --- Benchmark Globals ---
int num_atoms;
int num_time_steps;
double non_bonded_cutoff_radius;
double box_size;

double* positions_x;
double* positions_y;
double* positions_z;
double* velocities_x;
double* velocities_y;
double* velocities_z;
double* forces_x;
double* forces_y;
double* forces_z;
double* masses;

double final_potential_energy;


// Generate a random double between low and high
double rand_double(double low, double high) {
    return low + (mt_rand() / (double)UINT32_MAX) * (high - low);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_atoms> <num_time_steps> <non_bonded_cutoff_radius> <seed>\n", argv[0]);
        exit(1);
    }

    num_atoms = atoi(argv[1]);
    num_time_steps = atoi(argv[2]);
    non_bonded_cutoff_radius = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);
    
    box_size = cbrt(num_atoms * 100.0); // Simple heuristic for constant density

    positions_x  = (double*)malloc(num_atoms * sizeof(double));
    positions_y  = (double*)malloc(num_atoms * sizeof(double));
    positions_z  = (double*)malloc(num_atoms * sizeof(double));
    velocities_x = (double*)malloc(num_atoms * sizeof(double));
    velocities_y = (double*)malloc(num_atoms * sizeof(double));
    velocities_z = (double*)malloc(num_atoms * sizeof(double));
    forces_x     = (double*)malloc(num_atoms * sizeof(double));
    forces_y     = (double*)malloc(num_atoms * sizeof(double));
    forces_z     = (double*)malloc(num_atoms * sizeof(double));
    masses       = (double*)malloc(num_atoms * sizeof(double));

    if (!positions_x || !positions_y || !positions_z || !velocities_x || !velocities_y || !velocities_z || !forces_x || !forces_y || !forces_z || !masses) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_atoms; i++) {
        positions_x[i] = rand_double(0.0, box_size);
        positions_y[i] = rand_double(0.0, box_size);
        positions_z[i] = rand_double(0.0, box_size);
        velocities_x[i] = rand_double(-1.0, 1.0);
        velocities_y[i] = rand_double(-1.0, 1.0);
        velocities_z[i] = rand_double(-1.0, 1.0);
        masses[i] = 12.01; // Carbon atom mass
    }
}

void run_computation() {
    double dt = 0.001; // Time step
    double cutoff_radius_sq = non_bonded_cutoff_radius * non_bonded_cutoff_radius;
    double current_potential_energy = 0.0;

    for (int step = 0; step < num_time_steps; ++step) {
        // 1. Force Calculation (Lennard-Jones potential)
        current_potential_energy = 0.0;
        memset(forces_x, 0, num_atoms * sizeof(double));
        memset(forces_y, 0, num_atoms * sizeof(double));
        memset(forces_z, 0, num_atoms * sizeof(double));

        for (int i = 0; i < num_atoms; ++i) {
            for (int j = i + 1; j < num_atoms; ++j) {
                double dx = positions_x[i] - positions_x[j];
                double dy = positions_y[i] - positions_y[j];
                double dz = positions_z[i] - positions_z[j];
                
                dx -= box_size * round(dx / box_size);
                dy -= box_size * round(dy / box_size);
                dz -= box_size * round(dz / box_size);

                double r2 = dx * dx + dy * dy + dz * dz;

                if (r2 < cutoff_radius_sq && r2 > 1e-9) {
                    double r2_inv = 1.0 / r2;
                    double r6_inv = r2_inv * r2_inv * r2_inv;
                    
                    double lj_potential = 4.0 * (r6_inv * r6_inv - r6_inv);
                    current_potential_energy += lj_potential;

                    double force_magnitude = 24.0 * (2.0 * r6_inv * r6_inv - r6_inv) * r2_inv;

                    double fx = force_magnitude * dx;
                    double fy = force_magnitude * dy;
                    double fz = force_magnitude * dz;

                    forces_x[i] += fx;
                    forces_x[j] -= fx;
                    forces_y[i] += fy;
                    forces_y[j] -= fy;
                    forces_z[i] += fz;
                    forces_z[j] -= fz;
                }
            }
        }

        // 2. Integration (Simplified Velocity Verlet)
        for (int i = 0; i < num_atoms; ++i) {
            double mass_inv = 1.0 / masses[i];
            velocities_x[i] += forces_x[i] * mass_inv * dt;
            velocities_y[i] += forces_y[i] * mass_inv * dt;
            velocities_z[i] += forces_z[i] * mass_inv * dt;
            
            positions_x[i] += velocities_x[i] * dt;
            positions_y[i] += velocities_y[i] * dt;
            positions_z[i] += velocities_z[i] * dt;
        }
        
        // 3. Simple Thermostat (Velocity Rescaling)
        double current_kinetic_energy = 0.0;
        for (int i = 0; i < num_atoms; i++) {
            current_kinetic_energy += 0.5 * masses[i] * (velocities_x[i]*velocities_x[i] + velocities_y[i]*velocities_y[i] + velocities_z[i]*velocities_z[i]);
        }
        
        double target_temp = 1.0;
        double current_temp = (2.0 * current_kinetic_energy) / (3.0 * num_atoms);
        
        if (current_temp > 1e-6) {
            double scale_factor = sqrt(target_temp / current_temp);
            for (int i = 0; i < num_atoms; i++) {
                velocities_x[i] *= scale_factor;
                velocities_y[i] *= scale_factor;
                velocities_z[i] *= scale_factor;
            }
        }
    }
    
    final_potential_energy = current_potential_energy;
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
    free(masses);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_potential_energy);

    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
