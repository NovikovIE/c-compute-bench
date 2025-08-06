#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Benchmark parameters
int NUM_PARTICLES;
int NUM_TIME_STEPS;

// Particle data arrays (Global)
double* pos_x; 
double* pos_y;
double* pos_z;
double* vel_x;
double* vel_y;
double* vel_z;
double* acc_x;
double* acc_y;
double* acc_z;

// Final result accumulator
double final_result = 0.0;

// SPH simulation constants
const double G = -9.8;               // Gravity on Y axis
const double H = 0.05;               // Smoothing radius
const double H2 = H * H;             // Smoothing radius squared
const double MASS = 0.1;             // Particle mass
const double DT = 0.01;              // Time step
const double BOUND_DAMPING = -0.5;   // Wall collision damping

// Simulation domain boundaries
const double BOX_X = 1.0;
const double BOX_Y = 1.0;
const double BOX_Z = 1.0;

double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_particles> <num_time_steps> <seed>\n", argv[0]);
        exit(1);
    }
    NUM_PARTICLES = atoi(argv[1]);
    NUM_TIME_STEPS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_PARTICLES <= 0 || NUM_TIME_STEPS <= 0) {
        fprintf(stderr, "FATAL: Number of particles and time steps must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    pos_x = (double*)malloc(NUM_PARTICLES * sizeof(double));
    pos_y = (double*)malloc(NUM_PARTICLES * sizeof(double));
    pos_z = (double*)malloc(NUM_PARTICLES * sizeof(double));
    vel_x = (double*)malloc(NUM_PARTICLES * sizeof(double));
    vel_y = (double*)malloc(NUM_PARTICLES * sizeof(double));
    vel_z = (double*)malloc(NUM_PARTICLES * sizeof(double));
    acc_x = (double*)malloc(NUM_PARTICLES * sizeof(double));
    acc_y = (double*)malloc(NUM_PARTICLES * sizeof(double));
    acc_z = (double*)malloc(NUM_PARTICLES * sizeof(double));

    if (!pos_x || !pos_y || !pos_z || !vel_x || !vel_y || !vel_z || !acc_x || !acc_y || !acc_z) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
        // Place particles in a smaller cube in the top-left corner to see them drop
        pos_x[i] = rand_double() * 0.5 * BOX_X;
        pos_y[i] = 0.5 * BOX_Y + rand_double() * 0.5 * BOX_Y;
        pos_z[i] = rand_double() * 0.5 * BOX_Z;
        
        vel_x[i] = 0.0;
        vel_y[i] = 0.0;
        vel_z[i] = 0.0;
    }
}

void run_computation() {
    for (int t = 0; t < NUM_TIME_STEPS; t++) {
        // 1. Compute accelerations from pairwise particle interactions (O(N^2))
        for (int i = 0; i < NUM_PARTICLES; i++) {
            // Reset acceleration and apply gravity
            acc_x[i] = 0.0;
            acc_y[i] = G;
            acc_z[i] = 0.0;
            
            for (int j = 0; j < NUM_PARTICLES; j++) {
                if (i == j) continue;

                double dx = pos_x[i] - pos_x[j];
                double dy = pos_y[i] - pos_y[j];
                double dz = pos_z[i] - pos_z[j];
                double r2 = dx*dx + dy*dy + dz*dz;
                
                // If particles are within smoothing radius, apply repulsion force
                if (r2 < H2 && r2 > 1e-12) { 
                    double r = sqrt(r2);
                    // A simplified pressure-like force pushing particles apart
                    double force_scale = MASS * (H - r) / r;
                    
                    acc_x[i] += force_scale * dx;
                    acc_y[i] += force_scale * dy;
                    acc_z[i] += force_scale * dz;
                }
            }
        }

        // 2. Integrate to update velocities and positions
        for (int i = 0; i < NUM_PARTICLES; i++) {
            vel_x[i] += acc_x[i] * DT;
            vel_y[i] += acc_y[i] * DT;
            vel_z[i] += acc_z[i] * DT;

            pos_x[i] += vel_x[i] * DT;
            pos_y[i] += vel_y[i] * DT;
            pos_z[i] += vel_z[i] * DT;

            // 3. Handle boundary collisions
            if (pos_x[i] < 0.0)   { pos_x[i] = 0.0;   vel_x[i] *= BOUND_DAMPING; }
            if (pos_x[i] > BOX_X) { pos_x[i] = BOX_X; vel_x[i] *= BOUND_DAMPING; }
            if (pos_y[i] < 0.0)   { pos_y[i] = 0.0;   vel_y[i] *= BOUND_DAMPING; }
            if (pos_y[i] > BOX_Y) { pos_y[i] = BOX_Y; vel_y[i] *= BOUND_DAMPING; }
            if (pos_z[i] < 0.0)   { pos_z[i] = 0.0;   vel_z[i] *= BOUND_DAMPING; }
            if (pos_z[i] > BOX_Z) { pos_z[i] = BOX_Z; vel_z[i] *= BOUND_DAMPING; }
        }
    }

    // 4. Calculate a final checksum to prevent dead code elimination
    double checksum = 0.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        checksum += pos_x[i] + pos_y[i] + pos_z[i];
    }
    final_result = checksum;
}

void cleanup() {
    free(pos_x);
    free(pos_y);
    free(pos_z);
    free(vel_x);
    free(vel_y);
    free(vel_z);
    free(acc_x);
    free(acc_y);
    free(acc_z);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final checksum to stdout
    printf("%f\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
