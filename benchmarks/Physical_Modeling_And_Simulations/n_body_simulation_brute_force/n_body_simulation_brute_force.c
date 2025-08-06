#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (Do Not Modify) ---
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

// A utility function to generate a random double between -1.0 and 1.0
double rand_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
}

// Global struct to hold all benchmark data
typedef struct {
    int num_bodies;
    int num_time_steps;
    double time_delta;
    
    double *pos_x, *pos_y, *pos_z;
    double *vel_x, *vel_y, *vel_z;
    double *acc_x, *acc_y, *acc_z;
    double *mass;

    double final_result; // To store total kinetic energy
} BenchmarkData;

static BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_bodies num_time_steps time_delta seed\n", argv[0]);
        exit(1);
    }

    g_data.num_bodies = atoi(argv[1]);
    g_data.num_time_steps = atoi(argv[2]);
    g_data.time_delta = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory using Structure of Arrays (SoA) for better cache performance
    int n = g_data.num_bodies;
    g_data.pos_x = (double*)malloc(n * sizeof(double));
    g_data.pos_y = (double*)malloc(n * sizeof(double));
    g_data.pos_z = (double*)malloc(n * sizeof(double));
    g_data.vel_x = (double*)malloc(n * sizeof(double));
    g_data.vel_y = (double*)malloc(n * sizeof(double));
    g_data.vel_z = (double*)malloc(n * sizeof(double));
    g_data.acc_x = (double*)malloc(n * sizeof(double));
    g_data.acc_y = (double*)malloc(n * sizeof(double));
    g_data.acc_z = (double*)malloc(n * sizeof(double));
    g_data.mass  = (double*)malloc(n * sizeof(double));

    // Check for allocation failure
    if (!g_data.pos_x || !g_data.pos_y || !g_data.pos_z || !g_data.vel_x || !g_data.vel_y || !g_data.vel_z ||
        !g_data.acc_x || !g_data.acc_y || !g_data.acc_z || !g_data.mass) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize bodies with random positions, velocities, and masses
    const double mass_scale = 1000.0;
    const double pos_scale = 10.0;
    const double vel_scale = 1.0;

    for (int i = 0; i < n; ++i) {
        // Mass is positive
        g_data.mass[i] = ((double)mt_rand() / (double)UINT32_MAX) * mass_scale + 1.0;
        
        // Positions and velocities in a [-scale, +scale] cube
        g_data.pos_x[i] = rand_double() * pos_scale;
        g_data.pos_y[i] = rand_double() * pos_scale;
        g_data.pos_z[i] = rand_double() * pos_scale;
        g_data.vel_x[i] = rand_double() * vel_scale;
        g_data.vel_y[i] = rand_double() * vel_scale;
        g_data.vel_z[i] = rand_double() * vel_scale;
    }
}

void run_computation() {
    const int n = g_data.num_bodies;
    const double dt = g_data.time_delta;
    // Softening factor to prevent division by zero for close particles
    const double epsilon_sq = 1e-6; 

    for (int t = 0; t < g_data.num_time_steps; ++t) {
        // Step 1: Calculate forces/accelerations for all bodies (O(N^2))
        for (int i = 0; i < n; ++i) {
            double ax = 0.0, ay = 0.0, az = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;

                // Calculate distance vector
                double dx = g_data.pos_x[j] - g_data.pos_x[i];
                double dy = g_data.pos_y[j] - g_data.pos_y[i];
                double dz = g_data.pos_z[j] - g_data.pos_z[i];

                // Calculate squared distance with softening
                double dist_sq = dx * dx + dy * dy + dz * dz + epsilon_sq;

                // Calculate force magnitude based on Newton's law of gravitation
                // (G=1, m_i cancels out for acceleration)
                double inv_dist_cubed = 1.0 / (dist_sq * sqrt(dist_sq));
                double force_mag_over_m_i = g_data.mass[j] * inv_dist_cubed;

                // Accumulate acceleration
                ax += dx * force_mag_over_m_i;
                ay += dy * force_mag_over_m_i;
                az += dz * force_mag_over_m_i;
            }
            g_data.acc_x[i] = ax;
            g_data.acc_y[i] = ay;
            g_data.acc_z[i] = az;
        }

        // Step 2: Update velocities and positions using Euler-Cromer integration
        for (int i = 0; i < n; ++i) {
            // Update velocity
            g_data.vel_x[i] += g_data.acc_x[i] * dt;
            g_data.vel_y[i] += g_data.acc_y[i] * dt;
            g_data.vel_z[i] += g_data.acc_z[i] * dt;
            
            // Update position
            g_data.pos_x[i] += g_data.vel_x[i] * dt;
            g_data.pos_y[i] += g_data.vel_y[i] * dt;
            g_data.pos_z[i] += g_data.vel_z[i] * dt;
        }
    }

    // Calculate total kinetic energy at the end to prevent dead code elimination
    double total_kinetic_energy = 0.0;
    for (int i = 0; i < n; ++i) {
        double vel_sq = g_data.vel_x[i] * g_data.vel_x[i] +
                        g_data.vel_y[i] * g_data.vel_y[i] +
                        g_data.vel_z[i] * g_data.vel_z[i];
        total_kinetic_energy += 0.5 * g_data.mass[i] * vel_sq;
    }
    g_data.final_result = total_kinetic_energy;
}

void cleanup() {
    free(g_data.pos_x);
    free(g_data.pos_y);
    free(g_data.pos_z);
    free(g_data.vel_x);
    free(g_data.vel_y);
    free(g_data.vel_z);
    free(g_data.acc_x);
    free(g_data.acc_y);
    free(g_data.acc_z);
    free(g_data.mass);
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

    // Print final result to stdout
    printf("%.6f\n", g_data.final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
