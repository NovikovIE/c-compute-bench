#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// START MT19937
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
// END MT19937

// --- BENCHMARK DATA AND PARAMETERS ---

// Parameters
int num_links;
int num_simulation_steps;

// Robot State Data
double* q;      // Joint positions
double* q_dot;  // Joint velocities
double* q_ddot; // Joint accelerations
double* tau;    // Joint torques (input)

// Intermediate computation buffers to simulate dynamics calculations
double* pass1_data; // Simulates kinematic propagation (e.g., spatial velocities)
double* pass2_data; // Simulates inertia/force propagation

// Final result accumulator
double final_result = 0.0;

// Helper to generate a double between -1.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX * 2.0 - 1.0;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_links> <num_simulation_steps> <seed>\n", argv[0]);
        exit(1);
    }

    num_links = atoi(argv[1]);
    num_simulation_steps = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    q      = (double*)malloc(num_links * sizeof(double));
    q_dot  = (double*)malloc(num_links * sizeof(double));
    q_ddot = (double*)malloc(num_links * sizeof(double));
    tau    = (double*)malloc(num_links * sizeof(double));

    // Using a factor of 6 to simulate 6D spatial vectors
    pass1_data = (double*)malloc(num_links * 6 * sizeof(double));
    pass2_data = (double*)malloc(num_links * 6 * sizeof(double));

    if (!q || !q_dot || !q_ddot || !tau || !pass1_data || !pass2_data) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_links; ++i) {
        q[i] = rand_double() * M_PI; // Initial angle
        q_dot[i] = rand_double() * 0.5; // Initial velocity
        q_ddot[i] = 0.0;
        tau[i] = rand_double() * 10.0; // Input torque
    }
}

void run_computation() {
    const double dt = 0.01; // Simulation time step

    for (int step = 0; step < num_simulation_steps; ++step) {

        // --- Pass 1: Outward recursion (base to tip) ---
        // This pass simulates the propagation of kinematic information, like velocities,
        // from the base of the robot arm to its tip.
        for (int k = 0; k < 6; k++) {
            pass1_data[k] = sin(q[0] + (double)k) * q_dot[0] * 0.1;
        }
        for (int i = 1; i < num_links; ++i) {
            for (int k = 0; k < 6; ++k) {
                double prev_val = pass1_data[(i - 1) * 6 + k];
                pass1_data[i * 6 + k] = prev_val * cos(q[i]) + sin(q_dot[i] + (double)k) * 0.5;
            }
        }

        // --- Pass 2: Inward recursion (tip to base) ---
        // This pass simulates the calculation of articulated-body inertias and bias forces,
        // propagating from the tip back to the base.
        for (int k = 0; k < 6; ++k) {
            pass2_data[(num_links - 1) * 6 + k] = cos(q[num_links - 1] + tau[num_links-1]) * 0.2;
        }
        for (int i = num_links - 2; i >= 0; --i) {
            for (int k = 0; k < 6; ++k) {
                double child_val = pass2_data[(i + 1) * 6 + k];
                double p1_val = pass1_data[i * 6 + k];
                pass2_data[i * 6 + k] = child_val * fmod(tau[i] + 1.1, 2.0) + p1_val * sin(q[i]);
            }
        }

        // --- Pass 3: Outward recursion (base to tip) ---
        // This pass simulates solving for joint accelerations, starting from the base
        // (whose acceleration is now known from Pass 2) and propagating to the tip.
        q_ddot[0] = pass2_data[0] * 0.05 + tau[0];
        for (int k = 1; k < 6; k++) q_ddot[0] += pass2_data[k]; // Fold in other spatial components
        
        for (int i = 1; i < num_links; ++i) {
            double p1_influence = 0.0;
            double p2_influence = 0.0;
            for(int k = 0; k < 6; ++k) {
                p1_influence += pass1_data[i * 6 + k];
                p2_influence += pass2_data[i * 6 + k];
            }
            q_ddot[i] = q_ddot[i-1] * 0.7 + p1_influence * 0.1 + p2_influence * 0.2 + tau[i];
        }
        
        // --- Integration Step ---
        // Update joint positions and velocities using the calculated accelerations.
        for (int i = 0; i < num_links; ++i) {
            q_dot[i] += q_ddot[i] * dt;
            q[i] += q_dot[i] * dt;
            q[i] = fmod(q[i], 2.0 * M_PI);
        }
    }

    // Accumulate a result to prevent dead-code elimination
    double sum = 0.0;
    for (int i = 0; i < num_links; ++i) {
        sum += q[i] + q_dot[i];
    }
    final_result = sum;
}

void cleanup() {
    free(q);
    free(q_dot);
    free(q_ddot);
    free(tau);
    free(pass1_data);
    free(pass2_data);
    q = q_dot = q_ddot = tau = pass1_data = pass2_data = NULL;
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

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
