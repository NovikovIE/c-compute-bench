/*
 * Computational Chemistry: Free Energy Perturbation Benchmark
 *
 * This benchmark simulates a key component of a Free Energy Perturbation (FEP) calculation,
 * often used in drug discovery and materials science to compute the free energy difference
 * between two states (e.g., a ligand binding to a protein).
 *
 * The core of the simulation involves:
 * 1. Defining two end states of a system (State A and State B), characterized by different
 *    Lennard-Jones potential parameters.
 * 2. Creating a series of intermediate states ("lambda windows") that smoothly interpolate
 *    between State A and State B.
 * 3. For each lambda window, running a short molecular dynamics simulation to sample
 *    the system's configurations (atom positions).
 * 4. At each simulation step, calculating the potential energy difference between the next
 *    lambda state and the current one, given the current atom positions.
 * 5. This energy difference is then used in the Zwanzig equation (or related formulas) to
 *    calculate the free energy change. The benchmark accumulates a value based on this
 *    energy difference to produce a final result.
 *
 * The computational workload is dominated by the pairwise potential energy calculation,
 * which scales as O(N^2) with the number of atoms (N). The total complexity is approximately
 * O(num_lambda_windows * time_steps_per_window * num_atoms^2).
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) PRNG --- VERBATIM ---
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
// --- End of MT19937 ---

// --- Benchmark Data and Globals ---
typedef struct {
    int num_atoms;
    int num_lambda_windows;
    int time_steps_per_window;

    float *positions;          // size 3 * num_atoms (x, y, z)
    float *velocities;         // size 3 * num_atoms (vx, vy, vz)
    float *param_A;            // size num_atoms * num_atoms, for State A potential
    float *param_B;            // size num_atoms * num_atoms, for State B potential

    double accumulated_result; // Final result from computation
} BenchmarkData;

BenchmarkData g_data;

// Helper to generate a random float in [min, max]
float random_float(float min, float max) {
    return min + (max - min) * ((float)mt_rand() / (float)UINT32_MAX);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_atoms num_lambda_windows time_steps_per_window seed\n", argv[0]);
        exit(1);
    }

    g_data.num_atoms = atoi(argv[1]);
    g_data.num_lambda_windows = atoi(argv[2]);
    g_data.time_steps_per_window = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    if (g_data.num_atoms <= 1 || g_data.num_lambda_windows <= 1 || g_data.time_steps_per_window <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters. All must be positive, and atoms/lambdas > 1.\n");
        exit(1);
    }

    // Allocate memory
    g_data.positions = (float*)malloc(3 * g_data.num_atoms * sizeof(float));
    g_data.velocities = (float*)malloc(3 * g_data.num_atoms * sizeof(float));
    g_data.param_A = (float*)malloc(g_data.num_atoms * g_data.num_atoms * sizeof(float));
    g_data.param_B = (float*)malloc(g_data.num_atoms * g_data.num_atoms * sizeof(float));
    
    if (!g_data.positions || !g_data.velocities || !g_data.param_A || !g_data.param_B) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize atom positions and velocities
    for (int i = 0; i < g_data.num_atoms * 3; ++i) {
        g_data.positions[i] = random_float(-5.0f, 5.0f);
        g_data.velocities[i] = random_float(-0.1f, 0.1f);
    }

    // Initialize pairwise potential parameters for states A and B
    for (int i = 0; i < g_data.num_atoms; ++i) {
        for (int j = 0; j < g_data.num_atoms; ++j) {
            if (i == j) continue; // No self-interaction
            // Use symmetric parameters for pairs
            size_t idx = i * g_data.num_atoms + j;
            size_t rev_idx = j * g_data.num_atoms + i;
            float val_A = random_float(0.5f, 1.5f);
            float val_B = random_float(0.5f, 1.5f);
            g_data.param_A[idx] = g_data.param_A[rev_idx] = val_A;
            g_data.param_B[idx] = g_data.param_B[rev_idx] = val_B;
        }
    }
    g_data.accumulated_result = 0.0;
}

void run_computation() {
    const float dt = 0.01f; // Timestep for simple dynamics
    const float delta_lambda = 1.0f / (float)(g_data.num_lambda_windows - 1);
    g_data.accumulated_result = 0.0;

    // Loop over each lambda window to run a short simulation
    for (int l = 0; l < g_data.num_lambda_windows; ++l) {
        
        // Run a short MD simulation to sample configurations for this lambda
        for (int t = 0; t < g_data.time_steps_per_window; ++t) {
            
            // Calculate the potential energy difference (d_U) for the perturbation
            double d_U_total = 0.0;
            for (int i = 0; i < g_data.num_atoms; ++i) {
                for (int j = i + 1; j < g_data.num_atoms; ++j) {
                    float dx = g_data.positions[i*3 + 0] - g_data.positions[j*3 + 0];
                    float dy = g_data.positions[i*3 + 1] - g_data.positions[j*3 + 1];
                    float dz = g_data.positions[i*3 + 2] - g_data.positions[j*3 + 2];
                    
                    float r2 = dx*dx + dy*dy + dz*dz;
                    if (r2 < 1e-8f) continue; // Avoid singularity

                    float r6 = r2 * r2 * r2;
                    
                    float C_A = g_data.param_A[i * g_data.num_atoms + j];
                    float C_B = g_data.param_B[i * g_data.num_atoms + j];
                    
                    // dU = d(lambda) * (U_B - U_A)
                    // Simplified potential: U = -C / r^6
                    // U_B - U_A = (-C_B - (-C_A)) / r6 = (C_A - C_B) / r6
                    double d_U_pair = delta_lambda * (C_A - C_B) / r6;
                    d_U_total += d_U_pair;
                }
            }
            
            // Accumulate result using an exponential term (related to Zwanzig formula)
            g_data.accumulated_result += exp(-d_U_total / g_data.num_atoms);
            
            // Update system state using simple Euler integration
            for (int i = 0; i < g_data.num_atoms * 3; ++i) {
                g_data.positions[i] += g_data.velocities[i] * dt;
            }
        }
    }
}

void cleanup() {
    free(g_data.positions);
    free(g_data.velocities);
    free(g_data.param_A);
    free(g_data.param_B);
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
    printf("%f\n", g_data.accumulated_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
