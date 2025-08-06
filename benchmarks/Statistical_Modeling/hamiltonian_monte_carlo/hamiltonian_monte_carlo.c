#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK PARAMETERS AND DATA ---
typedef struct {
    int num_samples;
    int num_burn_in;
    int num_leapfrog_steps;
    double step_size;
    int num_parameters;

    double **samples;
    double *current_q; // Position vector (parameters)
    
    // Re-usable buffers for computation
    double *p_momentum; // Momentum vector
    double *proposed_q; // Proposed position

    double result_accumulator;
} BenchmarkData;

BenchmarkData g_data;

// --- BENCHMARK FUNCTIONS ---

// Potential Energy U(q) = 0.5 * q^T * q
// The gradient dU/dq_i is simply q_i.
// We use this simple form as our target distribution is a standard multivariate Normal.
double potential_energy(double* q, int n) {
    double U = 0.0;
    for (int i = 0; i < n; i++) {
        U += q[i] * q[i];
    }
    return 0.5 * U;
}

// Kinetic Energy K(p) = 0.5 * p^T * p
double kinetic_energy(double* p, int n) {
    double K = 0.0;
    for (int i = 0; i < n; i++) {
        K += p[i] * p[i];
    }
    return 0.5 * K;
}


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_samples num_burn_in num_leapfrog_steps step_size num_parameters seed\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_burn_in = atoi(argv[2]);
    g_data.num_leapfrog_steps = atoi(argv[3]);
    g_data.step_size = atof(argv[4]);
    g_data.num_parameters = atoi(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);
    mt_seed(seed);

    // Allocate memory for samples and initial state
    g_data.samples = (double**)malloc(g_data.num_samples * sizeof(double*));
    if (!g_data.samples) { perror("malloc failed"); exit(1); }
    for (int i = 0; i < g_data.num_samples; i++) {
        g_data.samples[i] = (double*)malloc(g_data.num_parameters * sizeof(double));
        if (!g_data.samples[i]) { perror("malloc failed"); exit(1); }
    }

    g_data.current_q = (double*)malloc(g_data.num_parameters * sizeof(double));
    if (!g_data.current_q) { perror("malloc failed"); exit(1); }
    
    // Allocate reusable buffers for computation
    g_data.p_momentum = (double*)malloc(g_data.num_parameters * sizeof(double));
    g_data.proposed_q = (double*)malloc(g_data.num_parameters * sizeof(double));
    if (!g_data.p_momentum || !g_data.proposed_q) { perror("malloc failed"); exit(1); }

    // Initialize starting position vector q to zeros
    for (int i = 0; i < g_data.num_parameters; i++) {
        g_data.current_q[i] = 0.0;
    }
}

void run_computation() {
    int sample_idx = 0;
    int total_iterations = g_data.num_samples + g_data.num_burn_in;

    for (int i = 0; i < total_iterations; i++) {
        // 1. Sample momentum p from standard normal N(0,1) using Box-Muller transform
        for (int j = 0; j < g_data.num_parameters; j += 2) {
            double u1 = (mt_rand() + 1.0) / (UINT32_MAX + 2.0); // Ensure u1 is in (0, 1]
            double u2 = mt_rand() / (double)UINT32_MAX;
            double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
            double z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
            g_data.p_momentum[j] = z0;
            if (j + 1 < g_data.num_parameters) {
                g_data.p_momentum[j + 1] = z1;
            }
        }

        // Store current position and compute current energy
        for(int j = 0; j < g_data.num_parameters; j++) {
            g_data.proposed_q[j] = g_data.current_q[j];
        }
        double current_H = potential_energy(g_data.current_q, g_data.num_parameters) + kinetic_energy(g_data.p_momentum, g_data.num_parameters);

        // 2. Leapfrog Integration
        // grad(U) = q. Half step for momentum.
        for (int j = 0; j < g_data.num_parameters; j++) {
            g_data.p_momentum[j] -= g_data.step_size * 0.5 * g_data.proposed_q[j];
        }

        for (int k = 0; k < g_data.num_leapfrog_steps; k++) {
            // Full step for position
            for (int j = 0; j < g_data.num_parameters; j++) {
                g_data.proposed_q[j] += g_data.step_size * g_data.p_momentum[j];
            }
            // Full step for momentum (except for the final step)
            if (k < g_data.num_leapfrog_steps - 1) {
                for (int j = 0; j < g_data.num_parameters; j++) {
                    g_data.p_momentum[j] -= g_data.step_size * g_data.proposed_q[j];
                }
            }
        }

        // Final half step for momentum
        for (int j = 0; j < g_data.num_parameters; j++) {
            g_data.p_momentum[j] -= g_data.step_size * 0.5 * g_data.proposed_q[j];
        }

        // 3. Metropolis-Hastings acceptance step
        double proposed_H = potential_energy(g_data.proposed_q, g_data.num_parameters) + kinetic_energy(g_data.p_momentum, g_data.num_parameters);
        double log_alpha = current_H - proposed_H;

        double u_accept = mt_rand() / (double)UINT32_MAX;
        if (log(u_accept) < log_alpha) {
            // Accepted: update current position
            for (int j = 0; j < g_data.num_parameters; j++) {
                g_data.current_q[j] = g_data.proposed_q[j];
            }
        }

        // 4. Store sample if we are past the burn-in phase
        if (i >= g_data.num_burn_in) {
            for (int j = 0; j < g_data.num_parameters; j++) {
                g_data.samples[sample_idx][j] = g_data.current_q[j];
            }
            sample_idx++;
        }
    }

    // Accumulate a final result to prevent dead code elimination
    g_data.result_accumulator = 0.0;
    for (int i = 0; i < g_data.num_samples; i++) {
        for (int j = 0; j < g_data.num_parameters; j++) {
            g_data.result_accumulator += g_data.samples[i][j];
        }
    }
}

void cleanup() {
    for (int i = 0; i < g_data.num_samples; i++) {
        free(g_data.samples[i]);
    }
    free(g_data.samples);
    free(g_data.current_q);
    free(g_data.p_momentum);
    free(g_data.proposed_q);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%f\n", g_data.result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
