#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- START of Mersenne Twister (DO NOT MODIFY) ---
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
// --- END of Mersenne Twister ---

// Benchmark parameters and data
int g_num_walkers;
int g_num_time_steps;
int g_num_dimensions;
const double g_step_size = 0.1; // Fixed step size for the random walk

// Data structures
double** g_walkers; // Positions of the walkers [num_walkers][num_dimensions]
double g_total_energy; // Accumulated result

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_walkers> <num_time_steps> <num_dimensions> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_walkers = atoi(argv[1]);
    g_num_time_steps = atoi(argv[2]);
    g_num_dimensions = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if(g_num_walkers <= 0 || g_num_time_steps <= 0 || g_num_dimensions <= 0) {
        fprintf(stderr, "FATAL: Parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    g_walkers = (double**)malloc(g_num_walkers * sizeof(double*));
    if (g_walkers == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for walkers array.\n");
        exit(1);
    }

    for (int i = 0; i < g_num_walkers; i++) {
        g_walkers[i] = (double*)malloc(g_num_dimensions * sizeof(double));
         if (g_walkers[i] == NULL) {
            fprintf(stderr, "FATAL: Memory allocation failed for walker %d.\n", i);
            exit(1);
        }
        // Initialize walkers randomly in a hypercube [-0.5, 0.5]^D
        for (int j = 0; j < g_num_dimensions; j++) {
            g_walkers[i][j] = (double)mt_rand() / (double)UINT32_MAX - 0.5;
        }
    }

    g_total_energy = 0.0;
}

void run_computation() {
    double temp_total_energy = 0.0;
    
    for (int step = 0; step < g_num_time_steps; step++) {
        for (int i = 0; i < g_num_walkers; i++) {
            double r_squared = 0.0;
            
            // Propose a random move for the walker. This is the core of the Monte Carlo method.
            for (int j = 0; j < g_num_dimensions; j++) {
                double displacement = ((double)mt_rand() / (double)UINT32_MAX - 0.5) * g_step_size;
                g_walkers[i][j] += displacement;
                r_squared += g_walkers[i][j] * g_walkers[i][j];
            }

            // Calculate the local energy (using a simple harmonic oscillator potential as a proxy).
            // This step represents the evaluation of physical properties in a simulation.
            // E_local = 0.5 * D + 0.5 * r^2
            double local_energy = 0.5 * g_num_dimensions + 0.5 * r_squared;
            
            temp_total_energy += local_energy;
        }
    }
    g_total_energy = temp_total_energy;
}

void cleanup() {
    if (g_walkers != NULL) {
        for (int i = 0; i < g_num_walkers; i++) {
            if (g_walkers[i] != NULL) {
                free(g_walkers[i]);
            }
        }
        free(g_walkers);
    }
}

int main(int argc, char* argv[]) {
    struct timespec start, end;
    
    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Calculate the final result (average energy per walker-step) to prevent dead code elimination.
    double average_energy = g_total_energy / ((double)g_num_walkers * g_num_time_steps);

    cleanup();

    // Print result to stdout
    printf("%f\n", average_energy);
    
    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
