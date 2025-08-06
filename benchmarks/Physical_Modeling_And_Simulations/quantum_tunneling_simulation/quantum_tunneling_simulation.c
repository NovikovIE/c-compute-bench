#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Global Data Structures ---
// Simulation Parameters
int grid_points;
int barrier_width;
int num_time_steps;

// Wavefunction arrays (using double buffering)
double *psi_real_current;
double *psi_imag_current;
double *psi_real_next;
double *psi_imag_next;

// Potential barrier
double *potential;

// Final result
double final_result_accumulator;

// --- Mersenne Twister (Verbatim as Required) ---
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

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_points barrier_width num_time_steps seed\n", argv[0]);
        exit(1);
    }

    grid_points = atoi(argv[1]);
    barrier_width = atoi(argv[2]);
    num_time_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if(grid_points <= 0 || barrier_width <= 0 || num_time_steps <= 0 || barrier_width >= grid_points / 2) {
        fprintf(stderr, "Invalid parameters. Ensure grid_points > 0, barrier_width > 0, num_time_steps > 0, and barrier_width < grid_points / 2.\n");
        exit(1);
    }
    
    mt_seed(seed);

    // Allocate memory
    psi_real_current = (double *)malloc(grid_points * sizeof(double));
    psi_imag_current = (double *)malloc(grid_points * sizeof(double));
    psi_real_next = (double *)malloc(grid_points * sizeof(double));
    psi_imag_next = (double *)malloc(grid_points * sizeof(double));
    potential = (double *)malloc(grid_points * sizeof(double));

    if (!psi_real_current || !psi_imag_current || !psi_real_next || !psi_imag_next || !potential) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // --- Initialize Potential Barrier ---
    double barrier_height = 10.0;
    int barrier_start = grid_points / 2 - barrier_width / 2;
    int barrier_end = grid_points / 2 + barrier_width / 2;
    for (int i = 0; i < grid_points; ++i) {
        if (i >= barrier_start && i < barrier_end) {
            potential[i] = barrier_height;
        } else {
            potential[i] = 0.0;
        }
    }

    // --- Initialize Wave Packet (Gaussian) ---
    // Use PRNG to slightly vary the initial momentum
    double k0_base = 25.0;
    double k0_rand_component = (double)mt_rand() / (double)UINT32_MAX;
    double k0 = k0_base + k0_rand_component; // Initial momentum
    
    double dx = 0.1;
    double sigma = 40.0 * dx; // width of wave packet
    double x0 = grid_points * dx / 5.0; // initial position (start of the grid)

    double norm_factor = 0.0;
    for (int i = 0; i < grid_points; ++i) {
        double x = i * dx;
        double exp_arg = -pow(x - x0, 2) / (4.0 * sigma * sigma);
        double cos_arg = k0 * (x - x0);
        double sin_arg = k0 * (x - x0);

        psi_real_current[i] = exp(exp_arg) * cos(cos_arg);
        psi_imag_current[i] = exp(exp_arg) * sin(sin_arg);
        norm_factor += psi_real_current[i] * psi_real_current[i] + psi_imag_current[i] * psi_imag_current[i];
    }
    
    // Normalize the wave function
    norm_factor = sqrt(norm_factor * dx);
    for (int i = 0; i < grid_points; ++i) {
        psi_real_current[i] /= norm_factor;
        psi_imag_current[i] /= norm_factor;
        psi_real_next[i] = 0.0;
        psi_imag_next[i] = 0.0;
    }
}

void run_computation() {
    // Simulation constants
    const double hbar = 1.0;
    const double mass = 2.0;
    const double dt = 0.0001; // Time step
    const double dx = 0.1;    // Spatial step

    const double C1 = -hbar * hbar / (2.0 * mass);
    const double C2 = dt / hbar;
    const double inv_dx2 = 1.0 / (dx*dx);
    
    // Time evolution loop (FDTD)
    for (int t = 0; t < num_time_steps; ++t) {
        // Space loop (update wave function)
        for (int i = 1; i < grid_points - 1; ++i) {
            // Calculate H * psi (Hamiltonian operator on the wave function)
            // Real part
            double laplacian_real = (psi_real_current[i+1] - 2.0 * psi_real_current[i] + psi_real_current[i-1]) * inv_dx2;
            double H_psi_real = C1 * laplacian_real + potential[i] * psi_real_current[i];
            
            // Imaginary part
            double laplacian_imag = (psi_imag_current[i+1] - 2.0 * psi_imag_current[i] + psi_imag_current[i-1]) * inv_dx2;
            double H_psi_imag = C1 * laplacian_imag + potential[i] * psi_imag_current[i];

            // Update using an explicit finite-difference scheme (Euler method)
            psi_real_next[i] = psi_real_current[i] - C2 * H_psi_imag;
            psi_imag_next[i] = psi_imag_current[i] + C2 * H_psi_real;
        }

        // --- Swap buffers for next time step ---
        double *temp_real = psi_real_current;
        psi_real_current = psi_real_next;
        psi_real_next = temp_real;

        double *temp_imag = psi_imag_current;
        psi_imag_current = psi_imag_next;
        psi_imag_next = temp_imag;
    }
    
    // --- Calculate final result to prevent dead code elimination ---
    // Calculate the probability of finding the particle beyond the barrier
    final_result_accumulator = 0.0;
    int barrier_end = grid_points / 2 + barrier_width / 2;
    for (int i = barrier_end; i < grid_points; ++i) {
        final_result_accumulator += (psi_real_current[i] * psi_real_current[i] + psi_imag_current[i] * psi_imag_current[i]) * dx;
    }
}

void cleanup() {
    free(psi_real_current);
    free(psi_imag_current);
    free(psi_real_next);
    free(psi_imag_next);
    free(potential);
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
    printf("%f\n", final_result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
