#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

//
// Mersenne Twister Generator (from problem description)
//
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

//
// Benchmark: 1D Schrodinger Equation Solver
//

// Global variables for benchmark data
int grid_points;
int num_time_steps;
double *psi_real;      // Real part of the wave function
double *psi_imag;      // Imaginary part of the wave function
double *psi_real_next;  // Next time step's real part
double *psi_imag_next;  // Next time step's imaginary part
double *potential;     // Potential energy field V(x)
double final_result;   // Accumulated result to prevent dead code elimination

// Physics and simulation constants
const double DX = 0.02; // Spatial step
const double DT = 1e-7; // Time step

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s grid_points num_time_steps seed\n", argv[0]);
        exit(1);
    }

    grid_points = atoi(argv[1]);
    num_time_steps = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (grid_points <= 1 || num_time_steps <= 0) {
        fprintf(stderr, "ERROR: grid_points must be > 1 and num_time_steps must be > 0.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    psi_real = (double *)malloc(grid_points * sizeof(double));
    psi_imag = (double *)malloc(grid_points * sizeof(double));
    psi_real_next = (double *)malloc(grid_points * sizeof(double));
    psi_imag_next = (double *)malloc(grid_points * sizeof(double));
    potential = (double *)malloc(grid_points * sizeof(double));

    if (!psi_real || !psi_imag || !psi_real_next || !psi_imag_next || !potential) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize potential (a double barrier)
    double barrier_pos_1 = grid_points * DX * 0.45;
    double barrier_pos_2 = grid_points * DX * 0.55;
    double barrier_width = grid_points * DX * 0.02;
    double barrier_height = (mt_rand() / (double)UINT32_MAX) * 1e5 + 1e5;
    for (int i = 0; i < grid_points; ++i) {
        double x = i * DX;
        if ((fabs(x - barrier_pos_1) < barrier_width / 2.0) || (fabs(x - barrier_pos_2) < barrier_width / 2.0)) {
            potential[i] = barrier_height;
        } else {
            potential[i] = 0.0;
        }
    }

    // Initialize wave function (a Gaussian wave packet)
    double x0 = grid_points * DX * 0.25;    // Initial position
    double sigma = grid_points * DX * 0.05;   // Packet width
    double p0 = 75.0;                         // Momentum
    double norm = 0.0;

    for (int i = 0; i < grid_points; ++i) {
        double x = i * DX;
        double exp_arg = -pow(x - x0, 2) / (2.0 * pow(sigma, 2));
        double real_part = exp(exp_arg) * cos(p0 * (x - x0));
        double imag_part = exp(exp_arg) * sin(p0 * (x - x0));
        psi_real[i] = real_part;
        psi_imag[i] = imag_part;
        norm += (real_part * real_part + imag_part * imag_part);
    }

    // Normalize the wave function
    norm = sqrt(norm * DX);
    for (int i = 0; i < grid_points; ++i) {
        psi_real[i] /= norm;
        psi_imag[i] /= norm;
    }

    // Ensure boundaries are zero
    psi_real[0] = psi_imag[0] = 0.0;
    psi_real[grid_points - 1] = psi_imag[grid_points - 1] = 0.0;
}

void run_computation() {
    double dt_dx2 = DT / (DX * DX);

    for (int t = 0; t < num_time_steps; ++t) {
        #pragma omp parallel for
        for (int i = 1; i < grid_points - 1; ++i) {
            // Compute Laplacian of real and imaginary parts
            double laplacian_real = (psi_real[i + 1] - 2 * psi_real[i] + psi_real[i - 1]);
            double laplacian_imag = (psi_imag[i + 1] - 2 * psi_imag[i] + psi_imag[i - 1]);

            // Update using Finite-Difference Time-Domain (FDTD)
            // i d(psi)/dt = -hbar^2/2m * d^2(psi)/dx^2 + V*psi
            // For hbar=1, 2m=1: d(psi_r)/dt = - H(psi_i), d(psi_i)/dt = H(psi_r)
            // H(psi) = -d^2(psi)/dx^2 + V*psi
            psi_real_next[i] = psi_real[i] - dt_dx2 * laplacian_imag + DT * potential[i] * psi_imag[i];
            psi_imag_next[i] = psi_imag[i] + dt_dx2 * laplacian_real - DT * potential[i] * psi_real[i];
        }

        // Swap buffers using memcpy (could also swap pointers)
        memcpy(psi_real, psi_real_next, grid_points * sizeof(double));
        memcpy(psi_imag, psi_imag_next, grid_points * sizeof(double));
    }

    // Calculate the total probability as the final result
    double total_probability = 0.0;
    for (int i = 0; i < grid_points; ++i) {
        total_probability += (psi_real[i] * psi_real[i] + psi_imag[i] * psi_imag[i]) * DX;
    }
    final_result = total_probability;
}

void cleanup() {
    free(psi_real);
    free(psi_imag);
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

    // Print result to stdout (total probability, should be ~1.0)
    printf("%.10f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
