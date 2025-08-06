#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Mersenne Twister (MT19937) generator - VERBATIM AS PROVIDED
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
// END of Mersenne Twister

// Benchmark parameters
static int GRID_X, GRID_Y, NUM_TIME_STEPS;
static double REYNOLDS_NUMBER;

// D2Q9 model constants
const int Q = 9;
const int ex[] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int ey[] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double w[] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
const int opposite[] = {0, 3, 4, 1, 2, 7, 8, 5, 6};


// Data structures
static double *f, *f_new; // Distribution functions (current and next step)
static int *is_obstacle;  // Boolean flags for obstacle nodes
static double final_result = 0.0;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s grid_dim_x grid_dim_y num_time_steps reynolds_number seed\n", argv[0]);
        exit(1);
    }

    GRID_X = atoi(argv[1]);
    GRID_Y = atoi(argv[2]);
    NUM_TIME_STEPS = atoi(argv[3]);
    REYNOLDS_NUMBER = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    size_t grid_size = (size_t)GRID_X * GRID_Y;
    
    f = (double*)malloc(grid_size * Q * sizeof(double));
    f_new = (double*)malloc(grid_size * Q * sizeof(double));
    is_obstacle = (int*)malloc(grid_size * sizeof(int));
    
    if (!f || !f_new || !is_obstacle) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize obstacles: a cylinder in a channel
    double obstacle_radius = GRID_Y / 9.0;
    double obstacle_cx = GRID_X / 5.0;
    double obstacle_cy = GRID_Y / 2.0;
    for (int y = 0; y < GRID_Y; ++y) {
        for (int x = 0; x < GRID_X; ++x) {
            double dx = x - obstacle_cx;
            double dy = y - obstacle_cy;
            is_obstacle[y * GRID_X + x] = (dx * dx + dy * dy <= obstacle_radius * obstacle_radius);
        }
    }
    
    // Initialize distribution functions to equilibrium (rho=1, u=0)
    for (int y = 0; y < GRID_Y; ++y) {
        for (int x = 0; x < GRID_X; ++x) {
            size_t idx = y * GRID_X + x;
            for (int q = 0; q < Q; ++q) {
                f[q * grid_size + idx] = w[q];
            }
        }
    }
}


void run_computation() {
    size_t grid_size = (size_t)GRID_X * GRID_Y;
    // LBM parameters
    const double u_max = 0.04; // Characteristic velocity, kept low for stability
    const double nu = u_max * (GRID_Y / 9.0 * 2.0) / REYNOLDS_NUMBER; // Viscosity from Re
    const double tau = 3.0 * nu + 0.5; // Relaxation time
    const double omega = 1.0 / tau; // Relaxation frequency
    const double force_x = 1.0e-5; // Body force to drive the flow

    for (int t = 0; t < NUM_TIME_STEPS; ++t) {
        // This loop performs a combined collision and push-streaming step.
        // It reads from `f`, and writes the streamed, post-collision values to `f_new`.
        for (int y = 0; y < GRID_Y; ++y) {
            for (int x = 0; x < GRID_X; ++x) {
                size_t idx = y * GRID_X + x;

                if (is_obstacle[idx]) {
                     // Obstacle nodes perform bounce-back. A particle arriving at an obstacle
                     // is reflected back to where it came from.
                    for (int q = 0; q < Q; ++q) {
                         // This is handled by neighbor fluid cells; the obstacle itself is passive.
                    }
                    continue;
                }

                // --- Collision ---
                double rho = 0.0, ux = 0.0, uy = 0.0;
                
                // 1. Calculate macroscopic properties (density, velocity) from `f`
                for (int q = 0; q < Q; ++q) {
                    double f_val = f[q * grid_size + idx];
                    rho += f_val;
                    ux += f_val * ex[q];
                    uy += f_val * ey[q];
                }
                ux /= rho;
                uy /= rho;
                
                // 2. Apply body force
                ux += force_x * tau / rho;

                // 3. For each direction, calculate equilibrium, collide, and stream
                for (int q = 0; q < Q; ++q) {
                    double u_dot_e = ex[q] * ux + ey[q] * uy;
                    double u_sq = ux * ux + uy * uy;
                    double feq = w[q] * rho * (1.0 + 3.0 * u_dot_e + 4.5 * u_dot_e * u_dot_e - 1.5 * u_sq);
                    
                    double f_current = f[q * grid_size + idx];
                    double f_post_collision = f_current - omega * (f_current - feq);

                    // --- Streaming (Push) ---
                    int dest_x = (x + ex[q] + GRID_X) % GRID_X; // Periodic boundaries
                    int dest_y = (y + ey[q] + GRID_Y) % GRID_Y;
                    size_t dest_idx = dest_y * GRID_X + dest_x;

                    // If the destination is an obstacle, reflect back to the current node.
                    if (is_obstacle[dest_idx]) {
                        f_new[opposite[q] * grid_size + idx] = f_post_collision;
                    } else {
                        f_new[q * grid_size + dest_idx] = f_post_collision;
                    }
                }
            }
        }
        
        // Swap pointers for next iteration
        double *temp = f;
        f = f_new;
        f_new = temp;
    }

    // Calculate final result: average kinetic energy of the fluid
    double total_ke = 0.0;
    long long fluid_cell_count = 0;
    for (int y = 0; y < GRID_Y; ++y) {
        for (int x = 0; x < GRID_X; ++x) {
            size_t idx = y * GRID_X + x;
            if (!is_obstacle[idx]) {
                double rho = 0.0, ux = 0.0, uy = 0.0;
                for (int q = 0; q < Q; ++q) {
                    double f_val = f[q * grid_size + idx];
                    rho += f_val;
                    ux += f_val * ex[q];
                    uy += f_val * ey[q];
                }
                ux /= rho;
                uy /= rho;
                total_ke += 0.5 * rho * (ux * ux + uy * uy);
                fluid_cell_count++;
            }
        }
    }
    final_result = (fluid_cell_count > 0) ? (total_ke / fluid_cell_count) : 0.0;
}

void cleanup() {
    free(f);
    free(f_new);
    free(is_obstacle);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
