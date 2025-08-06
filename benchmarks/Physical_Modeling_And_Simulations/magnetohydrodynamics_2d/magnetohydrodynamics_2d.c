#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of MT19937 ---

// --- Benchmark Globals ---
typedef struct {
    int grid_dim_x;
    int grid_dim_y;
    int num_time_steps;
    long grid_size;

    // Conserved variable arrays for the 2D grid
    // rho: density, mx/my: momentum, Bx/By: magnetic field, E: total energy
    double *rho, *mx, *my, *E, *Bx, *By; // Current state
    double *rho_new, *mx_new, *my_new, *E_new, *Bx_new, *By_new; // Next state

    double final_result;
} BenchmarkData;

static BenchmarkData B;

// Physical constants
const double GAMMA_CONST = 5.0 / 3.0; // Adiabatic index
const double DT = 0.01;      // Time step

// --- Function Prototypes ---
void setup_benchmark(int argc, char* argv[]);
void run_computation();
void cleanup();

// --- Main --- 
int main(int argc, char* argv[]) {
    struct timespec start, end;

    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_dim_x grid_dim_y num_time_steps seed\n", argv[0]);
        return 1;
    }

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", B.final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}

// --- Data Setup ---
void setup_benchmark(int argc, char* argv[]) {
    B.grid_dim_x = atoi(argv[1]);
    B.grid_dim_y = atoi(argv[2]);
    B.num_time_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    B.grid_size = (long)B.grid_dim_x * B.grid_dim_y;

    // Allocate all data arrays
    B.rho = (double*)malloc(B.grid_size * sizeof(double));
    B.mx = (double*)malloc(B.grid_size * sizeof(double));
    B.my = (double*)malloc(B.grid_size * sizeof(double));
    B.E = (double*)malloc(B.grid_size * sizeof(double));
    B.Bx = (double*)malloc(B.grid_size * sizeof(double));
    B.By = (double*)malloc(B.grid_size * sizeof(double));
    B.rho_new = (double*)malloc(B.grid_size * sizeof(double));
    B.mx_new = (double*)malloc(B.grid_size * sizeof(double));
    B.my_new = (double*)malloc(B.grid_size * sizeof(double));
    B.E_new = (double*)malloc(B.grid_size * sizeof(double));
    B.Bx_new = (double*)malloc(B.grid_size * sizeof(double));
    B.By_new = (double*)malloc(B.grid_size * sizeof(double));

    // Check for allocation failures
    if (!B.rho || !B.mx || !B.my || !B.E || !B.Bx || !B.By || !B.rho_new || !B.mx_new || !B.my_new || !B.E_new || !B.Bx_new || !B.By_new) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize grid with random data
    for (long i = 0; i < B.grid_size; ++i) {
        // Start with primitive variables: density, velocity, pressure, magnetic field
        double r  = (mt_rand() / (double)UINT32_MAX) * 0.5 + 0.5; // rho in [0.5, 1.0]
        double vx = (mt_rand() / (double)UINT32_MAX) * 0.2 - 0.1; // vx in [-0.1, 0.1]
        double vy = (mt_rand() / (double)UINT32_MAX) * 0.2 - 0.1; // vy in [-0.1, 0.1]
        double p  = (mt_rand() / (double)UINT32_MAX) * 0.8 + 0.2; // pressure in [0.2, 1.0]
        double b_x = (mt_rand() / (double)UINT32_MAX) * 0.2 - 0.1; // Bx in [-0.1, 0.1]
        double b_y = (mt_rand() / (double)UINT32_MAX) * 0.2 - 0.1; // By in [-0.1, 0.1]

        // Compute and store conserved variables
        B.rho[i] = r;
        B.mx[i] = r * vx;
        B.my[i] = r * vy;
        B.Bx[i] = b_x;
        B.By[i] = b_y;
        B.E[i] = p / (GAMMA_CONST - 1.0) + 0.5 * r * (vx*vx + vy*vy) + 0.5 * (b_x*b_x + b_y*b_y);
    }
}

// --- Computation ---
void run_computation() {
    for (int t = 0; t < B.num_time_steps; ++t) {
        // Update only interior grid points
        for (int j = 1; j < B.grid_dim_y - 1; ++j) {
            for (int i = 1; i < B.grid_dim_x - 1; ++i) {
                long idx = j * B.grid_dim_x + i;
                long ix_p1 = idx + 1;
                long ix_m1 = idx - 1;
                long iy_p1 = idx + B.grid_dim_x;
                long iy_m1 = idx - B.grid_dim_x;

                // Simplified stencil computation for MHD-like equations
                // Not physically accurate, but represents the computational load and data access patterns.
                double inv_rho = 1.0 / (B.rho[idx] + 1e-9); // Add epsilon to avoid division by zero
                double vx = B.mx[idx] * inv_rho;
                double vy = B.my[idx] * inv_rho;

                // 1. Update Density (rho) - Advection-diffusion
                double rho_advect = -vx * (B.rho[ix_p1] - B.rho[ix_m1]) - vy * (B.rho[iy_p1] - B.rho[iy_m1]);
                double rho_diffuse = 0.05 * (B.rho[ix_p1] + B.rho[ix_m1] + B.rho[iy_p1] + B.rho[iy_m1] - 4.0 * B.rho[idx]);
                B.rho_new[idx] = B.rho[idx] + DT * (rho_advect + rho_diffuse);
                if (B.rho_new[idx] < 1e-3) B.rho_new[idx] = 1e-3; // Clamp density to prevent instability

                // 2. Update Momentum (mx, my) - Advection, diffusion, and simplified Lorentz force
                double mx_advect = -vx * (B.mx[ix_p1] - B.mx[ix_m1]) - vy * (B.mx[iy_p1] - B.mx[iy_m1]);
                double my_advect = -vx * (B.my[ix_p1] - B.my[ix_m1]) - vy * (B.my[iy_p1] - B.my[iy_m1]);
                double lorentz_x = B.By[idx] * (B.By[iy_p1] - B.By[iy_m1]); // Simplified force terms
                double lorentz_y = B.Bx[idx] * (B.Bx[ix_p1] - B.Bx[ix_m1]);
                B.mx_new[idx] = B.mx[idx] + DT * (mx_advect + lorentz_x);
                B.my_new[idx] = B.my[idx] + DT * (my_advect - lorentz_y);

                // 3. Update Magnetic Field (Bx, By) - Advection and simplified induction
                double Bx_advect = -vx * (B.Bx[ix_p1] - B.Bx[ix_m1]) - vy * (B.Bx[iy_p1] - B.Bx[iy_m1]);
                double By_advect = -vx * (B.By[ix_p1] - B.By[ix_m1]) - vy * (B.By[iy_p1] - B.By[iy_m1]);
                B.Bx_new[idx] = B.Bx[idx] + DT * Bx_advect;
                B.By_new[idx] = B.By[idx] + DT * By_advect;

                // 4. Update Energy (E) - Advection and work term
                double E_advect = -vx * (B.E[ix_p1] - B.E[ix_m1]) - vy * (B.E[iy_p1] - B.E[iy_m1]);
                double work = 0.01 * (B.mx[idx] * lorentz_x - B.my[idx] * lorentz_y);
                B.E_new[idx] = B.E[idx] + DT * (E_advect + work);
            }
        }

        // Swap pointers for the next time step (more efficient than copying)
        double *temp_ptr;
        temp_ptr = B.rho; B.rho = B.rho_new; B.rho_new = temp_ptr;
        temp_ptr = B.mx; B.mx = B.mx_new; B.mx_new = temp_ptr;
        temp_ptr = B.my; B.my = B.my_new; B.my_new = temp_ptr;
        temp_ptr = B.E; B.E = B.E_new; B.E_new = temp_ptr;
        temp_ptr = B.Bx; B.Bx = B.Bx_new; B.Bx_new = temp_ptr;
        temp_ptr = B.By; B.By = B.By_new; B.By_new = temp_ptr;
    }

    // Calculate final checksum to prevent dead code elimination and provide a result
    double checksum = 0.0;
    for (long i = 0; i < B.grid_size; ++i) {
        checksum += B.rho[i] + B.E[i];
    }
    B.final_result = checksum;
}

// --- Cleanup ---
void cleanup() {
    free(B.rho); free(B.mx);  free(B.my); free(B.E); free(B.Bx); free(B.By);
    free(B.rho_new); free(B.mx_new); free(B.my_new); free(B.E_new); free(B.Bx_new); free(B.By_new); 
}
