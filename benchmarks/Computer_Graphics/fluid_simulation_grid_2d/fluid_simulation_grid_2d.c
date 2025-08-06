#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK-SPECIFIC GLOBALS ---
#define IX(i, j) ((i) + (grid.width + 2) * (j))
#define SWAP(x0, x) do { double *tmp = x0; x0 = x; x = tmp; } while (0)
#define ITER 16 // Solver iterations

typedef struct {
    int width;
    int height;
    int num_steps;
    double dt;      // time step
    double diff;    // diffusion
    double visc;    // viscosity

    double *u, *v;       // velocity fields (x, y)
    double *u_prev, *v_prev; // previous velocity fields
    double *dens, *dens_prev; // density and previous density
} FluidGrid;

static FluidGrid grid;
static double final_result = 0.0;

// --- BENCHMARK-SPECIFIC LOGIC ---

// Set boundary conditions
void set_bnd(int b, double *x) {
    int N = grid.width;
    int M = grid.height;
    for (int i = 1; i <= M; i++) {
        x[IX(0, i)]     = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
    }
    for (int i = 1; i <= N; i++) {
        x[IX(i, 0)]     = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, M + 1)] = b == 2 ? -x[IX(i, M)] : x[IX(i, M)];
    }

    x[IX(0, 0)]         = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, M + 1)]     = 0.5 * (x[IX(1, M + 1)] + x[IX(0, M)]);
    x[IX(N + 1, 0)]     = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
    x[IX(N + 1, M + 1)] = 0.5 * (x[IX(N, M + 1)] + x[IX(N + 1, M)]);
}

// Gauss-Seidel relaxation solver
void lin_solve(int b, double *x, double *x0, double a, double c) {
    double cRecip = 1.0 / c;
    for (int k = 0; k < ITER; k++) {
        for (int j = 1; j <= grid.height; j++) {
            for (int i = 1; i <= grid.width; i++) {
                x[IX(i, j)] =
                    (x0[IX(i, j)]
                        + a*(x[IX(i + 1, j)]
                            + x[IX(i - 1, j)]
                            + x[IX(i, j + 1)]
                            + x[IX(i, j - 1)]
                        )) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}

// Diffusion step
void diffuse(int b, double *x, double *x0, double diff, double dt) {
    double a = dt * diff * grid.width * grid.height;
    lin_solve(b, x, x0, a, 1 + 4 * a);
}

// Advection step
void advect(int b, double *d, double *d0, double *velocX, double *velocY, double dt) {
    int N = grid.width;
    int M = grid.height;
    double dtx = dt * N;
    double dty = dt * M;

    for (int j = 1; j <= M; j++) {
        for (int i = 1; i <= N; i++) {
            double x = i - dtx * velocX[IX(i, j)];
            double y = j - dty * velocY[IX(i, j)];

            if (x < 0.5) x = 0.5; if (x > N + 0.5) x = N + 0.5;
            int i0 = (int)x; int i1 = i0 + 1;
            if (y < 0.5) y = 0.5; if (y > M + 0.5) y = M + 0.5;
            int j0 = (int)y; int j1 = j0 + 1;

            double s1 = x - i0; double s0 = 1 - s1;
            double t1 = y - j0; double t0 = 1 - t1;

            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                          s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(b, d);
}

// Projection step to enforce mass conservation
void project(double *velocX, double *velocY, double *p, double *div) {
    int N = grid.width;
    int M = grid.height;
    for (int j = 1; j <= M; j++) {
        for (int i = 1; i <= N; i++) {
            div[IX(i, j)] = -0.5*(
                 velocX[IX(i + 1, j)]
                -velocX[IX(i - 1, j)]
                +velocY[IX(i, j + 1)]
                -velocY[IX(i, j - 1)]
                ) / sqrt(N*M);
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 4);

    for (int j = 1; j <= M; j++) {
        for (int i = 1; i <= N; i++) {
            velocX[IX(i, j)] -= 0.5 * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5 * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * M;
        }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

// Single simulation step
void fluid_step() {
    diffuse(1, grid.u_prev, grid.u, grid.visc, grid.dt);
    diffuse(2, grid.v_prev, grid.v, grid.visc, grid.dt);

    project(grid.u_prev, grid.v_prev, grid.u, grid.v);

    advect(1, grid.u, grid.u_prev, grid.u_prev, grid.v_prev, grid.dt);
    advect(2, grid.v, grid.v_prev, grid.u_prev, grid.v_prev, grid.dt);

    project(grid.u, grid.v, grid.u_prev, grid.v_prev);

    diffuse(0, grid.dens_prev, grid.dens, grid.diff, grid.dt);
    advect(0, grid.dens, grid.dens_prev, grid.u, grid.v, grid.dt);
}

void add_source(double *x, double *s, double dt) {
    int size = (grid.width + 2) * (grid.height + 2);
    for(int i = 0; i < size; i++) x[i] += dt * s[i];
}

// --- BENCHMARK FRAMEWORK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_width grid_height num_simulation_steps seed\n", argv[0]);
        exit(1);
    }

    grid.width = atoi(argv[1]);
    grid.height = atoi(argv[2]);
    grid.num_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    grid.dt = 0.1;
    grid.diff = 0.0000; // diffusion rate
    grid.visc = 0.0000; // viscosity

    int size = (grid.width + 2) * (grid.height + 2);

    grid.u = (double*)malloc(size * sizeof(double));
    grid.v = (double*)malloc(size * sizeof(double));
    grid.u_prev = (double*)malloc(size * sizeof(double));
    grid.v_prev = (double*)malloc(size * sizeof(double));
    grid.dens = (double*)malloc(size * sizeof(double));
    grid.dens_prev = (double*)malloc(size * sizeof(double));

    if (!grid.u || !grid.v || !grid.u_prev || !grid.v_prev || !grid.dens || !grid.dens_prev) {
        fprintf(stderr, "Failed to allocate memory for grids.\n");
        exit(1);
    }

    memset(grid.u, 0, size * sizeof(double));
    memset(grid.v, 0, size * sizeof(double));
    memset(grid.u_prev, 0, size * sizeof(double));
    memset(grid.v_prev, 0, size * sizeof(double));
    memset(grid.dens, 0, size * sizeof(double));
    memset(grid.dens_prev, 0, size * sizeof(double));

    // Initial random fluid state
    for (int i = 1; i <= grid.width; i++) {
        for (int j = 1; j <= grid.height; j++) {
            if (mt_rand() % 100 < 5) {
                grid.dens_prev[IX(i, j)] = (mt_rand() / (double)UINT32_MAX) * 500.0;
                grid.u_prev[IX(i, j)] = (mt_rand() / (double)UINT32_MAX - 0.5) * 2.0;
                grid.v_prev[IX(i, j)] = (mt_rand() / (double)UINT32_MAX - 0.5) * 2.0;
            }
        }
    }
}

void run_computation() {
    for (int t = 0; t < grid.num_steps; t++) {
        // Add continuous source in the middle
        int cx = grid.width / 2;
        int cy = grid.height / 2;
        grid.dens_prev[IX(cx, cy)] += 100.0;
        grid.u_prev[IX(cx, cy)] += 0.5;
        grid.v_prev[IX(cx, cy)] += 1.0;

        add_source(grid.u, grid.u_prev, grid.dt);
        add_source(grid.v, grid.v_prev, grid.dt);
        add_source(grid.dens, grid.dens_prev, grid.dt);

        SWAP(grid.u_prev, grid.u);
        SWAP(grid.v_prev, grid.v);
        SWAP(grid.dens_prev, grid.dens);
        
        fluid_step();
         
        // Clear previous states for next iteration's additions
        int size = (grid.width + 2) * (grid.height + 2);
        memset(grid.u_prev, 0, size * sizeof(double));
        memset(grid.v_prev, 0, size * sizeof(double));
        memset(grid.dens_prev, 0, size * sizeof(double));
    }

    double total_density = 0.0;
    int size = (grid.width + 2) * (grid.height + 2);
    for (int i = 0; i < size; i++) {
        total_density += grid.dens[i];
    }
    final_result = total_density;
}

void cleanup() {
    free(grid.u);
    free(grid.v);
    free(grid.u_prev);
    free(grid.v_prev);
    free(grid.dens);
    free(grid.dens_prev);
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
