#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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

// --- BENCHMARK DATA AND PARAMETERS ---
typedef struct {
    int grid_width, grid_height, grid_depth;
    int W, H, D; // Padded dimensions
    int size; // Total grid cells
    int num_simulation_steps;

    float dt; // time step
    float diff; // diffusion
    float visc; // viscosity

    float *vx, *vy, *vz;       // velocity
    float *vx0, *vy0, *vz0;    // previous velocity
    float *density, *density0; // density and previous density
    float *p, *div;            // pressure and divergence

    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

#define IX(i, j, k) ((i) + g_data.W * ((j) + g_data.H * (k)))
#define SWAP(x0, x) do { float *tmp = x0; x0 = x; x = tmp; } while (0)

// --- PROTOTYPES for internal functions ---
static void fluid_step(void);
static void set_bnd(int b, float *x);

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s grid_width grid_height grid_depth num_simulation_steps seed\n", argv[0]);
        exit(1);
    }

    g_data.grid_width = atoi(argv[1]);
    g_data.grid_height = atoi(argv[2]);
    g_data.grid_depth = atoi(argv[3]);
    g_data.num_simulation_steps = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    g_data.W = g_data.grid_width + 2;
    g_data.H = g_data.grid_height + 2;
    g_data.D = g_data.grid_depth + 2;
    g_data.size = g_data.W * g_data.H * g_data.D;

    g_data.dt = 0.1f;
    g_data.diff = 0.0f; // No diffusion for density
    g_data.visc = 0.00001f; // low viscosity

    g_data.vx = (float *)calloc(g_data.size, sizeof(float));
    g_data.vy = (float *)calloc(g_data.size, sizeof(float));
    g_data.vz = (float *)calloc(g_data.size, sizeof(float));
    g_data.vx0 = (float *)calloc(g_data.size, sizeof(float));
    g_data.vy0 = (float *)calloc(g_data.size, sizeof(float));
    g_data.vz0 = (float *)calloc(g_data.size, sizeof(float));
    g_data.density = (float *)calloc(g_data.size, sizeof(float));
    g_data.density0 = (float *)calloc(g_data.size, sizeof(float));
    g_data.p = (float *)calloc(g_data.size, sizeof(float));
    g_data.div = (float *)calloc(g_data.size, sizeof(float));

    if (!g_data.vx || !g_data.vy || !g_data.vz || !g_data.vx0 || !g_data.vy0 || !g_data.vz0 || !g_data.density || !g_data.density0 || !g_data.p || !g_data.div) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Seed the random number generator
    mt_seed(seed);

    // Initial state: a block of density in the center with an upward velocity
    int cx_start = g_data.grid_width / 2 - 2;
    int cx_end = g_data.grid_width / 2 + 2;
    int cy_start = g_data.grid_height / 4;
    int cy_end = g_data.grid_height / 4 + 4;
    int cz_start = g_data.grid_depth / 2 - 2;
    int cz_end = g_data.grid_depth / 2 + 2;

    for (int k = cz_start; k <= cz_end; ++k) {
        for (int j = cy_start; j <= cy_end; ++j) {
            for (int i = cx_start; i <= cx_end; ++i) {
                g_data.density[IX(i, j, k)] = 100.0f;
                g_data.vy[IX(i, j, k)] = 2.0f; 
            }
        }
    }

    // Add small random fluctuations to velocity to satisfy MT19937 usage requirement
    for(int i = 0; i < g_data.size; ++i) {
        g_data.vx[i] += (mt_rand() / (float)UINT32_MAX - 0.5f) * 0.01f;
        g_data.vy[i] += (mt_rand() / (float)UINT32_MAX - 0.5f) * 0.01f;
        g_data.vz[i] += (mt_rand() / (float)UINT32_MAX - 0.5f) * 0.01f;
    }
}

void run_computation() {
    for (int t = 0; t < g_data.num_simulation_steps; ++t) {
        fluid_step();
    }

    double total_density = 0.0;
    for (int i = 0; i < g_data.size; ++i) {
        total_density += g_data.density[i];
    }
    g_data.final_result = total_density;
}

void cleanup() {
    free(g_data.vx);
    free(g_data.vy);
    free(g_data.vz);
    free(g_data.vx0);
    free(g_data.vy0);
    free(g_data.vz0);
    free(g_data.density);
    free(g_data.density0);
    free(g_data.p);
    free(g_data.div);
}

// --- INTERNAL ALGORITHMS ---

static void lin_solve(int b, float *x, float *x0, float a, float c, int iter) {
    float c_recip = 1.0f / c;
    for (int iter_k = 0; iter_k < iter; ++iter_k) {
        for (int k = 1; k <= g_data.grid_depth; ++k) {
            for (int j = 1; j <= g_data.grid_height; ++j) {
                for (int i = 1; i <= g_data.grid_width; ++i) {
                    x[IX(i, j, k)] = 
                        (x0[IX(i, j, k)] + a * (x[IX(i-1, j, k)] + x[IX(i+1, j, k)] +
                                                 x[IX(i, j-1, k)] + x[IX(i, j+1, k)] +
                                                 x[IX(i, j, k-1)] + x[IX(i, j, k+1)])) * c_recip;
                }
            }
        }
        set_bnd(b, x);
    }
}

static void diffuse(int b, float *x, float *x0, float diff, float dt, int iter) {
    float a = dt * diff * g_data.grid_width * g_data.grid_height * g_data.grid_depth;
    lin_solve(b, x, x0, a, 1 + 6 * a, iter);
}

static void advect(int b, float *d, float *d0, float *velocX, float *velocY, float *velocZ, float dt) {
    float dtx = dt * g_data.grid_width;
    float dty = dt * g_data.grid_height;
    float dtz = dt * g_data.grid_depth;

    for (int k = 1; k <= g_data.grid_depth; ++k) {
        for (int j = 1; j <= g_data.grid_height; ++j) {
            for (int i = 1; i <= g_data.grid_width; ++i) {
                float x = (float)i - dtx * velocX[IX(i,j,k)];
                float y = (float)j - dty * velocY[IX(i,j,k)];
                float z = (float)k - dtz * velocZ[IX(i,j,k)];

                if (x < 0.5f) x = 0.5f; if (x > g_data.grid_width + 0.5f) x = g_data.grid_width + 0.5f;
                int i0 = (int)x; int i1 = i0 + 1;
                if (y < 0.5f) y = 0.5f; if (y > g_data.grid_height + 0.5f) y = g_data.grid_height + 0.5f;
                int j0 = (int)y; int j1 = j0 + 1;
                if (z < 0.5f) z = 0.5f; if (z > g_data.grid_depth + 0.5f) z = g_data.grid_depth + 0.5f;
                int k0 = (int)z; int k1 = k0 + 1;

                float s1 = x - i0; float s0 = 1 - s1;
                float t1 = y - j0; float t0 = 1 - t1;
                float u1 = z - k0; float u0 = 1 - u1;

                d[IX(i,j,k)] =
                    s0 * (t0 * (u0*d0[IX(i0,j0,k0)] + u1*d0[IX(i0,j0,k1)]) +
                          t1 * (u0*d0[IX(i0,j1,k0)] + u1*d0[IX(i0,j1,k1)])) +
                    s1 * (t0 * (u0*d0[IX(i1,j0,k0)] + u1*d0[IX(i1,j0,k1)]) +
                          t1 * (u0*d0[IX(i1,j1,k0)] + u1*d0[IX(i1,j1,k1)]));
            }
        }
    }
    set_bnd(b, d);
}

static void project(float *velocX, float *velocY, float *velocZ, float *p, float *div, int iter) {
    for (int k = 1; k <= g_data.grid_depth; ++k) {
        for (int j = 1; j <= g_data.grid_height; ++j) {
            for (int i = 1; i <= g_data.grid_width; ++i) {
                div[IX(i,j,k)] = -0.5f * (velocX[IX(i+1,j,k)] - velocX[IX(i-1,j,k)] +
                                         velocY[IX(i,j+1,k)] - velocY[IX(i,j-1,k)] +
                                         velocZ[IX(i,j,k+1)] - velocZ[IX(i,j,k-1)]);
                p[IX(i,j,k)] = 0;
            }
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6, iter);

    for (int k = 1; k <= g_data.grid_depth; ++k) {
        for (int j = 1; j <= g_data.grid_height; ++j) {
            for (int i = 1; i <= g_data.grid_width; ++i) {
                velocX[IX(i,j,k)] -= 0.5f * (p[IX(i+1,j,k)] - p[IX(i-1,j,k)]);
                velocY[IX(i,j,k)] -= 0.5f * (p[IX(i,j+1,k)] - p[IX(i,j-1,k)]);
                velocZ[IX(i,j,k)] -= 0.5f * (p[IX(i,j,k+1)] - p[IX(i,j,k-1)]);
            }
        }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
    set_bnd(3, velocZ);
}

static void fluid_step(void) {
    int N = g_data.grid_width;
    float visc = g_data.visc;
    float dt = g_data.dt;
    float *vx = g_data.vx, *vy = g_data.vy, *vz = g_data.vz;
    float *vx0 = g_data.vx0, *vy0 = g_data.vy0, *vz0 = g_data.vz0;
    float *dens = g_data.density, *dens0 = g_data.density0;
    int iter = 20; // Number of solver iterations

    // Velocity diffusion/viscosity
    diffuse(1, vx0, vx, visc, dt, iter);
    diffuse(2, vy0, vy, visc, dt, iter);
    diffuse(3, vz0, vz, visc, dt, iter);
    project(vx0, vy0, vz0, g_data.p, g_data.div, iter);
    
    // Advect velocity
    advect(1, vx, vx0, vx0, vy0, vz0, dt);
    advect(2, vy, vy0, vx0, vy0, vz0, dt);
    advect(3, vz, vz0, vx0, vy0, vz0, dt);
    project(vx, vy, vz, g_data.p, g_data.div, iter);
    
    // Advect density
    SWAP(dens0, dens);
    diffuse(0, dens, dens0, g_data.diff, dt, iter);
    SWAP(dens0, dens);
    advect(0, dens, dens0, vx, vy, vz, dt);
}


static void set_bnd(int b, float *x) {
    // Walls are boundaries
    for (int k = 1; k <= g_data.grid_depth; ++k) {
        for (int j = 1; j <= g_data.grid_height; ++j) {
            x[IX(0, j, k)] = b == 1 ? -x[IX(1, j, k)] : x[IX(1, j, k)];
            x[IX(g_data.grid_width + 1, j, k)] = b == 1 ? -x[IX(g_data.grid_width, j, k)] : x[IX(g_data.grid_width, j, k)];
        }
    }
    for (int k = 1; k <= g_data.grid_depth; ++k) {
        for (int i = 1; i <= g_data.grid_width; ++i) {
            x[IX(i, 0, k)] = b == 2 ? -x[IX(i, 1, k)] : x[IX(i, 1, k)];
            x[IX(i, g_data.grid_height + 1, k)] = b == 2 ? -x[IX(i, g_data.grid_height, k)] : x[IX(i, g_data.grid_height, k)];
        }
    }
    for (int j = 1; j <= g_data.grid_height; ++j) {
        for (int i = 1; i <= g_data.grid_width; ++i) {
            x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
            x[IX(i, j, g_data.grid_depth + 1)] = b == 3 ? -x[IX(i, j, g_data.grid_depth)] : x[IX(i, j, g_data.grid_depth)];
        }
    }
    // Edges and corners
    x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
    x[IX(0, g_data.H-1, 0)] = 0.33f * (x[IX(1, g_data.H-1, 0)] + x[IX(0, g_data.H-2, 0)] + x[IX(0, g_data.H-1, 1)]);
    x[IX(g_data.W-1, 0, 0)] = 0.33f * (x[IX(g_data.W-2, 0, 0)] + x[IX(g_data.W-1, 1, 0)] + x[IX(g_data.W-1, 0, 1)]);
    x[IX(0, 0, g_data.D-1)] = 0.33f * (x[IX(1, 0, g_data.D-1)] + x[IX(0, 1, g_data.D-1)] + x[IX(0, 0, g_data.D-2)]);
    x[IX(g_data.W-1, g_data.H-1, 0)] = 0.33f * (x[IX(g_data.W-2, g_data.H-1, 0)] + x[IX(g_data.W-1, g_data.H-2, 0)] + x[IX(g_data.W-1, g_data.H-1, 1)]);
    x[IX(0, g_data.H-1, g_data.D-1)] = 0.33f * (x[IX(1, g_data.H-1, g_data.D-1)] + x[IX(0, g_data.H-2, g_data.D-1)] + x[IX(0, g_data.H-1, g_data.D-2)]);
    x[IX(g_data.W-1, 0, g_data.D-1)] = 0.33f * (x[IX(g_data.W-2, 0, g_data.D-1)] + x[IX(g_data.W-1, 1, g_data.D-1)] + x[IX(g_data.W-1, 0, g_data.D-2)]);
    x[IX(g_data.W-1, g_data.H-1, g_data.D-1)] = 0.33f * (x[IX(g_data.W-2, g_data.H-1, g_data.D-1)] + x[IX(g_data.W-1, g_data.H-2, g_data.D-1)] + x[IX(g_data.W-1, g_data.H-1, g_data.D-2)]);
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
    printf("%f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
