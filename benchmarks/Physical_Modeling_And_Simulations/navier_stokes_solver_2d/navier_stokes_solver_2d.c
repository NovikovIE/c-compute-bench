#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

// --- BENCHMARK DATA AND PARAMETERS ---

// Parameters
int grid_dim_x, grid_dim_y;
int N_x, N_y; // Padded grid dimensions (N+2)
int num_time_steps;
float viscosity;

// Data arrays
float *u, *v;       // velocity fields
float *u_prev, *v_prev; // previous velocity fields
float *p, *div_v;   // pressure and divergence fields (reused as scratch space)

// Result accumulator
double final_result;

// --- HELPER FUNCTIONS FOR COMPUTATION ---

#define IX(i, j) ((i) + (N_x) * (j))
#define SWAP(x0, x) {float *tmp=x0; x0=x; x=tmp;}

// Forward declarations
void set_bnd(int b, float *x);
void lin_solve(int b, float *x, float *x0, float a, float c);
void diffuse(int b, float *x, float *x0, float diff, float dt);
void advect(int b, float *d, float *d0, float *u, float *v, float dt);
void project(float *u, float *v, float *p, float *div);

// --- BENCHMARK CORE FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s grid_dim_x grid_dim_y num_time_steps viscosity seed\n", argv[0]);
        exit(1);
    }

    grid_dim_x = atoi(argv[1]);
    grid_dim_y = atoi(argv[2]);
    num_time_steps = atoi(argv[3]);
    viscosity = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    N_x = grid_dim_x + 2;
    N_y = grid_dim_y + 2;
    size_t size = (size_t)N_x * N_y * sizeof(float);

    u       = (float*)malloc(size);
    v       = (float*)malloc(size);
    u_prev  = (float*)malloc(size);
    v_prev  = (float*)malloc(size);
    p       = (float*)malloc(size);
    div_v   = (float*)malloc(size);

    // Zero out all fields initially
    memset(u, 0, size);
    memset(v, 0, size);
    memset(u_prev, 0, size);
    memset(v_prev, 0, size);
    memset(p, 0, size);
    memset(div_v, 0, size);

    // Add initial random velocities in a central region to start the simulation
    int center_x = N_x / 2;
    int center_y = N_y / 4;
    int region_size = 5;
    for (int i = -region_size; i <= region_size; i++) {
        for (int j = -region_size; j <= region_size; j++) {
            if (center_x + i > 0 && center_x + i < N_x -1 && center_y + j > 0 && center_y + j < N_y - 1) {
                u[IX(center_x + i, center_y + j)] = (mt_rand() / (float)UINT32_MAX - 0.5f) * 50.0f;
                v[IX(center_x + i, center_y + j)] = (mt_rand() / (float)UINT32_MAX) * 100.0f;
            }
        }
    }
}

void run_computation() {
    float dt = 0.2f; // Simulation time step length

    for (int t = 0; t < num_time_steps; ++t) {
        // Diffuse velocity fields
        SWAP(u_prev, u);
        diffuse(1, u, u_prev, viscosity, dt);
        SWAP(v_prev, v);
        diffuse(2, v, v_prev, viscosity, dt);

        // Project to enforce incompressibility
        project(u, v, p, div_v);

        // Advect velocity fields through themselves
        SWAP(u_prev, u);
        SWAP(v_prev, v);
        advect(1, u, u_prev, u_prev, v_prev, dt);
        advect(2, v, v_prev, u_prev, v_prev, dt);

        // Project again to clean up any errors
        project(u, v, p, div_v);
    }

    double sum = 0.0;
    for (int i = 1; i <= grid_dim_x; i++) {
        for (int j = 1; j <= grid_dim_y; j++) {
            sum += u[IX(i, j)] + v[IX(i, j)];
        }
    }
    final_result = sum;
}

void cleanup() {
    free(u);
    free(v);
    free(u_prev);
    free(v_prev);
    free(p);
    free(div_v);
}

// --- MAIN AND HELPER IMPLEMENTATIONS ---

void set_bnd(int b, float *x) {
    for (int i = 1; i <= grid_dim_y; i++) {
        x[IX(0, i)]     = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
        x[IX(N_x-1, i)] = b == 1 ? -x[IX(N_x-2, i)] : x[IX(N_x-2, i)];
    }
    for (int i = 1; i <= grid_dim_x; i++) {
        x[IX(i, 0)]     = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N_y-1)] = b == 2 ? -x[IX(i, N_y-2)] : x[IX(i, N_y-2)];
    }
    x[IX(0, 0)]         = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N_y-1)]     = 0.5f * (x[IX(1, N_y-1)] + x[IX(0, N_y-2)]);
    x[IX(N_x-1, 0)]     = 0.5f * (x[IX(N_x-2, 0)] + x[IX(N_x-1, 1)]);
    x[IX(N_x-1, N_y-1)] = 0.5f * (x[IX(N_x-2, N_y-1)] + x[IX(N_x-1, N_y-2)]);
}

void lin_solve(int b, float *x, float *x0, float a, float c) {
    float cRecip = 1.0f / c;
    for (int k = 0; k < 20; k++) {
        for (int j = 1; j <= grid_dim_y; j++) {
            for (int i = 1; i <= grid_dim_x; i++) {
                x[IX(i, j)] = (x0[IX(i, j)] 
                                + a*(x[IX(i+1, j)] + x[IX(i-1, j)] 
                                + x[IX(i, j+1)] + x[IX(i, j-1)])) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}

void diffuse(int b, float *x, float *x0, float diff, float dt) {
    float a = dt * diff * (grid_dim_x * grid_dim_y);
    lin_solve(b, x, x0, a, 1 + 4 * a);
}

void advect(int b, float *d, float *d0, float *u, float *v, float dt) {
    float dt_x = dt * (N_x - 2);
    float dt_y = dt * (N_y - 2);
    for (int j = 1; j <= grid_dim_y; j++) {
        for (int i = 1; i <= grid_dim_x; i++) {
            float x = i - dt_x * u[IX(i, j)];
            float y = j - dt_y * v[IX(i, j)];
            if (x < 0.5f) x = 0.5f; if (x > (N_x-2) + 0.5f) x = (N_x-2) + 0.5f;
            int i0 = (int)x; int i1 = i0 + 1;
            if (y < 0.5f) y = 0.5f; if (y > (N_y-2) + 0.5f) y = (N_y-2) + 0.5f;
            int j0 = (int)y; int j1 = j0 + 1;
            float s1 = x - i0; float s0 = 1 - s1;
            float t1 = y - j0; float t0 = 1 - t1;
            d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                          s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
        }
    }
    set_bnd(b, d);
}

void project(float *u, float *v, float *p, float *div) {
    for (int j = 1; j <= grid_dim_y; j++) {
        for (int i = 1; i <= grid_dim_x; i++) {
            div[IX(i, j)] = -0.5f * (u[IX(i+1, j)] - u[IX(i-1, j)] +
                                     v[IX(i, j+1)] - v[IX(i, j-1)]) / grid_dim_x;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 4);
    for (int j = 1; j <= grid_dim_y; j++) {
        for (int i = 1; i <= grid_dim_x; i++) {
            u[IX(i, j)] -= 0.5f * (p[IX(i+1, j)] - p[IX(i-1, j)]) * grid_dim_x;
            v[IX(i, j)] -= 0.5f * (p[IX(i, j+1)] - p[IX(i, j-1)]) * grid_dim_y;
        }
    }
    set_bnd(1, u);
    set_bnd(2, v);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}