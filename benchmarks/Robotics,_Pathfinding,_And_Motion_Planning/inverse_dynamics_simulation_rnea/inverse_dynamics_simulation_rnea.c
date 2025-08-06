#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Data Structures ---
typedef struct {
    int num_links;
    int num_trajectory_points;

    // Robot Model (constants)
    int* parent;
    double** joint_axis;
    double** link_com;
    double* link_mass;
    double*** link_inertia;

    // Trajectory (input)
    double** q;
    double** qd;
    double** qdd;

    // Trajectory (output)
    double** tau;

    // Intermediate/temporary variables for RNEA
    double** w;
    double** wd;
    double** vd;
    double** vdc;
    double** f;
    double** n;

    // Final result accumulator
    double result_accumulator;

} BenchmarkData;

BenchmarkData g_data;

// --- Utility and Memory Functions ---
double rand_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

double** allocate_2d(int rows, int cols) {
    double** arr = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        arr[i] = (double*)calloc(cols, sizeof(double));
    }
    return arr;
}

double*** allocate_3d(int d1, int d2, int d3) {
    double*** arr = (double***)malloc(d1 * sizeof(double**));
    for (int i = 0; i < d1; i++) {
        arr[i] = (double**)malloc(d2 * sizeof(double*));
        for (int j = 0; j < d2; j++) {
            arr[i][j] = (double*)calloc(d3, sizeof(double));
        }
    }
    return arr;
}

void free_2d(double** arr, int rows) {
    for (int i = 0; i < rows; i++) {
        free(arr[i]);
    }
    free(arr);
}

void free_3d(double*** arr, int d1, int d2) {
    for (int i = 0; i < d1; i++) {
        for (int j = 0; j < d2; j++) {
            free(arr[i][j]);
        }
        free(arr[i]);
    }
    free(arr);
}

// --- 3D Vector/Matrix Math Helpers ---
void cross_product(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

void mat_vec_mul(const double M[3][3], const double v[3], double out[3]) {
    out[0] = M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2];
    out[1] = M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2];
    out[2] = M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2];
}

void mat_transpose_vec_mul(const double M[3][3], const double v[3], double out[3]) {
    out[0] = M[0][0] * v[0] + M[1][0] * v[1] + M[2][0] * v[2];
    out[1] = M[0][1] * v[0] + M[1][1] * v[1] + M[2][1] * v[2];
    out[2] = M[0][2] * v[0] + M[1][2] * v[1] + M[2][2] * v[2];
}

void vec_add(const double a[3], const double b[3], double out[3]) {
    out[0] = a[0] + b[0];
    out[1] = a[1] + b[1];
    out[2] = a[2] + b[2];
}

void vec_sub(const double a[3], const double b[3], double out[3]) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

double vec_dot(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

void vec_scale(const double v[3], double s, double out[3]) {
    out[0] = v[0] * s;
    out[1] = v[1] * s;
    out[2] = v[2] * s;
}

void angle_axis_to_rot_mat(const double axis[3], double angle, double R[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    double t = 1.0 - c;
    double x = axis[0], y = axis[1], z = axis[2];
    
    R[0][0] = t*x*x + c;   R[0][1] = t*x*y - s*z; R[0][2] = t*x*z + s*y;
    R[1][0] = t*x*y + s*z; R[1][1] = t*y*y + c;   R[1][2] = t*y*z - s*x;
    R[2][0] = t*x*z - s*y; R[2][1] = t*y*z + s*x; R[2][2] = t*z*z + c;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_links> <num_trajectory_points> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_links = atoi(argv[1]);
    g_data.num_trajectory_points = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);
    mt_seed(seed);

    int nl = g_data.num_links;
    int nt = g_data.num_trajectory_points;

    // Allocate memory
    g_data.parent = (int*)malloc(nl * sizeof(int));
    g_data.link_mass = (double*)malloc(nl * sizeof(double));
    g_data.joint_axis = allocate_2d(nl, 3);
    g_data.link_com = allocate_2d(nl, 3);
    g_data.link_inertia = allocate_3d(nl, 3, 3);
    g_data.q = allocate_2d(nt, nl);
    g_data.qd = allocate_2d(nt, nl);
    g_data.qdd = allocate_2d(nt, nl);
    g_data.tau = allocate_2d(nt, nl);
    
    g_data.w = allocate_2d(nl, 3);
    g_data.wd = allocate_2d(nl, 3);
    g_data.vd = allocate_2d(nl, 3);
    g_data.vdc = allocate_2d(nl, 3);
    g_data.f = allocate_2d(nl, 3);
    g_data.n = allocate_2d(nl, 3);

    // Generate robot model data
    for (int i = 0; i < nl; ++i) {
        g_data.parent[i] = i - 1; // Simple serial chain
        g_data.link_mass[i] = rand_double(1.0, 5.0);

        double norm = 0.0;
        for(int j=0; j<3; ++j) {
            g_data.joint_axis[i][j] = rand_double(-1.0, 1.0);
            norm += g_data.joint_axis[i][j] * g_data.joint_axis[i][j];
        }
        norm = sqrt(norm);
        for(int j=0; j<3; ++j) g_data.joint_axis[i][j] /= norm;

        for(int j=0; j<3; ++j) g_data.link_com[i][j] = rand_double(-0.1, 0.1);

        for(int j=0; j<3; ++j) {
            for(int k=0; k<3; ++k) {
                g_data.link_inertia[i][j][k] = (j == k) ? rand_double(0.01, 0.1) : 0.0;
            }
        }
    }

    // Generate trajectory data
    for (int t = 0; t < nt; ++t) {
        for (int i = 0; i < nl; ++i) {
            g_data.q[t][i] = rand_double(-M_PI, M_PI);
            g_data.qd[t][i] = rand_double(-2.0, 2.0);
            g_data.qdd[t][i] = rand_double(-5.0, 5.0);
        }
    }
    
    g_data.result_accumulator = 0.0;
}

void run_computation() {
    int nl = g_data.num_links;
    int nt = g_data.num_trajectory_points;
    double grav[3] = {0.0, 0.0, -9.81};

    double temp_v1[3], temp_v2[3], temp_v3[3];
    double R[3][3];
    
    for (int t = 0; t < nt; ++t) {
        double w_p[3] = {0,0,0}, wd_p[3] = {0,0,0}, vd_p[3];
        vec_scale(grav, -1.0, vd_p); // Base acceleration opposes gravity

        // Forward pass: compute velocities and accelerations
        for (int i = 0; i < nl; ++i) {
            angle_axis_to_rot_mat(g_data.joint_axis[i], g_data.q[t][i], R);

            mat_transpose_vec_mul(R, w_p, temp_v1);
            vec_scale(g_data.joint_axis[i], g_data.qd[t][i], temp_v2);
            vec_add(temp_v1, temp_v2, g_data.w[i]);

            mat_transpose_vec_mul(R, wd_p, temp_v1);
            cross_product(g_data.w[i], temp_v2, temp_v3);
            vec_add(temp_v1, temp_v3, g_data.wd[i]);
            vec_scale(g_data.joint_axis[i], g_data.qdd[t][i], temp_v2);
            vec_add(g_data.wd[i], temp_v2, g_data.wd[i]);

            mat_transpose_vec_mul(R, vd_p, g_data.vd[i]);

            cross_product(g_data.wd[i], g_data.link_com[i], temp_v1);
            cross_product(g_data.w[i], g_data.link_com[i], temp_v2);
            cross_product(g_data.w[i], temp_v2, temp_v3);
            vec_add(g_data.vd[i], temp_v1, g_data.vdc[i]);
            vec_add(g_data.vdc[i], temp_v3, g_data.vdc[i]);
            
            for(int j=0; j<3; ++j) {
                w_p[j] = g_data.w[i][j];
                wd_p[j] = g_data.wd[i][j];
                vd_p[j] = g_data.vd[i][j];
            }
        }

        // Backward pass: compute forces and torques
        double f_child[3] = {0,0,0}, n_child[3] = {0,0,0};
        for (int i = nl - 1; i >= 0; --i) {
            // Inertial force and torque on the link
            vec_scale(g_data.vdc[i], g_data.link_mass[i], temp_v1); // F_i = m_i * a_c_i
            
            double I_mat[3][3];
            for(int r=0; r<3; ++r) for(int c=0; c<3; ++c) I_mat[r][c] = g_data.link_inertia[i][r][c];
            
            double Iw[3], w_x_Iw[3];
            mat_vec_mul(I_mat, g_data.w[i], Iw);
            cross_product(g_data.w[i], Iw, w_x_Iw);
            
            double Iwd[3];
            mat_vec_mul(I_mat, g_data.wd[i], Iwd);
            vec_add(Iwd, w_x_Iw, temp_v2); // N_i = I_i * wd_i + w_i x (I_i * w_i)

            vec_add(temp_v1, f_child, g_data.f[i]);
            
            cross_product(g_data.link_com[i], g_data.f[i], temp_v3);
            vec_add(n_child, temp_v2, g_data.n[i]);
            vec_add(g_data.n[i], temp_v3, g_data.n[i]);
            
            g_data.tau[t][i] = vec_dot(g_data.n[i], g_data.joint_axis[i]);
            
            if (i > 0) {
              angle_axis_to_rot_mat(g_data.joint_axis[i], g_data.q[t][i], R);
              mat_vec_mul(R, g_data.f[i], f_child);
              mat_vec_mul(R, g_data.n[i], n_child);
            }
        }
        g_data.result_accumulator += g_data.tau[t][0];
    }
}

void cleanup() {
    int nl = g_data.num_links;
    int nt = g_data.num_trajectory_points;
    
    free(g_data.parent);
    free(g_data.link_mass);
    free_2d(g_data.joint_axis, nl);
    free_2d(g_data.link_com, nl);
    free_3d(g_data.link_inertia, nl, 3);
    
    free_2d(g_data.q, nt);
    free_2d(g_data.qd, nt);
    free_2d(g_data.qdd, nt);
    free_2d(g_data.tau, nt);
    
    free_2d(g_data.w, nl);
    free_2d(g_data.wd, nl);
    free_2d(g_data.vd, nl);
    free_2d(g_data.vdc, nl);
    free_2d(g_data.f, nl);
    free_2d(g_data.n, nl);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    
    setup_benchmark(argc, argv);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%f\n", g_data.result_accumulator);
    fprintf(stderr, "%.6f", time_taken);
    
    cleanup();
    
    return 0;
}
