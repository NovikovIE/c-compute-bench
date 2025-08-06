#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

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
// --- End of Mersenne Twister ---

// Benchmark parameters
int num_control_points;
int num_queries;

// Global data structures for the benchmark
double *control_x;      // X coordinates of control points
double *control_y;      // Y coordinates of control points
double *query_x;        // X coordinates of query points
double *query_y;        // Y coordinates of query points
double *tps_coeffs;      // Solved TPS coefficients (N_control weights + 3 affine)
double final_result;    // Accumulator for results

// Generates a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / UINT32_MAX;
}

// Thin Plate Spline radial basis function: U(r) = r^2 * log(r)
double tps_kernel(double r) {
    if (r < 1e-9) { // Avoid log(0)
        return 0.0;
    }
    return r * r * log(r);
}

// Solves the linear system Ax=b for x using Gaussian elimination with partial pivoting.
// A is modified in place. The solution x is stored in b.
void solve_linear_system(double **A, double *b, int n) {
    for (int i = 0; i < n; i++) {
        // Find pivot row
        int max = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[max][i])) {
                max = k;
            }
        }

        // Swap rows in A and b
        double *temp_row = A[i]; A[i] = A[max]; A[max] = temp_row;
        double temp_b = b[i]; b[i] = b[max]; b[max] = temp_b;

        if (fabs(A[i][i]) < 1e-12) {
            fprintf(stderr, "FATAL: Singular or near-singular matrix.\n");
            exit(1);
        }

        // Forward elimination
        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            b[k] -= factor * b[i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * b[i + (j - i)]; // Using b as solution vector
        }
        b[i] = (b[i] - sum) / A[i][i];
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_control_points num_queries seed\n", argv[0]);
        exit(1);
    }

    num_control_points = atoi(argv[1]);
    num_queries = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    // Allocate global arrays
    control_x = (double*)malloc(num_control_points * sizeof(double));
    control_y = (double*)malloc(num_control_points * sizeof(double));
    double *control_z = (double*)malloc(num_control_points * sizeof(double));
    query_x = (double*)malloc(num_queries * sizeof(double));
    query_y = (double*)malloc(num_queries * sizeof(double));
    
    int system_size = num_control_points + 3;
    tps_coeffs = (double*)malloc(system_size * sizeof(double));

    // Generate random control and query points
    for (int i = 0; i < num_control_points; ++i) {
        control_x[i] = rand_double();
        control_y[i] = rand_double();
        control_z[i] = rand_double();
    }
    for (int i = 0; i < num_queries; ++i) {
        query_x[i] = rand_double();
        query_y[i] = rand_double();
    }

    // --- Setup and solve the linear system for TPS coefficients ---
    // Matrix L and vector V for the system L*w = V
    double **L = (double**)malloc(system_size * sizeof(double*));
    for(int i = 0; i < system_size; ++i) {
        L[i] = (double*)malloc(system_size * sizeof(double));
        memset(L[i], 0, system_size * sizeof(double));
    }
    double *V = tps_coeffs; // Use tps_coeffs as a temporary V vector, then it becomes the solution

    // 1. Fill K part of L (upper-left NxN submatrix)
    for (int i = 0; i < num_control_points; ++i) {
        for (int j = i; j < num_control_points; ++j) {
            double dx = control_x[i] - control_x[j];
            double dy = control_y[i] - control_y[j];
            double dist = sqrt(dx * dx + dy * dy);
            double val = tps_kernel(dist);
            L[i][j] = val;
            L[j][i] = val;
        }
    }

    // 2. Fill P and P' parts of L
    for (int i = 0; i < num_control_points; ++i) {
        // P part (upper-right)
        L[i][num_control_points] = 1.0;
        L[i][num_control_points + 1] = control_x[i];
        L[i][num_control_points + 2] = control_y[i];
        
        // P' part (lower-left)
        L[num_control_points][i] = 1.0;
        L[num_control_points + 1][i] = control_x[i];
        L[num_control_points + 2][i] = control_y[i];
    }

    // 3. Fill V vector
    for (int i = 0; i < num_control_points; ++i) {
        V[i] = control_z[i];
    }
    V[num_control_points] = 0.0;
    V[num_control_points + 1] = 0.0;
    V[num_control_points + 2] = 0.0;

    // 4. Solve L * w = V. Solution w is stored in V (tps_coeffs).
    solve_linear_system(L, V, system_size);

    // Cleanup temporary structures
    free(control_z);
    for(int i = 0; i < system_size; ++i) {
        free(L[i]);
    }
    free(L);
}

void run_computation() {
    final_result = 0.0;
    double total_interpolated_value = 0.0;

    double *weights = tps_coeffs;
    double a1 = tps_coeffs[num_control_points];
    double a2 = tps_coeffs[num_control_points + 1];
    double a3 = tps_coeffs[num_control_points + 2];

    for (int i = 0; i < num_queries; ++i) {
        double qx = query_x[i];
        double qy = query_y[i];

        // Affine part
        double interpolated_value = a1 + a2 * qx + a3 * qy;

        // Radial basis part
        for (int j = 0; j < num_control_points; ++j) {
            double dx = qx - control_x[j];
            double dy = qy - control_y[j];
            double dist = sqrt(dx * dx + dy * dy);
            interpolated_value += weights[j] * tps_kernel(dist);
        }
        total_interpolated_value += interpolated_value;
    }
    final_result = total_interpolated_value;
}

void cleanup() {
    free(control_x);
    free(control_y);
    free(query_x);
    free(query_y);
    free(tps_coeffs);
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
