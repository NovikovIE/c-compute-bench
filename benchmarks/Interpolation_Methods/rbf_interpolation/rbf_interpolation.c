#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) --- Start
// Do Not Modify - Include This Verbatim
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
// --- Mersenne Twister --- End

// --- Benchmark Globals --- 
typedef struct {
    int num_points;
    int num_queries;
    int rbf_type; // 0: Gaussian, 1: Multiquadric, 2: Inverse Multiquadric
    double epsilon;

    // Data arrays
    double* points_x;
    double* points_y;
    double* points_z; // Known values at (x,y)
    double* weights;
    double* query_x;
    double* query_y;

    double result_accumulator;
} BenchmarkData;

BenchmarkData g_data;

// Helper to get a random double between 0 and 1
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

double rbf_kernel(double r_sq) {
    double r = sqrt(r_sq);
    double eps_r = g_data.epsilon * r;
    switch (g_data.rbf_type) {
        case 0: // Gaussian
            return exp(-(eps_r * eps_r));
        case 1: // Multiquadric
            return sqrt(1.0 + (eps_r * eps_r));
        case 2: // Inverse Multiquadric
            return 1.0 / sqrt(1.0 + (eps_r * eps_r));
        default: // Fallback to Gaussian
            return exp(-(eps_r * eps_r));
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_points num_queries rbf_type epsilon seed\n", argv[0]);
        exit(1);
    }

    g_data.num_points = atoi(argv[1]);
    g_data.num_queries = atoi(argv[2]);
    g_data.rbf_type = atoi(argv[3]);
    g_data.epsilon = atof(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    g_data.points_x = (double*)malloc(g_data.num_points * sizeof(double));
    g_data.points_y = (double*)malloc(g_data.num_points * sizeof(double));
    g_data.points_z = (double*)malloc(g_data.num_points * sizeof(double));
    g_data.weights = (double*)malloc(g_data.num_points * sizeof(double));
    g_data.query_x = (double*)malloc(g_data.num_queries * sizeof(double));
    g_data.query_y = (double*)malloc(g_data.num_queries * sizeof(double));

    if (!g_data.points_x || !g_data.points_y || !g_data.points_z || !g_data.weights || !g_data.query_x || !g_data.query_y) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate random known points and values
    for (int i = 0; i < g_data.num_points; ++i) {
        g_data.points_x[i] = rand_double();
        g_data.points_y[i] = rand_double();
        g_data.points_z[i] = sin(g_data.points_x[i] * 2 * M_PI) * cos(g_data.points_y[i] * 2 * M_PI);
    }

    // Generate random query points
    for (int i = 0; i < g_data.num_queries; ++i) {
        g_data.query_x[i] = rand_double();
        g_data.query_y[i] = rand_double();
    }

    // --- RBF Weight Calculation --- 
    // Solve the linear system A*w = z to find weights 'w'.
    int n = g_data.num_points;
    double* A = (double*)malloc(n * n * sizeof(double));
    if (!A) { fprintf(stderr, "Matrix allocation failed.\n"); exit(1); }
    
    // Initialize RHS vector b (with z values) and matrix A
    for(int i = 0; i < n; ++i) g_data.weights[i] = g_data.points_z[i];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double dx = g_data.points_x[i] - g_data.points_x[j];
            double dy = g_data.points_y[i] - g_data.points_y[j];
            A[i * n + j] = rbf_kernel(dx * dx + dy * dy);
        }
    }

    // Gaussian elimination with partial pivoting to solve A*w = z
    for (int i = 0; i < n; i++) {
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k * n + i]) > fabs(A[max_row * n + i])) {
                max_row = k;
            }
        }

        for (int k = i; k < n; k++) {
            double temp = A[i * n + k];
            A[i * n + k] = A[max_row * n + k];
            A[max_row * n + k] = temp;
        }
        double temp_b = g_data.weights[i];
        g_data.weights[i] = g_data.weights[max_row];
        g_data.weights[max_row] = temp_b;

        for (int k = i + 1; k < n; k++) {
            double factor = A[k * n + i] / A[i * n + i];
            g_data.weights[k] -= factor * g_data.weights[i];
            for (int j = i; j < n; j++) {
                A[k * n + j] -= factor * A[i * n + j];
            }
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i * n + j] * g_data.weights[j];
        }
        g_data.weights[i] = (g_data.weights[i] - sum) / A[i * n + i];
    }

    free(A);
}

void run_computation() {
    g_data.result_accumulator = 0.0;
    for (int i = 0; i < g_data.num_queries; i++) {
        double qx = g_data.query_x[i];
        double qy = g_data.query_y[i];
        double interpolated_value = 0.0;

        for (int j = 0; j < g_data.num_points; j++) {
            double dx = qx - g_data.points_x[j];
            double dy = qy - g_data.points_y[j];
            double r_sq = dx * dx + dy * dy;
            interpolated_value += g_data.weights[j] * rbf_kernel(r_sq);
        }
        g_data.result_accumulator += interpolated_value;
    }
}

void cleanup() {
    free(g_data.points_x);
    free(g_data.points_y);
    free(g_data.points_z);
    free(g_data.weights);
    free(g_data.query_x);
    free(g_data.query_y);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.result_accumulator);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
