#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- START: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// --- START: Benchmark Globals and Structs ---
typedef struct {
    double x;
    double y;
    double value;
} KnownPoint;

typedef struct {
    // Parameters
    int num_known_points;
    int output_grid_width;
    int output_grid_height;
    double nugget;
    double sill;
    double range;

    // Data
    KnownPoint *known_points;
    double *output_grid;
    double *lhs_matrix_inv; // The inverted (N+1)x(N+1) matrix

    // Result
    double accumulated_result;
} BenchmarkData;

static BenchmarkData g_data;

// --- END: Benchmark Globals and Structs ---

// --- START: Helper Functions ---

double random_double() {
    return mt_rand() / (double)UINT32_MAX;
}

double euclidean_distance(double x1, double y1, double x2, double y2) {
    double dx = x1 - x2;
    double dy = y1 - y2;
    return sqrt(dx * dx + dy * dy);
}

// Spherical semivariogram model covariance function C(h) = sill - gamma(h)
double covariance_model(double h, double nugget, double sill, double range) {
    if (h == 0.0) return sill;
    if (h > range) return 0.0; // sill - sill = 0
    double variogram_part = (sill - nugget) * (1.5 * (h / range) - 0.5 * pow(h / range, 3.0));
    return sill - (nugget + variogram_part);
}

// In-place matrix inversion using LU decomposition with partial pivoting.
// Returns 0 on success, -1 on failure (singular matrix).
// The input matrix A is destroyed and replaced with its inverse.
static int matrix_invert(double *A, int n) {
    int *pivot = (int *)malloc(n * sizeof(int));
    if (!pivot) return -1;

    // LU decomposition with partial pivoting
    for (int i = 0; i < n; i++) {
        pivot[i] = i;
    }
    for (int i = 0; i < n; i++) {
        int max_j = i;
        for (int j = i + 1; j < n; j++) {
            if (fabs(A[j * n + i]) > fabs(A[max_j * n + i])) {
                max_j = j;
            }
        }
        if (max_j != i) {
            for (int k = 0; k < n; k++) {
                double tmp = A[i * n + k];
                A[i * n + k] = A[max_j * n + k];
                A[max_j * n + k] = tmp;
            }
            int tmp_p = pivot[i];
            pivot[i] = pivot[max_j];
            pivot[max_j] = tmp_p;
        }
        if (fabs(A[i * n + i]) < 1e-12) { // Singular matrix?
            free(pivot);
            return -1;
        }
        for (int j = i + 1; j < n; j++) {
            A[j * n + i] /= A[i * n + i];
            for (int k = i + 1; k < n; k++) {
                A[j * n + k] -= A[j * n + i] * A[i * n + k];
            }
        }
    }

    // LU inversion
    double *invA = (double*)malloc(n * n * sizeof(double));
    if (!invA) {
        free(pivot);
        return -1;
    }

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            int p_i = pivot[i];
            invA[i*n+j] = (p_i == j) ? 1.0 : 0.0;
            for (int k = 0; k < i; k++) {
                invA[i*n+j] -= A[i*n+k] * invA[k*n+j];
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            for (int k = i + 1; k < n; k++) {
                invA[i*n+j] -= A[i*n+k] * invA[k*n+j];
            }
            invA[i*n+j] /= A[i*n+i];
        }
    }
    memcpy(A, invA, n * n * sizeof(double));
    free(invA);
    free(pivot);
    return 0;
}

// --- END: Helper Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 8) {
        fprintf(stderr, "Usage: %s num_known_points output_grid_width output_grid_height nugget sill range seed\n", argv[0]);
        exit(1);
    }

    g_data.num_known_points = atoi(argv[1]);
    g_data.output_grid_width = atoi(argv[2]);
    g_data.output_grid_height = atoi(argv[3]);
    g_data.nugget = atof(argv[4]);
    g_data.sill = atof(argv[5]);
    g_data.range = atof(argv[6]);
    uint32_t seed = (uint32_t)atoi(argv[7]);
    mt_seed(seed);

    g_data.known_points = (KnownPoint *)malloc(g_data.num_known_points * sizeof(KnownPoint));
    g_data.output_grid = (double *)malloc(g_data.output_grid_width * g_data.output_grid_height * sizeof(double));
    
    int matrix_size = g_data.num_known_points + 1;
    g_data.lhs_matrix_inv = (double *)malloc(matrix_size * matrix_size * sizeof(double));

    if (!g_data.known_points || !g_data.output_grid || !g_data.lhs_matrix_inv) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_known_points; ++i) {
        g_data.known_points[i].x = random_double();
        g_data.known_points[i].y = random_double();
        g_data.known_points[i].value = random_double() * 100.0;
    }

    // Build the LHS matrix for Kriging system: A * w = b
    double *lhs_matrix = g_data.lhs_matrix_inv; // Use the final destination buffer temporarily
    for (int i = 0; i < g_data.num_known_points; ++i) {
        for (int j = 0; j < g_data.num_known_points; ++j) {
            double dist = euclidean_distance(g_data.known_points[i].x, g_data.known_points[i].y,
                                               g_data.known_points[j].x, g_data.known_points[j].y);
            lhs_matrix[i * matrix_size + j] = covariance_model(dist, g_data.nugget, g_data.sill, g_data.range);
        }
    }

    // Add Lagrange multiplier row/column for unbiasedness constraint
    for (int i = 0; i < g_data.num_known_points; ++i) {
        lhs_matrix[i * matrix_size + g_data.num_known_points] = 1.0;
        lhs_matrix[g_data.num_known_points * matrix_size + i] = 1.0;
    }
    lhs_matrix[matrix_size * matrix_size - 1] = 0.0;

    // Invert the LHS matrix
    if (matrix_invert(lhs_matrix, matrix_size) != 0) {
        fprintf(stderr, "Matrix inversion failed. The matrix might be singular.\n");
        exit(1);
    }
}

void run_computation() {
    g_data.accumulated_result = 0.0;
    int n = g_data.num_known_points;
    int matrix_size = n + 1;

    double *rhs_vector = (double *)malloc(matrix_size * sizeof(double));
    double *weights = (double *)malloc(matrix_size * sizeof(double));
    if (!rhs_vector || !weights) {
        fprintf(stderr, "Memory allocation failed in computation\n");
        exit(1);
    }

    for (int gy = 0; gy < g_data.output_grid_height; ++gy) {
        for (int gx = 0; gx < g_data.output_grid_width; ++gx) {
            double px = (double)gx / (g_data.output_grid_width > 1 ? g_data.output_grid_width - 1 : 1);
            double py = (double)gy / (g_data.output_grid_height > 1 ? g_data.output_grid_height - 1 : 1);

            // Build the RHS vector
            for (int i = 0; i < n; ++i) {
                double dist = euclidean_distance(px, py, g_data.known_points[i].x, g_data.known_points[i].y);
                rhs_vector[i] = covariance_model(dist, g_data.nugget, g_data.sill, g_data.range);
            }
            rhs_vector[n] = 1.0;

            // Solve for weights: w = A_inv * b
            for (int i = 0; i < matrix_size; ++i) {
                weights[i] = 0.0;
                for (int j = 0; j < matrix_size; ++j) {
                    weights[i] += g_data.lhs_matrix_inv[i * matrix_size + j] * rhs_vector[j];
                }
            }

            // Estimate value
            double estimated_value = 0.0;
            for (int i = 0; i < n; ++i) {
                estimated_value += weights[i] * g_data.known_points[i].value;
            }

            g_data.output_grid[gy * g_data.output_grid_width + gx] = estimated_value;
            g_data.accumulated_result += estimated_value;
        }
    }

    free(rhs_vector);
    free(weights);
}

void cleanup() {
    free(g_data.known_points);
    free(g_data.output_grid);
    free(g_data.lhs_matrix_inv);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.accumulated_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
