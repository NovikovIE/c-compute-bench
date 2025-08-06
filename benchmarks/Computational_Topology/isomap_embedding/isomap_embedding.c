#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

// --- Mersenne Twister (verbatim) ---
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

// Benchmark parameters
int num_points;
int num_neighbors;
int target_dimension; // Parameter is present but not used in this simplified computation

// Data structures
double **points;         // num_points x 3 (original high-dim data - Swiss Roll)
double **dist_matrix;    // num_points x num_points (geodesic distances)
double final_result_sum; // To prevent dead code elimination

// Helper for qsort
int compare_doubles(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

// Function to allocate a 2D double matrix
double** allocate_matrix(int rows, int cols) {
    double **mat = (double **)malloc(rows * sizeof(double *));
    if (!mat) return NULL;
    for (int i = 0; i < rows; i++) {
        mat[i] = (double *)malloc(cols * sizeof(double));
        if (!mat[i]) {
            for (int j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
    }
    return mat;
}

// Function to free a 2D double matrix
void free_matrix(double **mat, int rows) {
    if (!mat) return;
    for (int i = 0; i < rows; i++) {
        free(mat[i]);
    }
    free(mat);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_points num_neighbors target_dimension seed\n", argv[0]);
        exit(1);
    }
    num_points = atoi(argv[1]);
    num_neighbors = atoi(argv[2]);
    target_dimension = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (num_points <= 0 || num_neighbors <= 0 || target_dimension <= 0 || num_neighbors >= num_points) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }
    
    mt_seed(seed);

    points = allocate_matrix(num_points, 3);
    dist_matrix = allocate_matrix(num_points, num_points);

    if (!points || !dist_matrix){
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate points on a Swiss Roll manifold
    for (int i = 0; i < num_points; i++) {
        double t = 1.5 * M_PI * (1.0 + 2.0 * (mt_rand() / (double)UINT32_MAX));
        double y = 21.0 * (mt_rand() / (double)UINT32_MAX);
        points[i][0] = t * cos(t);
        points[i][1] = y;
        points[i][2] = t * sin(t);
    }
}

void run_computation() {
    double** euclidean_dists = allocate_matrix(num_points, num_points);
    if (!euclidean_dists) {
        fprintf(stderr, "Failed to allocate euclidean distance matrix\n");
        exit(1);
    }

    // 1. Compute pairwise Euclidean distances
    for (int i = 0; i < num_points; i++) {
        for (int j = i; j < num_points; j++) {
            if (i == j) {
                euclidean_dists[i][j] = 0.0;
            } else {
                double dx = points[i][0] - points[j][0];
                double dy = points[i][1] - points[j][1];
                double dz = points[i][2] - points[j][2];
                double dist = sqrt(dx*dx + dy*dy + dz*dz);
                euclidean_dists[i][j] = dist;
                euclidean_dists[j][i] = dist;
            }
        }
    }

    // 2. Construct neighborhood graph (adjacency matrix stored in dist_matrix)
    double* sorted_dists = (double*)malloc(num_points * sizeof(double));
    if (!sorted_dists) {
        fprintf(stderr, "Failed to allocate sorted distances array\n");
        exit(1);
    }
    for (int i = 0; i < num_points; i++) {
        for(int j=0; j < num_points; j++) sorted_dists[j] = euclidean_dists[i][j];
        qsort(sorted_dists, num_points, sizeof(double), compare_doubles);
        double threshold = sorted_dists[num_neighbors];

        for (int j = 0; j < num_points; j++) {
            if (i == j) {
                dist_matrix[i][j] = 0.0;
            } else if (euclidean_dists[i][j] <= threshold) {
                dist_matrix[i][j] = euclidean_dists[i][j];
            } else {
                dist_matrix[i][j] = DBL_MAX;
            }
        }
    }
    free(sorted_dists);
    free_matrix(euclidean_dists, num_points);

    // 3. Compute all-pairs shortest paths using Floyd-Warshall
    for (int k = 0; k < num_points; k++) {
        for (int i = 0; i < num_points; i++) {
            for (int j = 0; j < num_points; j++) {
                if (dist_matrix[i][k] != DBL_MAX && dist_matrix[k][j] != DBL_MAX) {
                    double new_dist = dist_matrix[i][k] + dist_matrix[k][j];
                    if (new_dist < dist_matrix[i][j]) {
                        dist_matrix[i][j] = new_dist;
                    }
                }
            }
        }
    }

    // 4. Calculate a final result to prevent dead-code elimination.
    double sum = 0.0;
    for (int i = 0; i < num_points; i++) {
        for (int j = 0; j < num_points; j++) {
            if (dist_matrix[i][j] != DBL_MAX) {
                sum += dist_matrix[i][j];
            }
        }
    }
    final_result_sum = sum;
}

void cleanup() {
    free_matrix(points, num_points);
    free_matrix(dist_matrix, num_points);
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
    printf("%f\n", final_result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
