#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Benchmark Globals ---

// Parameters from argv
int num_dimensions;
int num_iterations;

// Nelder-Mead constants
const double ALPHA = 1.0;
const double GAMMA = 2.0;
const double RHO = 0.5;
const double SIGMA = 0.5;

// Global data structures for the benchmark
double** simplex;      // Simplex vertices: (N+1) x N
double* f_values;      // Function value for each vertex: (N+1)
double* centroid;      // Centroid of best N points: N
double* p_refl;        // Reflected point: N
double* p_exp;         // Expanded point: N
double* p_contr;       // Contracted point: N
double* target_min;    // The "true" minimum for the objective function: N
int* indices;          // Indices for sorting points: (N+1)

double final_result; // To prevent dead-code elimination

// --- Objective Function & Helpers ---

// Objective function: Shifted Sphere function
// f(x) = sum_{i=0}^{N-1} (x_i - target_min_i)^2
double objective_function(const double* x) {
    double sum = 0.0;
    for (int i = 0; i < num_dimensions; ++i) {
        double diff = x[i] - target_min[i];
        sum += diff * diff;
    }
    return sum;
}

// qsort comparator for sorting vertex indices based on their function values
int compare_indices(const void* a, const void* b) {
    int idx1 = *(const int*)a;
    int idx2 = *(const int*)b;
    if (f_values[idx1] < f_values[idx2]) return -1;
    if (f_values[idx1] > f_values[idx2]) return 1;
    return 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char** argv) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_dimensions> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }
    num_dimensions = atoi(argv[1]);
    num_iterations = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_dimensions <= 0 || num_iterations <= 0) {
        fprintf(stderr, "FATAL: num_dimensions and num_iterations must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    int n = num_dimensions;
    int n_plus_1 = num_dimensions + 1;

    simplex = (double**)malloc(n_plus_1 * sizeof(double*));
    for (int i = 0; i < n_plus_1; ++i) {
        simplex[i] = (double*)malloc(n * sizeof(double));
    }
    f_values = (double*)malloc(n_plus_1 * sizeof(double));
    indices = (int*)malloc(n_plus_1 * sizeof(int));

    centroid = (double*)malloc(n * sizeof(double));
    p_refl = (double*)malloc(n * sizeof(double));
    p_exp = (double*)malloc(n * sizeof(double));
    p_contr = (double*)malloc(n * sizeof(double));
    target_min = (double*)malloc(n * sizeof(double));

    // Initialize the target minimum for the objective function
    for (int i = 0; i < n; ++i) {
        target_min[i] = ((double)mt_rand() / 4294967296.0) * 20.0 - 10.0; // Random in [-10, 10]
    }

    // Initialize the simplex with random points
    for (int i = 0; i < n_plus_1; ++i) {
        for (int j = 0; j < n; ++j) {
            simplex[i][j] = ((double)mt_rand() / 4294967296.0) * 10.0 - 5.0; // Random in [-5, 5]
        }
    }
}

void run_computation() {
    int n = num_dimensions;
    int n_plus_1 = n + 1;

    for (int iter = 0; iter < num_iterations; ++iter) {
        // 1. Evaluate objective function for all points and set up initial indices
        for (int i = 0; i < n_plus_1; ++i) {
            f_values[i] = objective_function(simplex[i]);
            indices[i] = i;
        }

        // 2. Sort points by function value (best to worst)
        qsort(indices, n_plus_1, sizeof(int), compare_indices);
        int best_idx = indices[0];
        int second_worst_idx = indices[n - 1];
        int worst_idx = indices[n];

        // 3. Calculate centroid of the best n points
        for (int j = 0; j < n; ++j) centroid[j] = 0.0;
        for (int i = 0; i < n; ++i) { // Excludes the worst point
            for (int j = 0; j < n; ++j) {
                centroid[j] += simplex[indices[i]][j];
            }
        }
        for (int j = 0; j < n; ++j) centroid[j] /= n;

        // 4. Reflection
        for (int j = 0; j < n; ++j) {
            p_refl[j] = centroid[j] + ALPHA * (centroid[j] - simplex[worst_idx][j]);
        }
        double f_refl = objective_function(p_refl);

        if (f_values[best_idx] <= f_refl && f_refl < f_values[second_worst_idx]) {
            // Accept reflected point
            for (int j = 0; j < n; ++j) simplex[worst_idx][j] = p_refl[j];
        }
        // 5. Expansion
        else if (f_refl < f_values[best_idx]) {
            for (int j = 0; j < n; ++j) {
                p_exp[j] = centroid[j] + GAMMA * (p_refl[j] - centroid[j]);
            }
            double f_exp = objective_function(p_exp);
            for (int j = 0; j < n; ++j) {
                simplex[worst_idx][j] = (f_exp < f_refl) ? p_exp[j] : p_refl[j];
            }
        }
        // 6. Contraction
        else {
            double f_contr;
            int perform_shrink = 0;
            if (f_refl < f_values[worst_idx]) {
                // Outside contraction
                for (int j = 0; j < n; ++j) p_contr[j] = centroid[j] + RHO * (p_refl[j] - centroid[j]);
                f_contr = objective_function(p_contr);
                if (f_contr <= f_refl) {
                    for (int j = 0; j < n; ++j) simplex[worst_idx][j] = p_contr[j];
                } else {
                    perform_shrink = 1;
                }
            } else {
                // Inside contraction
                for (int j = 0; j < n; ++j) p_contr[j] = centroid[j] - RHO * (simplex[worst_idx][j] - centroid[j]);
                f_contr = objective_function(p_contr);
                if (f_contr < f_values[worst_idx]) {
                    for (int j = 0; j < n; ++j) simplex[worst_idx][j] = p_contr[j];
                } else {
                    perform_shrink = 1;
                }
            }
            // 7. Shrink
            if (perform_shrink) {
                for (int i = 1; i < n_plus_1; ++i) {
                    int current_idx = indices[i];
                    for (int j = 0; j < n; ++j) {
                        simplex[current_idx][j] = simplex[best_idx][j] + SIGMA * (simplex[current_idx][j] - simplex[best_idx][j]);
                    }
                }
            }
        }
    }

    // To determine the final result, re-evaluate and find the best point's coordinate sum.
    for (int i = 0; i < n_plus_1; ++i) {
        f_values[i] = objective_function(simplex[i]);
        indices[i] = i;
    }
    qsort(indices, n_plus_1, sizeof(int), compare_indices);
    int best_idx = indices[0];

    double sum = 0.0;
    for (int j = 0; j < num_dimensions; ++j) {
        sum += simplex[best_idx][j];
    }
    final_result = sum;
}

void cleanup() {
    int n_plus_1 = num_dimensions + 1;
    for (int i = 0; i < n_plus_1; ++i) {
        free(simplex[i]);
    }
    free(simplex);
    free(f_values);
    free(indices);
    free(centroid);
    free(p_refl);
    free(p_exp);
    free(p_contr);
    free(target_min);
}

// --- Main Execution ---

int main(int argc, char** argv) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Final result to stdout
    printf("%f\n", final_result);

    // Time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
