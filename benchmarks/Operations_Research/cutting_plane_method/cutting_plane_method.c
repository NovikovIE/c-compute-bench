#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA AND GLOBALS ---
typedef struct {
    int num_variables;
    int num_constraints_initial;
    int max_cutting_planes;
    int current_num_constraints;
    int max_total_constraints;

    double** A; // Constraint matrix: A*x <= b
    double* b;  // Constraint limits
    double* x;  // Solution vector
} ProblemData;

static ProblemData* data;
static int final_result; // Accumulated result to prevent dead code elimination

// --- BENCHMARK FUNCTIONS ---

// Helper to generate a random double between -1.0 and 1.0
double random_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_variables num_constraints max_cutting_planes seed\n", argv[0]);
        exit(1);
    }

    data = (ProblemData*)malloc(sizeof(ProblemData));
    if (!data) {
        perror("Failed to allocate memory for ProblemData");
        exit(1);
    }

    data->num_variables = atoi(argv[1]);
    data->num_constraints_initial = atoi(argv[2]);
    data->max_cutting_planes = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);

    mt_seed(seed);

    data->current_num_constraints = data->num_constraints_initial;
    data->max_total_constraints = data->num_constraints_initial + data->max_cutting_planes;

    // Allocate matrix A for the maximum possible number of constraints
    data->A = (double**)malloc(data->max_total_constraints * sizeof(double*));
    if (!data->A) { perror("malloc failed"); exit(1); }
    for (int i = 0; i < data->max_total_constraints; ++i) {
        data->A[i] = (double*)malloc(data->num_variables * sizeof(double));
        if (!data->A[i]) { perror("malloc failed"); exit(1); }
    }

    // Allocate vectors b and x
    data->b = (double*)malloc(data->max_total_constraints * sizeof(double));
    data->x = (double*)malloc(data->num_variables * sizeof(double));
    if (!data->b || !data->x) { perror("malloc failed"); exit(1); }

    // Initialize the initial problem with random data
    for (int i = 0; i < data->num_constraints_initial; ++i) {
        for (int j = 0; j < data->num_variables; ++j) {
            data->A[i][j] = random_double();
        }
        data->b[i] = random_double() * data->num_variables;
    }

    // Initialize solution vector
    for (int i = 0; i < data->num_variables; ++i) {
        data->x[i] = random_double();
    }
}

void run_computation() {
    double accumulated_sum = 0.0;

    for (int k = 0; k < data->max_cutting_planes; ++k) {
        // 1. Simulate solving the LP relaxation. This is the main workload.
        // We do a few fixed iterations of work that resembles a solver.
        for (int iter = 0; iter < 5; ++iter) {
            for (int i = 0; i < data->current_num_constraints; ++i) {
                double dot_product = 0.0;
                for (int j = 0; j < data->num_variables; ++j) {
                    dot_product += data->A[i][j] * data->x[j];
                }
                // If a constraint is violated, "nudge" the solution vector.
                if (dot_product > data->b[i]) {
                    for (int j = 0; j < data->num_variables; ++j) {
                        data->x[j] -= 0.001 * data->A[i][j];
                    }
                }
            }
        }

        // 2. Simulate finding a fractional variable and adding a cutting plane.
        // In a real solver, a cut is derived from the simplex tableau.
        // Here, we simulate it by adding a new random constraint.
        if (data->current_num_constraints < data->max_total_constraints) {
            int new_idx = data->current_num_constraints;
            for (int j = 0; j < data->num_variables; ++j) {
                data->A[new_idx][j] = random_double();
            }
            data->b[new_idx] = random_double() * data->num_variables;
            data->current_num_constraints++;
        }

        // 3. Accumulate a value to prevent compiler optimization.
        accumulated_sum += data->x[k % data->num_variables];
    }

    final_result = (int)accumulated_sum;
}

void cleanup() {
    for (int i = 0; i < data->max_total_constraints; ++i) {
        free(data->A[i]);
    }
    free(data->A);
    free(data->b);
    free(data->x);
    free(data);
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
    printf("%d\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
