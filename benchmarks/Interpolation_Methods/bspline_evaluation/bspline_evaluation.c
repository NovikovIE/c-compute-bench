#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

/**
 * @file bspline_evaluation.c
 * @brief Benchmark for B-spline curve evaluation using De Boor's algorithm.
 *
 * This program evaluates a B-spline curve at numerous points. The B-spline
 * is defined by a set of control points, a degree, and a knot vector.
 * The core computation uses De Boor's algorithm, which is a numerically stable
 * and efficient method for this task.
 */

// --- Mersenne Twister (MT19937) --- Do Not Modify ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
static int num_control_points;
static int degree;
static int num_eval_points;

static double* control_points; // Y-values of control points
static double* knot_vector;
static double* eval_params;    // t-values for evaluation

static double accumulated_result;

// --- Function Prototypes for helper functions ---
static double evaluate_spline_at(double t);

/**
 * @brief Sets up the benchmark data structures and parameters.
 *
 * Parses command-line arguments, allocates memory for control points,
 * the knot vector, and evaluation parameters. It initializes these
 * structures with random data using the Mersenne Twister generator.
 * An open uniform (clamped) knot vector is created, which is a standard
 * choice for B-spline applications.
 */
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_control_points degree num_eval_points seed\n", argv[0]);
        exit(1);
    }

    num_control_points = atoi(argv[1]);
    degree = atoi(argv[2]);
    num_eval_points = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (num_control_points <= degree) {
        fprintf(stderr, "Error: num_control_points must be greater than degree.\n");
        exit(1);
    }
    if (degree < 1) {
        fprintf(stderr, "Error: degree must be at least 1.\n");
        exit(1);
    }

    mt_seed(seed);

    control_points = (double*)malloc(num_control_points * sizeof(double));
    knot_vector = (double*)malloc((num_control_points + degree + 1) * sizeof(double));
    eval_params = (double*)malloc(num_eval_points * sizeof(double));

    if (!control_points || !knot_vector || !eval_params) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Initialize control points with random values
    for (int i = 0; i < num_control_points; ++i) {
        control_points[i] = (double)mt_rand() / (double)UINT32_MAX * 100.0 - 50.0;
    }

    // Create a clamped (open uniform) knot vector
    for (int i = 0; i <= degree; ++i) {
        knot_vector[i] = 0.0;
    }
    for (int i = degree + 1; i < num_control_points; ++i) {
        knot_vector[i] = (double)(i - degree);
    }
    double max_knot_val = (double)(num_control_points - degree);
    for (int i = num_control_points; i <= num_control_points + degree; ++i) {
        knot_vector[i] = max_knot_val;
    }

    // Initialize evaluation points with random t-values within the valid domain
    for (int i = 0; i < num_eval_points; ++i) {
        eval_params[i] = ((double)mt_rand() / (double)UINT32_MAX) * max_knot_val;
    }

    accumulated_result = 0.0;
}

/**
 * @brief Frees all memory allocated by setup_benchmark.
 */
void cleanup() {
    free(control_points);
    free(knot_vector);
    free(eval_params);
}

/**
 * @brief Finds the knot interval for a given parameter t using binary search.
 * @param t The parameter value.
 * @return The index `k` such that knot_vector[k] <= t < knot_vector[k+1].
 */
static int find_knot_interval(double t) {
    if (t >= knot_vector[num_control_points]) {
        return num_control_points - 1;
    }

    int low = degree;
    int high = num_control_points;
    int mid;

    while (low < high) {
        mid = low + (high - low) / 2;
        if (t >= knot_vector[mid]) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return low - 1;
}

/**
 * @brief Evaluates the B-spline at parameter t using De Boor's algorithm.
 * @param t The parameter value at which to evaluate the spline.
 * @return The evaluated point on the spline.
 */
static double evaluate_spline_at(double t) {
    int k = find_knot_interval(t);

    // VLA is a C99 feature. It avoids heap allocation in the main loop.
    double d[degree + 1];

    for (int j = 0; j <= degree; ++j) {
        d[j] = control_points[j + k - degree];
    }

    for (int r = 1; r <= degree; ++r) {
        for (int j = degree; j >= r; --j) {
            int idx1 = j + k - degree;
            int idx2 = j + 1 + k - r;
            double denom = knot_vector[idx2] - knot_vector[idx1];
            double alpha = (denom == 0.0) ? 0.0 : (t - knot_vector[idx1]) / denom;
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
        }
    }

    return d[degree];
}


/**
 * @brief Executes the core computation of the benchmark.
 *
 * Iterates through all evaluation points, calculates the spline value at each point,
 * and accumulates the results.
 */
void run_computation() {
    double sum = 0.0;
    for (int i = 0; i < num_eval_points; i++) {
        sum += evaluate_spline_at(eval_params[i]);
    }
    accumulated_result = sum;
}

/**
 * @brief Main function to drive the benchmark.
 *
 * It follows the structure: setup, time(computation), cleanup, report.
 */
int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the accumulated result to stdout to prevent dead code elimination
    printf("%f\n", accumulated_result);

    // Print the time taken to stderr, without a newline
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
