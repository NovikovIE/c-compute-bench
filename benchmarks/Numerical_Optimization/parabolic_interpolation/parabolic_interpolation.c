#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim as Required) ---
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

// --- Benchmark Globals and Setup ---

// Benchmark parameters
static int NUM_ITERATIONS;
static int MAX_STEPS_PER_PROBLEM; // Controlled by 'precision' argument

// Data structure for a single optimization problem
typedef struct {
    double x1, x2, x3; // Initial 3-point bracket for the minimum
} ProblemData;

static ProblemData* problems;
static double final_result_sum; // To prevent dead code elimination

// The unimodal function to be minimized by the algorithm.
// A function with a non-trivial minimum: f(x) = (x - 2)^2 + 5*cos(x)
// Minimum is near x = 1.33
static inline double objective_function(double x) {
    double term1 = x - 2.0;
    return term1 * term1 + 5.0 * cos(x);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_iterations> <precision> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_ITERATIONS = atoi(argv[1]);
    MAX_STEPS_PER_PROBLEM = atoi(argv[2]); // 'precision' argument sets the number of steps
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    problems = (ProblemData*)malloc(NUM_ITERATIONS * sizeof(ProblemData));
    if (problems == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    const double approx_min = 1.33; // Approximate minimum of objective_function

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        // Generate a random center for the initial bracket, near the true minimum
        double center = approx_min + ((mt_rand() / (double)UINT32_MAX) - 0.5) * 2.0; // center in [0.33, 2.33]
        double width = 0.5 + (mt_rand() / (double)UINT32_MAX); // width in [0.5, 1.5]
        
        // Create an initial 3-point bracket [x1, x2, x3] that contains the minimum.
        // By setting x2 to the random center (close to the actual min) and x1/x3 further
        // away, we ensure f(x2) is likely the lowest, forming a valid initial condition.
        problems[i].x1 = center - width;
        problems[i].x2 = center;
        problems[i].x3 = center + width;
    }
    
    final_result_sum = 0.0;
}

void run_computation() {
    double total_min_x = 0.0;

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        double x1 = problems[i].x1;
        double x2 = problems[i].x2;
        double x3 = problems[i].x3;

        for (int j = 0; j < MAX_STEPS_PER_PROBLEM; ++j) {
            double f1 = objective_function(x1);
            double f2 = objective_function(x2);
            double f3 = objective_function(x3);

            // Calculate the vertex of the parabola passing through the three points.
            // This vertex is the next estimate for the minimum.
            double num = (x2 - x1) * (x2 - x1) * (f2 - f3) - (x2 - x3) * (x2 - x3) * (f2 - f1);
            double den = 2.0 * ((x2 - x1) * (f2 - f3) - (x2 - x3) * (f2 - f1));

            if (fabs(den) < 1e-12) {
                // Denominator is near zero; points are likely collinear. Stop iterating.
                break;
            }

            double x_new = x2 - num / den;

            // Simple update scheme: replace the point with the highest function value.
            // This drives the trio of points towards the minimum.
            if (f1 >= f2 && f1 >= f3) {
                x1 = x_new;
            } else if (f3 >= f1 && f3 >= f2) {
                x3 = x_new;
            } else {
                x2 = x_new;
            }
        }

        // After iterating, find the best minimum from the final three points
        double fx1 = objective_function(x1);
        double fx2 = objective_function(x2);
        double fx3 = objective_function(x3);
        
        double current_min_x;
        if (fx1 <= fx2 && fx1 <= fx3) {
            current_min_x = x1;
        } else if (fx2 <= fx1 && fx2 <= fx3) {
            current_min_x = x2;
        } else {
            current_min_x = x3;
        }
        total_min_x += current_min_x;
    }

    final_result_sum = total_min_x;
}

void cleanup() {
    free(problems);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%f\n", final_result_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
