/*
 * Benchmark: Nonlinear Conjugate Gradient for Rosenbrock Function
 * 
 * Description: This program uses the Fletcher-Reeves nonlinear conjugate gradient
 * method to find the minimum of the multi-dimensional Rosenbrock function.
 * The Rosenbrock function is a classic non-convex function used for testing
 * the performance of optimization algorithms.
 * 
 * f(x) = sum_{i=0}^{N-2} [100 * (x_{i+1} - x_i^2)^2 + (1 - x_i)^2]
 * 
 * The benchmark measures the time taken to perform a fixed number of iterations
 * of the optimization algorithm.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) PRNG --- VERBATIM ---
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
// --- END MT19937 --- 

// --- Benchmark Globals ---
int NUM_DIMENSIONS;
int NUM_ITERATIONS;
double* x;      // Current solution vector
double* g;      // Gradient vector
double* d;      // Search direction vector
double* x_temp; // Temporary vector for line search
double result_accumulator; // To prevent dead code elimination

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_dimensions> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_DIMENSIONS = atoi(argv[1]);
    NUM_ITERATIONS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_DIMENSIONS <= 1 || NUM_ITERATIONS <= 0) {
        fprintf(stderr, "FATAL: num_dimensions and num_iterations must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    x = (double*)malloc(NUM_DIMENSIONS * sizeof(double));
    g = (double*)malloc(NUM_DIMENSIONS * sizeof(double));
    d = (double*)malloc(NUM_DIMENSIONS * sizeof(double));
    x_temp = (double*)malloc(NUM_DIMENSIONS * sizeof(double));

    if (!x || !g || !d || !x_temp) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize starting point x with random values in [-2, 2]
    for (int i = 0; i < NUM_DIMENSIONS; ++i) {
        x[i] = ((double)mt_rand() / (double)UINT32_MAX) * 4.0 - 2.0;
    }

    result_accumulator = 0.0;
}

void cleanup() {
    free(x);
    free(g);
    free(d);
    free(x_temp);
}

// --- Computation Helpers ---

// Calculates the value of the Rosenbrock function
static double rosenbrock(const double* vec) {
    double f = 0.0;
    for (int i = 0; i < NUM_DIMENSIONS - 1; ++i) {
        double term1 = vec[i+1] - vec[i] * vec[i];
        double term2 = 1.0 - vec[i];
        f += 100.0 * term1 * term1 + term2 * term2;
    }
    return f;
}

// Calculates the gradient of the Rosenbrock function
static void calculate_gradient(double* grad, const double* vec) {
    for(int i = 0; i < NUM_DIMENSIONS; ++i) grad[i] = 0.0;

    for (int i = 0; i < NUM_DIMENSIONS - 1; ++i) {
        double term_sq = vec[i+1] - vec[i] * vec[i];
        grad[i]   += -400.0 * vec[i] * term_sq - 2.0 * (1.0 - vec[i]);
        grad[i+1] += 200.0 * term_sq;
    }
}

// Backtracking line search to find a suitable step size alpha
static double line_search(const double* current_x, const double* current_g, const double* search_d) {
    double alpha = 1.0;            // Initial step size
    const double rho = 0.5;        // Contraction factor
    const double c1 = 1e-4;        // Armijo condition constant

    const double f_x = rosenbrock(current_x);
    
    double g_dot_d = 0.0;
    for (int i = 0; i < NUM_DIMENSIONS; ++i) {
        g_dot_d += current_g[i] * search_d[i];
    }
    
    // Try up to 20 backtracking steps
    for (int i = 0; i < 20; ++i) {
        for(int j=0; j<NUM_DIMENSIONS; ++j) {
            x_temp[j] = current_x[j] + alpha * search_d[j];
        }
        
        const double f_x_new = rosenbrock(x_temp);

        if (f_x_new <= f_x + c1 * alpha * g_dot_d) {
            return alpha;
        }
        alpha *= rho;
    }
    return alpha; // Return last alpha if condition never met
}


// --- Core Computation ---
void run_computation() {
    // Initial state calculation
    calculate_gradient(g, x);
    for(int i = 0; i < NUM_DIMENSIONS; ++i) {
        d[i] = -g[i];
    }

    double g_dot_g_old = 0.0;
    for (int i = 0; i < NUM_DIMENSIONS; ++i) {
        g_dot_g_old += g[i] * g[i];
    }

    // Main optimization loop
    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        // 1. Find step size alpha using line search
        double alpha = line_search(x, g, d);
        
        // 2. Update position: x_{k+1} = x_k + alpha * d_k
        for (int i = 0; i < NUM_DIMENSIONS; ++i) {
            x[i] += alpha * d[i];
        }

        // 3. Compute new gradient: g_{k+1} = grad(f(x_{k+1}))
        calculate_gradient(g, x);
        
        // If gradient is close to zero, we might have converged
        if (g_dot_g_old < 1e-12) {
            break;
        }

        // 4. Compute beta using Fletcher-Reeves formula
        double g_dot_g_new = 0.0;
        for (int i = 0; i < NUM_DIMENSIONS; ++i) {
            g_dot_g_new += g[i] * g[i];
        }
        double beta = g_dot_g_new / g_dot_g_old;

        // 5. Update search direction: d_{k+1} = -g_{k+1} + beta * d_k
        for (int i = 0; i < NUM_DIMENSIONS; ++i) {
            d[i] = -g[i] + beta * d[i];
        }

        // 6. Update old gradient norm for next iteration
        g_dot_g_old = g_dot_g_new;
    }

    // Accumulate a result to prevent dead code elimination
    double sum = 0.0;
    for (int i = 0; i < NUM_DIMENSIONS; ++i) {
        sum += x[i];
    }
    result_accumulator = sum;
}


// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
