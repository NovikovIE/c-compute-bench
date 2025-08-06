#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

// --- BENCHMARK DATA AND PARAMETERS ---

// Simplified model for a nonlinear element (e.g., a diode)
typedef struct {
    int node_p;   // Positive terminal node index
    int node_n;   // Negative terminal node index
    double is;    // Saturation current
    double nvt;   // Ideality factor * thermal voltage
} NonlinearElement;

// Global structure to hold all benchmark data
typedef struct {
    int num_nodes;
    int num_nonlinear_elements;
    int max_iterations;
    uint32_t seed;

    double *g_matrix;         // Conductance matrix (flattened)
    double *i_vector;         // Current source vector
    double *x_vector;         // Node voltage vector (the solution)
    NonlinearElement *nonlinear_elements;

    double final_result; // To prevent dead code elimination
} BenchmarkData;

static BenchmarkData g_data;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_nodes> <num_nonlinear_elements> <max_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_nodes = atoi(argv[1]);
    g_data.num_nonlinear_elements = atoi(argv[2]);
    g_data.max_iterations = atoi(argv[3]);
    g_data.seed = (uint32_t)strtoul(argv[4], NULL, 10);

    mt_seed(g_data.seed);

    // Allocate memory
    size_t matrix_size = (size_t)g_data.num_nodes * g_data.num_nodes;
    g_data.g_matrix = (double*)malloc(matrix_size * sizeof(double));
    g_data.i_vector = (double*)malloc(g_data.num_nodes * sizeof(double));
    g_data.x_vector = (double*)malloc(g_data.num_nodes * sizeof(double));
    g_data.nonlinear_elements = (NonlinearElement*)malloc(g_data.num_nonlinear_elements * sizeof(NonlinearElement));

    if (!g_data.g_matrix || !g_data.i_vector || !g_data.x_vector || !g_data.nonlinear_elements) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize G matrix (linear conductances)
    // Create a diagonally dominant matrix to aid solver convergence
    for (int i = 0; i < g_data.num_nodes; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < g_data.num_nodes; ++j) {
            if (i != j) {
                double val = (double)mt_rand() / (double)UINT32_MAX * 1e-4; // Small off-diagonal conductance
                g_data.g_matrix[i * g_data.num_nodes + j] = -val;
                row_sum += val;
            }
        }
        // Diagonal element is the sum of others plus a connection to ground
        g_data.g_matrix[i * g_data.num_nodes + i] = row_sum + (double)mt_rand() / (double)UINT32_MAX * 1e-3;
    }

    // Initialize I vector (current sources) and initial X vector (voltages)
    for (int i = 0; i < g_data.num_nodes; ++i) {
        g_data.i_vector[i] = ((double)mt_rand() / (double)UINT32_MAX - 0.5) * 1e-3; // Small random currents
        g_data.x_vector[i] = 0.0; // Initial guess for voltages
    }

    // Initialize nonlinear elements
    for (int i = 0; i < g_data.num_nonlinear_elements; ++i) {
        int n1 = mt_rand() % g_data.num_nodes;
        int n2;
        do {
            n2 = mt_rand() % g_data.num_nodes;
        } while (n1 == n2);
        
        g_data.nonlinear_elements[i].node_p = n1;
        g_data.nonlinear_elements[i].node_n = n2;
        g_data.nonlinear_elements[i].is = 1e-12 + ((double)mt_rand() / (double)UINT32_MAX) * 1e-13; // Typical diode saturation current
        g_data.nonlinear_elements[i].nvt = 0.026 + ((double)mt_rand() / (double)UINT32_MAX) * 0.01; // n*Vt
    }
}

void run_computation() {
    int N = g_data.num_nodes;
    int M = g_data.num_nonlinear_elements;
    const int JACOBI_ITERATIONS = 10;

    // Temporary storage for Newton-Raphson method
    double *jacobian = (double*)malloc((size_t)N * N * sizeof(double));
    double *residual = (double*)malloc(N * sizeof(double));
    double *delta_x = (double*)malloc(N * sizeof(double));
    double *next_delta_x = (double*)malloc(N * sizeof(double));
    if (!jacobian || !residual || !delta_x || !next_delta_x) {
        fprintf(stderr, "FATAL: Temp memory allocation failed in computation.\n");
        exit(1);
    }
    
    for (int iter = 0; iter < g_data.max_iterations; ++iter) {
        // --- 1. Formulate the system of equations f(x) = J*x - b = 0 ---
        
        // Start with linear part: residual = G*x - I, jacobian = G
        for (int i = 0; i < N; ++i) {
            double gx_sum = 0.0;
            for (int j = 0; j < N; ++j) {
                gx_sum += g_data.g_matrix[i * N + j] * g_data.x_vector[j];
                jacobian[i * N + j] = g_data.g_matrix[i * N + j];
            }
            residual[i] = gx_sum - g_data.i_vector[i];
        }

        // Add contributions from nonlinear elements
        for (int i = 0; i < M; ++i) {
            NonlinearElement *elem = &g_data.nonlinear_elements[i];
            int p = elem->node_p;
            int n = elem->node_n;
            double vp = g_data.x_vector[p];
            double vn = g_data.x_vector[n];
            double vd = vp - vn;

            // Clamp vd to avoid floating point overflow in exp()
            if (vd > 0.7) vd = 0.7;
            if (vd < -0.7) vd = -0.7;
            
            double exp_term = exp(vd / elem->nvt);
            double id = elem->is * (exp_term - 1.0);
            double gd = (elem->is / elem->nvt) * exp_term;

            // Stamp into residual vector
            residual[p] += id;
            residual[n] -= id;

            // Stamp into Jacobian matrix
            jacobian[p * N + p] += gd;
            jacobian[n * N + n] += gd;
            jacobian[p * N + n] -= gd;
            jacobian[n * N + p] -= gd;
        }

        // --- 2. Solve the linear system J * delta_x = -residual ---
        // Using a few iterations of the Jacobi method.
        for(int i = 0; i < N; i++) delta_x[i] = 0.0; // Initial guess for delta_x

        for (int j_iter = 0; j_iter < JACOBI_ITERATIONS; ++j_iter) {
            for (int i = 0; i < N; ++i) {
                double sigma = 0.0;
                for (int j = 0; j < N; ++j) {
                    if (i != j) {
                        sigma += jacobian[i * N + j] * delta_x[j];
                    }
                }
                double diag = jacobian[i * N + i];
                if (fabs(diag) < 1e-12) diag = 1e-12; // Avoid division by zero
                next_delta_x[i] = (1.0 / diag) * (-residual[i] - sigma);
            }
            // Swap pointers for next iteration
            double *temp = delta_x;
            delta_x = next_delta_x;
            next_delta_x = temp;
        }

        // --- 3. Update the solution vector ---
        for (int i = 0; i < N; ++i) {
            g_data.x_vector[i] += delta_x[i];
        }
    }

    // Calculate final result checksum to prevent optimization
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += g_data.x_vector[i];
    }
    g_data.final_result = sum;

    free(jacobian);
    free(residual);
    free(delta_x);
    free(next_delta_x);
}

void cleanup() {
    free(g_data.g_matrix);
    free(g_data.i_vector);
    free(g_data.x_vector);
    free(g_data.nonlinear_elements);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
