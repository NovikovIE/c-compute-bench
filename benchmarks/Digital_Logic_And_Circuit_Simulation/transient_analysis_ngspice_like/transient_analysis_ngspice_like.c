#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <math.h>

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

// Represents a simple non-linear active element (e.g., a Voltage-Controlled Current Source)
typedef struct {
    int node_out;
    int control_node1;
    int control_node2;
    double transconductance_factor;
} ActiveElement;

// Global structure to hold all benchmark data
struct {
    int num_nodes;
    int num_active_elements;
    double simulation_time_end;
    double time_step;

    double **G_static;     // Static part of the Modified Nodal Analysis matrix (resistors)
    double **G_transient;  // Full MNA matrix for the current time step
    double *x;             // Vector of node voltages (the unknowns)
    double *b;             // Vector of sources (currents)
    ActiveElement *active_elements;

    double final_result;
} g_data;

double** allocate_matrix(int rows, int cols) {
    double** matrix = (double**)malloc(rows * sizeof(double*));
    if (!matrix) return NULL;
    for (int i = 0; i < rows; ++i) {
        matrix[i] = (double*)malloc(cols * sizeof(double));
        if (!matrix[i]) {
            // Free previously allocated rows
            for (int k = 0; k < i; ++k) free(matrix[k]);
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_nodes> <num_active_elements> <simulation_time_end> <time_step> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_nodes = atoi(argv[1]);
    g_data.num_active_elements = atoi(argv[2]);
    g_data.simulation_time_end = atof(argv[3]);
    g_data.time_step = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);
    mt_seed(seed);

    // Allocate matrices and vectors
    g_data.G_static = allocate_matrix(g_data.num_nodes, g_data.num_nodes);
    g_data.G_transient = allocate_matrix(g_data.num_nodes, g_data.num_nodes);
    g_data.x = (double*)malloc(g_data.num_nodes * sizeof(double));
    g_data.b = (double*)malloc(g_data.num_nodes * sizeof(double));
    g_data.active_elements = (ActiveElement*)malloc(g_data.num_active_elements * sizeof(ActiveElement));

    if (!g_data.G_static || !g_data.G_transient || !g_data.x || !g_data.b || !g_data.active_elements) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize G_static with random conductances, ensuring diagonal dominance for stability
    for (int i = 0; i < g_data.num_nodes; ++i) {
        double diag_sum = 0.0;
        for (int j = 0; j < g_data.num_nodes; ++j) {
            if (i != j) {
                double conductance = (double)mt_rand() / UINT32_MAX * 0.1;
                g_data.G_static[i][j] = -conductance;
                diag_sum += conductance;
            }
        }
        g_data.G_static[i][i] = diag_sum + (double)mt_rand() / UINT32_MAX * 0.5 + 1.0; // Add to diagonal
    }

    // Initialize vectors x and b
    memset(g_data.x, 0, g_data.num_nodes * sizeof(double));
    for (int i = 0; i < g_data.num_nodes; ++i) {
        g_data.b[i] = ((double)mt_rand() / UINT32_MAX - 0.5) * 0.1; // Small random initial currents
    }

    // Initialize active elements with random connections
    for (int i = 0; i < g_data.num_active_elements; ++i) {
        g_data.active_elements[i].node_out = mt_rand() % g_data.num_nodes;
        g_data.active_elements[i].control_node1 = mt_rand() % g_data.num_nodes;
        g_data.active_elements[i].control_node2 = mt_rand() % g_data.num_nodes;
        g_data.active_elements[i].transconductance_factor = ((double)mt_rand() / UINT32_MAX) * 0.05;
    }
    
    g_data.final_result = 0.0;
}

void run_computation() {
    const int SOLVER_ITERATIONS = 5; // Fixed number of iterations for the linear solver
    double accumulator = 0.0;
    double current_time = 0.0;

    while (current_time < g_data.simulation_time_end) {
        // 1. Copy static matrix to transient matrix for this time step
        for (int i = 0; i < g_data.num_nodes; ++i) {
            memcpy(g_data.G_transient[i], g_data.G_static[i], g_data.num_nodes * sizeof(double));
        }

        // 2. Update source vector b (e.g., a sinusoidal source at node 0)
        g_data.b[0] = 5.0 * sin(2.0 * M_PI * 50.0 * current_time);

        // 3. Add non-linear element contributions to the transient G matrix
        // This models elements whose behavior depends on the state (voltages) from the previous step
        for (int i = 0; i < g_data.num_active_elements; ++i) {
            ActiveElement* elem = &g_data.active_elements[i];
            double V_control = g_data.x[elem->control_node1] - g_data.x[elem->control_node2];
            double gm = elem->transconductance_factor * tanh(V_control); // Use tanh for a soft non-linearity
            
            g_data.G_transient[elem->node_out][elem->control_node1] += gm;
            g_data.G_transient[elem->node_out][elem->control_node2] -= gm;
        }

        // 4. Solve the linear system G*x = b using Gauss-Seidel iterative method
        for (int iter = 0; iter < SOLVER_ITERATIONS; ++iter) {
            for (int i = 0; i < g_data.num_nodes; ++i) {
                double sigma = 0.0;
                for (int j = 0; j < g_data.num_nodes; ++j) {
                    if (i != j) {
                        sigma += g_data.G_transient[i][j] * g_data.x[j];
                    }
                }
                // Avoid division by zero, though static matrix is built to be non-singular
                if (fabs(g_data.G_transient[i][i]) > 1e-12) {
                    double new_xi = (g_data.b[i] - sigma) / g_data.G_transient[i][i];
                    // Apply damping to help convergence
                    g_data.x[i] = g_data.x[i] * 0.5 + new_xi * 0.5;
                }
            }
        }

        // 5. Accumulate a result to prevent dead code elimination
        for (int i = 0; i < g_data.num_nodes; ++i) {
            accumulator += g_data.x[i];
        }

        current_time += g_data.time_step;
    }

    g_data.final_result = accumulator;
}

void cleanup() {
    for (int i = 0; i < g_data.num_nodes; ++i) {
        free(g_data.G_static[i]);
        free(g_data.G_transient[i]);
    }
    free(g_data.G_static);
    free(g_data.G_transient);
    free(g_data.x);
    free(g_data.b);
    free(g_data.active_elements);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
