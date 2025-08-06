/**
 * @file fea_truss_analysis.c
 * @brief A benchmark simulating Finite Element Analysis (FEA) of a 2D truss structure.
 *
 * This program generates a random 2D truss structure with a specified number of nodes
 * and elements. The core computation involves:
 * 1. Assembling the global stiffness matrix (K) from the properties of each element.
 * 2. Applying boundary conditions to constrain the structure (fixing nodes 0 and 1).
 * 3. Solving the system of linear equations [K]{U} = {F} for the nodal displacements {U},
 *    where {F} is the vector of applied external forces.
 * 4. The system is solved using the Conjugate Gradient (CG) iterative method, which is
 *    well-suited for the symmetric positive-definite stiffness matrices found in FEA.
 *
 * The benchmark's computational load is dominated by the matrix-vector multiplications
 * within the CG solver, which has a complexity of O(iterations * (num_nodes*2)^2) for the
 * dense matrix representation used here.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

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

// --- Benchmark Data Structures ---

int N_NODES;
int N_ELEMENTS;
int N_DOFS; // Degrees of Freedom (2 * N_NODES)

typedef struct {
    double x, y;
} Node;

typedef struct {
    int n1, n2;
    double E, A; // Young's Modulus, Area
} Element;

// Global pointers for benchmark data
Node* nodes;
Element* elements;
double* force_vector;
double* displacement_vector;
double* global_stiffness_matrix;

double final_displacement_sum;

// --- Function Prototypes ---

void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

static double dot_product(int n, const double* v1, const double* v2);
static void mat_vec_mult(int n, const double* mat, const double* vec, double* result);

// --- Benchmark Implementation ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_nodes> <num_elements> <seed>\n", argv[0]);
        exit(1);
    }
    N_NODES = atoi(argv[1]);
    N_ELEMENTS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (N_NODES <= 2 || N_ELEMENTS <= 0) {
        fprintf(stderr, "Number of nodes must be > 2 and elements > 0\n");
        exit(1);
    }

    N_DOFS = N_NODES * 2;
    mt_seed(seed);

    // Memory allocation
    nodes = (Node*)malloc(N_NODES * sizeof(Node));
    elements = (Element*)malloc(N_ELEMENTS * sizeof(Element));
    force_vector = (double*)malloc(N_DOFS * sizeof(double));
    displacement_vector = (double*)malloc(N_DOFS * sizeof(double));
    global_stiffness_matrix = (double*)malloc((size_t)N_DOFS * N_DOFS * sizeof(double));

    if (!nodes || !elements || !force_vector || !displacement_vector || !global_stiffness_matrix) {
        fprintf(stderr, "Failed to allocate memory\n");
        exit(1);
    }

    // Generate nodes randomly in a 100x100 box
    for (int i = 0; i < N_NODES; i++) {
        nodes[i].x = (double)(mt_rand() % 10000) / 100.0;
        nodes[i].y = (double)(mt_rand() % 10000) / 100.0;
    }

    // Generate elements connecting random nodes
    for (int i = 0; i < N_ELEMENTS; i++) {
        int n1, n2;
        do {
            n1 = mt_rand() % N_NODES;
            n2 = mt_rand() % N_NODES;
        } while (n1 == n2);
        elements[i].n1 = n1;
        elements[i].n2 = n2;
        // Realistic-ish properties for steel
        elements[i].E = 200e9; // Pa
        elements[i].A = (double)(mt_rand() % 100 + 1) / 10000.0; // m^2 (1-100 cm^2)
    }

    // Initialize vectors
    memset(displacement_vector, 0, N_DOFS * sizeof(double));
    memset(force_vector, 0, N_DOFS * sizeof(double));

    // Apply some random forces to non-constrained nodes (nodes > 1)
    for (int i = 2; i < N_NODES; ++i) {
        if ((mt_rand() % 10) < 2) { // apply force to ~20% of nodes
            force_vector[i * 2] = (double)(mt_rand() % 20000) - 10000.0;     // Force in X (N)
            force_vector[i * 2 + 1] = (double)(mt_rand() % 20000) - 10000.0; // Force in Y (N)
        }
    }
}

void run_computation() {
    // 1. Assemble Global Stiffness Matrix K
    memset(global_stiffness_matrix, 0, (size_t)N_DOFS * N_DOFS * sizeof(double));

    for (int i = 0; i < N_ELEMENTS; i++) {
        Element el = elements[i];
        Node n1 = nodes[el.n1];
        Node n2 = nodes[el.n2];

        double dx = n2.x - n1.x;
        double dy = n2.y - n1.y;
        double L = sqrt(dx * dx + dy * dy);

        if (L < 1e-9) continue;

        double c = dx / L;
        double s = dy / L;
        double EAL = (el.E * el.A) / L;

        double k_local[4][4] = {
            { c * c,  c * s, -c * c, -c * s },
            { c * s,  s * s, -c * s, -s * s },
            {-c * c, -c * s,  c * c,  c * s },
            {-c * s, -s * s,  c * s,  s * s }
        };

        int dof_indices[4] = { el.n1 * 2, el.n1 * 2 + 1, el.n2 * 2, el.n2 * 2 + 1 };

        for (int r = 0; r < 4; r++) {
            for (int c_idx = 0; c_idx < 4; c_idx++) {
                int global_row = dof_indices[r];
                int global_col = dof_indices[c_idx];
                global_stiffness_matrix[global_row * N_DOFS + global_col] += EAL * k_local[r][c_idx];
            }
        }
    }

    // 2. Apply Boundary Conditions (pin-joint at node 0, roller at node 1 in x-dir)
    // This makes the system statically determinate.
    // Fix node 0 completely (DOFs 0, 1) and node 1 in y-direction (DOF 3)
    int constrained_dofs[] = {0, 1, 3}; // Node 0: x,y ; Node 1: y
    for (int i = 0; i < 3; ++i) {
        int dof = constrained_dofs[i];
        for (int j = 0; j < N_DOFS; ++j) {
            global_stiffness_matrix[dof * N_DOFS + j] = 0.0;
            global_stiffness_matrix[j * N_DOFS + dof] = 0.0;
        }
        global_stiffness_matrix[dof * N_DOFS + dof] = 1.0;
        force_vector[dof] = 0.0;
    }

    // 3. Solve [K]{U} = {F} using Conjugate Gradient
    int max_iter = 100; // Fixed iterations for consistent benchmarking
    double tol = 1e-9;

    double* r = (double*)malloc(N_DOFS * sizeof(double));
    double* p = (double*)malloc(N_DOFS * sizeof(double));
    double* Ap = (double*)malloc(N_DOFS * sizeof(double));

    memcpy(r, force_vector, N_DOFS * sizeof(double));
    memcpy(p, r, N_DOFS * sizeof(double));

    double rs_old = dot_product(N_DOFS, r, r);
    if (rs_old > tol * tol) {
        for (int iter = 0; iter < max_iter; iter++) {
            mat_vec_mult(N_DOFS, global_stiffness_matrix, p, Ap);
            double p_dot_Ap = dot_product(N_DOFS, p, Ap);
            double alpha = (p_dot_Ap == 0) ? 0 : rs_old / p_dot_Ap;
            
            for(int i = 0; i < N_DOFS; ++i) displacement_vector[i] += alpha * p[i];
            for(int i = 0; i < N_DOFS; ++i) r[i] -= alpha * Ap[i];

            double rs_new = dot_product(N_DOFS, r, r);
            if (rs_new < tol * tol) break;

            double beta = rs_new / rs_old;
            for(int i = 0; i < N_DOFS; ++i) p[i] = r[i] + beta * p[i];
            
            rs_old = rs_new;
        }
    }
    free(r); free(p); free(Ap);

    // 4. Calculate final result (sum of absolute displacements)
    final_displacement_sum = 0.0;
    for (int i = 0; i < N_DOFS; i++) {
        final_displacement_sum += fabs(displacement_vector[i]);
    }
}

void cleanup() {
    free(nodes);
    free(elements);
    free(force_vector);
    free(displacement_vector);
    free(global_stiffness_matrix);
}

double dot_product(int n, const double* v1, const double* v2) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

void mat_vec_mult(int n, const double* mat, const double* vec, double* result) {
    for (int i = 0; i < n; i++) {
        result[i] = 0.0;
        for (int j = 0; j < n; j++) {
            result[i] += mat[i * n + j] * vec[j];
        }
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", final_displacement_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}