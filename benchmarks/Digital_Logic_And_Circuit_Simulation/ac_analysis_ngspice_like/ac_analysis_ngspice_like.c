#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- BENCHMARK SPECIFIC CODE ---

#define JACOBI_ITERATIONS 5
#define OUTPUT_NODE 0

// Complex number representation
typedef struct {
    double real;
    double imag;
} Complex;

// Global benchmark parameters and data structures
int num_nodes;
int num_active_elements;
int num_frequency_points;

// MNA matrices (G: conductance, C: susceptance components)
double** G_matrix;
double** C_matrix;

// Source vector and frequencies
Complex* I_vector;
double* frequencies;

double final_result; // Accumulated result to prevent dead code elimination

// --- COMPLEX NUMBER ARITHMETIC ---

Complex c_add(Complex a, Complex b) {
    return (Complex){a.real + b.real, a.imag + b.imag};
}

Complex c_sub(Complex a, Complex b) {
    return (Complex){a.real - b.real, a.imag - b.imag};
}

Complex c_mul(Complex a, Complex b) {
    return (Complex){
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real
    };
}

Complex c_div(Complex a, Complex b) {
    double denom = b.real * b.real + b.imag * b.imag;
    return (Complex){
        (a.real * b.real + a.imag * b.imag) / denom,
        (a.imag * b.real - a.real * b.imag) / denom
    };
}

double c_abs(Complex a) {
    return sqrt(a.real * a.real + a.imag * a.imag);
}

// Generates a random double between 0 and 1
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_nodes num_active_elements num_frequency_points seed\n", argv[0]);
        exit(1);
    }

    num_nodes = atoi(argv[1]);
    num_active_elements = atoi(argv[2]);
    num_frequency_points = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    mt_seed(seed);

    // Allocate matrices and vectors
    G_matrix = (double**)malloc(num_nodes * sizeof(double*));
    C_matrix = (double**)malloc(num_nodes * sizeof(double*));
    for (int i = 0; i < num_nodes; i++) {
        G_matrix[i] = (double*)calloc(num_nodes, sizeof(double));
        C_matrix[i] = (double*)calloc(num_nodes, sizeof(double));
    }

    I_vector = (Complex*)calloc(num_nodes, sizeof(Complex));
    frequencies = (double*)malloc(num_frequency_points * sizeof(double));

    // Populate G and C matrices with random elements to simulate a circuit netlist
    // This simulates stamping resistors, capacitors, and VCCSs
    for (int i = 0; i < num_active_elements; i++) {
        int n1 = mt_rand() % num_nodes;
        int n2 = mt_rand() % num_nodes;
        int n3 = mt_rand() % num_nodes;
        int n4 = mt_rand() % num_nodes;
        double r_val = 1.0 / (1.0 + rand_double() * 10000.0); // conductance
        double c_val = rand_double() * 1e-9; // capacitance
        double gm_val = rand_double() * 0.1;   // transconductance

        // Stamp resistor
        if (n1 != n2) {
            G_matrix[n1][n1] += r_val;
            G_matrix[n2][n2] += r_val;
            G_matrix[n1][n2] -= r_val;
            G_matrix[n2][n1] -= r_val;
        }

        // Stamp capacitor
        if (n1 != n2) {
            C_matrix[n1][n1] += c_val;
            C_matrix[n2][n2] += c_val;
            C_matrix[n1][n2] -= c_val;
            C_matrix[n2][n1] -= c_val;
        }
        
        // Stamp VCCS
        if (n1 != n2 && n3 != n4) {
            G_matrix[n1][n3] += gm_val;
            G_matrix[n2][n4] += gm_val;
            G_matrix[n1][n4] -= gm_val;
            G_matrix[n2][n3] -= gm_val;
        }
    }

    // Ensure diagonal is non-zero for solver stability
    for (int i = 0; i < num_nodes; i++) {
        G_matrix[i][i] += 1e-9;
    }

    // Set up a single current source
    I_vector[0] = (Complex){1.0, 0.0};

    // Generate logarithmically spaced frequencies (e.g., 1kHz to 1GHz)
    double freq_start = 1e3;
    double freq_stop = 1e9;
    double log_ratio = pow(freq_stop / freq_start, 1.0 / (num_frequency_points - 1));
    for (int i = 0; i < num_frequency_points; i++) {
        frequencies[i] = freq_start * pow(log_ratio, i);
    }
}

void run_computation() {
    final_result = 0.0;

    // Allocate work arrays for computation
    Complex* V_solution = (Complex*)calloc(num_nodes, sizeof(Complex));
    Complex* V_next = (Complex*)malloc(num_nodes * sizeof(Complex));
    Complex** Y_matrix = (Complex**)malloc(num_nodes * sizeof(Complex*));
    for(int i = 0; i < num_nodes; ++i) {
        Y_matrix[i] = (Complex*)malloc(num_nodes * sizeof(Complex));
    }

    // Main loop: iterate over each frequency point
    for (int f = 0; f < num_frequency_points; f++) {
        double omega = 2.0 * M_PI * frequencies[f];

        // 1. Build the complex admittance matrix Y = G + jw*C
        for (int i = 0; i < num_nodes; i++) {
            for (int j = 0; j < num_nodes; j++) {
                Y_matrix[i][j] = (Complex){G_matrix[i][j], omega * C_matrix[i][j]};
            }
        }

        // 2. Solve YV=I using Jacobi iterative method
        for (int i = 0; i < num_nodes; i++) {
            V_solution[i] = (Complex){0.0, 0.0};
        }

        for (int iter = 0; iter < JACOBI_ITERATIONS; iter++) {
            for (int i = 0; i < num_nodes; i++) {
                Complex sigma = {0.0, 0.0};
                for (int j = 0; j < num_nodes; j++) {
                    if (i != j) {
                        sigma = c_add(sigma, c_mul(Y_matrix[i][j], V_solution[j]));
                    }
                }
                V_next[i] = c_div(c_sub(I_vector[i], sigma), Y_matrix[i][i]);
            }
            // Copy V_next to V_solution for next iteration
            for(int i = 0; i < num_nodes; ++i) {
                V_solution[i] = V_next[i];
            }
        }

        // 3. Accumulate a result from the solution vector
        final_result += c_abs(V_solution[OUTPUT_NODE]);
    }
    
    // Free work arrays
    free(V_solution);
    free(V_next);
    for(int i=0; i<num_nodes; ++i) {
        free(Y_matrix[i]);
    }
    free(Y_matrix);
}

void cleanup() {
    for (int i = 0; i < num_nodes; i++) {
        free(G_matrix[i]);
        free(C_matrix[i]);
    }
    free(G_matrix);
    free(C_matrix);
    free(I_vector);
    free(frequencies);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
