#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator --- Do Not Modify ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
int NUM_POINTS;         // Number of points to evaluate the polynomial at.
int DEGREE;             // Degree of the interpolating polynomial.
int NUM_DATA_POINTS;    // Number of data points for interpolation (degree + 1).

double *x_data;         // x-coordinates of data points
double *y_data;         // y-coordinates of data points
double *eval_points;    // Points where the polynomial is evaluated

double final_result;    // Accumulated result to prevent dead code elimination

// Generate a random double in [-100, 100]
double random_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) * 200.0 - 100.0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_points degree seed\n", argv[0]);
        exit(1);
    }

    NUM_POINTS = atoi(argv[1]);
    DEGREE = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_POINTS <= 0 || DEGREE <= 0) {
        fprintf(stderr, "Error: num_points and degree must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    NUM_DATA_POINTS = DEGREE + 1;

    x_data = (double *)malloc(NUM_DATA_POINTS * sizeof(double));
    y_data = (double *)malloc(NUM_DATA_POINTS * sizeof(double));
    eval_points = (double *)malloc(NUM_POINTS * sizeof(double));

    if (!x_data || !y_data || !eval_points) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate distinct data points (x_i, y_i) for interpolation
    for (int i = 0; i < NUM_DATA_POINTS; i++) {
        // To ensure x_i are distinct and avoid division by zero
        x_data[i] = (double)i + (double)mt_rand() / (double)UINT32_MAX * 0.5;
        y_data[i] = random_double();
    }

    // Generate points at which to evaluate the polynomial
    for (int i = 0; i < NUM_POINTS; i++) {
        eval_points[i] = random_double();
    }

    final_result = 0.0;
}

void run_computation() {
    double total_sum = 0.0;

    // For each point we want to evaluate...
    for (int i = 0; i < NUM_POINTS; i++) {
        double eval_x = eval_points[i];
        double interpolated_y = 0.0;

        // ...calculate the interpolated value using Lagrange polynomials.
        // L(x) = sum_{j=0 to n} y_j * l_j(x)
        for (int j = 0; j < NUM_DATA_POINTS; j++) {
            // Calculate the j-th Lagrange basis polynomial, l_j(x).
            // l_j(x) = product_{k=0 to n, k!=j} (x - x_k) / (x_j - x_k)
            double basis_poly = 1.0;
            for (int k = 0; k < NUM_DATA_POINTS; k++) {
                if (j != k) {
                    basis_poly *= (eval_x - x_data[k]) / (x_data[j] - x_data[k]);
                }
            }
            interpolated_y += y_data[j] * basis_poly;
        }
        total_sum += interpolated_y;
    }

    final_result = total_sum;
}

void cleanup() {
    free(x_data);
    free(y_data);
    free(eval_points);
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

    // Print final result to stdout
    printf("%f\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
