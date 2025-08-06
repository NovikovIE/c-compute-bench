#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator - Do Not Modify ---
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
int poly_degree;
double lower_bound;
double upper_bound;
int max_steps;
double *poly_coeffs = NULL;
double final_result;

// --- Function to be integrated ---
double polynomial_func(double x) {
    double result = 0.0;
    // Horner's method for efficient evaluation
    for (int i = poly_degree; i >= 0; --i) {
        result = result * x + poly_coeffs[i];
    }
    return result;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s poly_degree lower_bound upper_bound max_steps seed\n", argv[0]);
        exit(1);
    }

    poly_degree = atoi(argv[1]);
    lower_bound = atof(argv[2]);
    upper_bound = atof(argv[3]);
    max_steps = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if (poly_degree <= 0 || max_steps <= 1 || max_steps > 30) {
        fprintf(stderr, "Invalid parameters: poly_degree > 0, 1 < max_steps <= 30\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for polynomial coefficients
    poly_coeffs = (double*)malloc((poly_degree + 1) * sizeof(double));
    if (!poly_coeffs) {
        fprintf(stderr, "Memory allocation failed for coefficients.\n");
        exit(1);
    }

    // Generate random coefficients between -1.0 and 1.0
    for (int i = 0; i <= poly_degree; i++) {
        poly_coeffs[i] = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }
}

void run_computation() {
    // Allocate Romberg table on the heap to avoid stack overflow
    double** R = (double**)malloc(max_steps * sizeof(double*));
    for (int i = 0; i < max_steps; i++) {
        R[i] = (double*)malloc(max_steps * sizeof(double));
    }

    double h = upper_bound - lower_bound;
    R[0][0] = (polynomial_func(lower_bound) + polynomial_func(upper_bound)) * h / 2.0;

    for (int i = 1; i < max_steps; i++) {
        h /= 2.0;
        double sum = 0.0;
        long n_points = 1L << (i - 1); // 2^(i-1)
        for (long k = 1; k <= n_points; k++) {
            sum += polynomial_func(lower_bound + (2 * k - 1) * h);
        }
        // More accurate trapezoidal rule for R[i][0]
        R[i][0] = 0.5 * R[i-1][0] + sum * h; 

        // Richardson extrapolation for R[i][j]
        for (int j = 1; j <= i; j++) {
            double p4 = pow(4.0, j);
            R[i][j] = (p4 * R[i][j-1] - R[i-1][j-1]) / (p4 - 1.0);
        }
    }

    final_result = R[max_steps - 1][max_steps - 1];

    // Free the Romberg table
    for (int i = 0; i < max_steps; i++) {
        free(R[i]);
    }
    free(R);
}

void cleanup() {
    if (poly_coeffs != NULL) {
        free(poly_coeffs);
        poly_coeffs = NULL;
    }
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

    // Print the final result to stdout
    printf("%f\n", final_result);

    // Print the timing info to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
