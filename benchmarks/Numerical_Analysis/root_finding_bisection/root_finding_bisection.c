#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Global structure to hold benchmark data
typedef struct {
    int poly_degree;
    double lower_bound;
    double upper_bound;
    double tolerance;
    double *coeffs;
    double result_root;
} BenchmarkData;

static BenchmarkData g_data;

// Helper function to evaluate the polynomial P(x) using Horner's method
// P(x) = c_n*x^n + c_{n-1}*x^{n-1} + ... + c_1*x + c_0
static double evaluate_poly(double x) {
    double y = 0.0;
    for (int i = g_data.poly_degree; i >= 0; i--) {
        y = y * x + g_data.coeffs[i];
    }
    return y;
}


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s poly_degree lower_bound upper_bound tolerance seed\n", argv[0]);
        exit(1);
    }

    g_data.poly_degree = atoi(argv[1]);
    g_data.lower_bound = atof(argv[2]);
    g_data.upper_bound = atof(argv[3]);
    g_data.tolerance = atof(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);
    
    if (g_data.poly_degree <= 0) {
        fprintf(stderr, "FATAL: poly_degree must be positive.\n");
        exit(1);
    }
    if (g_data.lower_bound >= g_data.upper_bound) {
        fprintf(stderr, "FATAL: lower_bound must be less than upper_bound.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.coeffs = (double*)malloc((g_data.poly_degree + 1) * sizeof(double));
    if (g_data.coeffs == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for polynomial coefficients.\n");
        exit(1);
    }
    
    // Generate random polynomial coefficients until f(lower) and f(upper) have opposite signs,
    // which is a precondition for the bisection method.
    int max_attempts = 1000;
    int attempt = 0;
    while (attempt < max_attempts) {
        for (int i = 0; i <= g_data.poly_degree; i++) {
            // Generate coefficients in a range [-1.0, 1.0] to avoid large values
            g_data.coeffs[i] = ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        }

        if (evaluate_poly(g_data.lower_bound) * evaluate_poly(g_data.upper_bound) < 0) {
            break; // Found a valid polynomial
        }
        attempt++;
    }

    if (attempt == max_attempts) {
        fprintf(stderr, "FATAL: Could not generate a polynomial with a root in [%f, %f] after %d attempts.\n", 
                g_data.lower_bound, g_data.upper_bound, max_attempts);
        free(g_data.coeffs);
        exit(1);
    }

    g_data.result_root = 0.0; // Initialize result
}

void run_computation() {
    double a = g_data.lower_bound;
    double b = g_data.upper_bound;
    
    double fa = evaluate_poly(a);

    // The bisection method requires f(a) and f(b) to have opposite signs.
    // This is ensured by the setup function.
    if (fa * evaluate_poly(b) >= 0) {
        // This case should not be reached due to setup logic. As a failsafe, execution stops.
        g_data.result_root = (a + b) / 2.0;
        return;
    }

    // Calculate a safe maximum number of iterations to guarantee termination.
    int max_iter = (int)ceil(log2((b - a) / g_data.tolerance)) + 2;

    for (int i = 0; i < max_iter; ++i) {
        if ((b - a) < g_data.tolerance) {
            break;
        }
        
        double c = a + (b - a) / 2.0; // Avoid (a+b)/2 to prevent potential overflow
        double fc = evaluate_poly(c);

        if (fc == 0.0) { // Found an exact root
            a = b = c;
            break;
        }

        if (fa * fc < 0) {
            b = c;
        } else {
            a = c;
            fa = fc;
        }
    }
    
    g_data.result_root = (a + b) / 2.0;
}

void cleanup() {
    free(g_data.coeffs);
    g_data.coeffs = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final root to stdout
    printf("%.12f\n", g_data.result_root);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
