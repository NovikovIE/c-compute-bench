#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
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
// --- END MERSENNE TWISTER ---

#define NUM_COEFFS 6

// Global struct to hold all benchmark data
typedef struct {
    double (*f)(double, double);
    double initial_y;
    double start_x;
    double end_x;
    double step_size;
    long long num_steps;
    double* coeffs; // For randomized functions
    double accumulated_result;
} BenchmarkData;

static BenchmarkData g_data;

// ODE function: dy/dx = f(x, y)
// Polynomial form with random coefficients
double func_poly(double x, double y) {
    return g_data.coeffs[0] + g_data.coeffs[1] * x + g_data.coeffs[2] * y +
           g_data.coeffs[3] * x * y + g_data.coeffs[4] * x * x + g_data.coeffs[5] * y * y;
}

// ODE function: dy/dx = f(x, y)
// Trigonometric form with random coefficients
double func_trig(double x, double y) {
    return g_data.coeffs[0] * sin(g_data.coeffs[1] * x) + 
           g_data.coeffs[2] * cos(g_data.coeffs[3] * y);
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s function_str initial_y start_x end_x step_size seed\n", argv[0]);
        exit(1);
    }

    const char* function_str = argv[1];
    g_data.initial_y = atof(argv[2]);
    g_data.start_x = atof(argv[3]);
    g_data.end_x = atof(argv[4]);
    g_data.step_size = atof(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);
    
    mt_seed(seed);

    if (g_data.step_size <= 0 || g_data.end_x <= g_data.start_x) {
        fprintf(stderr, "Invalid range or step size.\n");
        exit(1);
    }
    
    g_data.num_steps = (long long)((g_data.end_x - g_data.start_x) / g_data.step_size);
    g_data.accumulated_result = 0.0;

    // Allocate and generate random coefficients for the ODE function
    g_data.coeffs = (double*)malloc(NUM_COEFFS * sizeof(double));
    if (g_data.coeffs == NULL) {
        fprintf(stderr, "Memory allocation failed for coefficients.\n");
        exit(1);
    }
    for (int i = 0; i < NUM_COEFFS; i++) {
        // Generate coefficients in a small range, e.g., [-0.1, 0.1]
        g_data.coeffs[i] = (mt_rand() / (double)UINT32_MAX) * 0.2 - 0.1;
    }

    // Select the ODE function based on the input string
    if (strcmp(function_str, "poly") == 0) {
        g_data.f = &func_poly;
    } else if (strcmp(function_str, "trig") == 0) {
        g_data.f = &func_trig;
    } else {
        fprintf(stderr, "Unknown function string: %s. Use 'poly' or 'trig'.\n", function_str);
        free(g_data.coeffs);
        exit(1);
    }
}

void run_computation() {
    double y = g_data.initial_y;
    double h = g_data.step_size;
    double h_half = h * 0.5;
    double accumulator = 0.0;

    for (long long i = 0; i < g_data.num_steps; ++i) {
        double x = g_data.start_x + i * h;

        // 4th Order Runge-Kutta method
        double k1 = g_data.f(x, y);
        double k2 = g_data.f(x + h_half, y + h_half * k1);
        double k3 = g_data.f(x + h_half, y + h_half * k2);
        double k4 = g_data.f(x + h, y + h * k3);

        y += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

        // To prevent dead-code elimination and produce a result, we accumulate values.
        // Accumulating every step is slow, so we do it periodically.
        if ((i & 1023) == 0) { // i % 1024 == 0
            accumulator += y;
        }
    }
    g_data.accumulated_result = accumulator;
}

void cleanup() {
    free(g_data.coeffs);
    g_data.coeffs = NULL;
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%f\n", g_data.accumulated_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
