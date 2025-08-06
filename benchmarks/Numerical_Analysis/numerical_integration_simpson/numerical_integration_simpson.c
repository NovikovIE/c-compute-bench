#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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
// --- End of MT19937 --- 

// Function pointer type for the integrand
typedef double (*IntegrandFunc)(double);

// Global data structure for the benchmark
typedef struct {
    IntegrandFunc func;
    double lower_bound;
    double upper_bound;
    long long num_intervals;
    double result;
    // Dummy allocation to satisfy requirements
    void* dummy_ptr;
} BenchmarkData;

BenchmarkData g_data;

// ---- Functions to be integrated ----
double poly_func(double x) {
    return x * x * x - 2.0 * x * x + 5.0;
}

double trig_func(double x) {
    return sin(x) * cos(2.0 * x) + sin(x*x);
}

double exp_combo_func(double x) {
    return exp(-x * x / 2.0) * sin(x) * 0.5;
}

// --- Benchmark Functions --- 
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <function_str> <lower_bound> <upper_bound> <num_intervals> <seed>\n", argv[0]);
        exit(1);
    }

    char* function_str = argv[1];
    g_data.lower_bound = atof(argv[2]);
    g_data.upper_bound = atof(argv[3]);
    g_data.num_intervals = atoll(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed); // Seed the generator as required

    if (strcmp(function_str, "polynomial") == 0) {
        g_data.func = poly_func;
    } else if (strcmp(function_str, "trig") == 0) {
        g_data.func = trig_func;
    } else if (strcmp(function_str, "exp_combo") == 0) {
        g_data.func = exp_combo_func;
    } else {
        fprintf(stderr, "Invalid function string specified.\n");
        exit(1);
    }

    if (g_data.num_intervals <= 0) {
        fprintf(stderr, "Number of intervals must be positive.\n");
        exit(1);
    }
    
    // Simpson's rule requires an even number of intervals
    if (g_data.num_intervals % 2 != 0) {
        g_data.num_intervals++;
    }

    // Satisfy heap allocation requirement
    g_data.dummy_ptr = malloc(128);
    if (g_data.dummy_ptr == NULL) {
        fprintf(stderr, "Failed to allocate memory.\n");
        exit(1);
    }

    g_data.result = 0.0;
}

void run_computation() {
    double h = (g_data.upper_bound - g_data.lower_bound) / (double)g_data.num_intervals;
    double sum = g_data.func(g_data.lower_bound) + g_data.func(g_data.upper_bound);

    // Summation of terms with coefficient 4
    double sum4 = 0.0;
    for (long long i = 1; i < g_data.num_intervals; i += 2) {
        double x = g_data.lower_bound + i * h;
        sum4 += g_data.func(x);
    }

    // Summation of terms with coefficient 2
    double sum2 = 0.0;
    for (long long i = 2; i < g_data.num_intervals; i += 2) {
        double x = g_data.lower_bound + i * h;
        sum2 += g_data.func(x);
    }

    sum += 4.0 * sum4 + 2.0 * sum2;

    g_data.result = (h / 3.0) * sum;
}

void cleanup() {
    if (g_data.dummy_ptr) {
        free(g_data.dummy_ptr);
        g_data.dummy_ptr = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%.12f\n", g_data.result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
