#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Configuration ---
typedef struct {
    double p1;          // Parameter 1 for the objective function
    double p2;          // Parameter 2 for the objective function
    double lower_bound; // 'a' bracket for the minimum
    double upper_bound; // 'b' bracket for the minimum
} ProblemInstance;

// Global variables for benchmark data and parameters
static int NUM_ITERATIONS;
static double PRECISION;
static ProblemInstance *problems;
static double final_result_accumulator; 

// --- Helper Functions ---
double rand_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

double objective_function(double x, double p1, double p2) {
    double term1 = x - p1;
    // A function with a global minimum and several local minima
    return term1 * term1 + p2 * cos(10.0 * x);
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_iterations> <precision> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_ITERATIONS = atoi(argv[1]);
    PRECISION = atof(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    problems = (ProblemInstance *)malloc(NUM_ITERATIONS * sizeof(ProblemInstance));
    if (problems == NULL) {
        fprintf(stderr, "Failed to allocate memory for problems\n");
        exit(1);
    }

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        problems[i].p1 = rand_double(-10.0, 10.0);
        problems[i].p2 = rand_double(0.1, 2.0);
        problems[i].lower_bound = rand_double(-20.0, 0.0);
        problems[i].upper_bound = rand_double(0.0, 20.0);
    }
}

void run_computation() {
    final_result_accumulator = 0.0;

    const double C = 0.5 * (3.0 - sqrt(5.0)); // Golden ratio
    const int MAX_BRENT_ITERS = 100; // Safeguard max iterations per problem
    const double ZEPS = 1.0e-11; // Small number to prevent division by zero

    for (int i = 0; i < NUM_ITERATIONS; ++i) {
        double a = problems[i].lower_bound;
        double b = problems[i].upper_bound;
        double p1 = problems[i].p1;
        double p2 = problems[i].p2;

        double d = 0.0, e = 0.0, u, v, w, x;
        double fu, fv, fw, fx;

        x = w = v = a + C * (b - a);
        fx = fw = fv = objective_function(x, p1, p2);

        for (int iter = 0; iter < MAX_BRENT_ITERS; ++iter) {
            double m = 0.5 * (a + b);
            double tol1 = PRECISION * fabs(x) + ZEPS;
            double tol2 = 2.0 * tol1;

            if (fabs(x - m) <= (tol2 - 0.5 * (b - a))) {
                break; // Met precision
            }

            if (fabs(e) > tol1) { // Try parabolic fit
                double r = (x - w) * (fx - fv);
                double q = (x - v) * (fx - fw);
                double p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                if (q > 0.0) p = -p;
                q = fabs(q);
                double etemp = e;
                e = d;
                if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                    d = C * (e = (x >= m ? a - x : b - x)); // Fallback to golden section
                } else {
                    d = p / q; // Parabolic step
                    u = x + d;
                    if (u - a < tol2 || b - u < tol2) {
                        d = (x < m) ? tol1 : -tol1;
                    }
                }
            } else { // Golden section step
                d = C * (e = (x >= m ? a - x : b - x));
            }
            u = x + ((fabs(d) >= tol1) ? d : (d > 0.0 ? tol1 : -tol1));
            fu = objective_function(u, p1, p2);

            if (fu <= fx) {
                if (u >= x) a = x; else b = x;
                v = w; w = x; x = u;
                fv = fw; fw = fx; fx = fu;
            } else {
                if (u < x) a = u; else b = u;
                if (fu <= fw || w == x) {
                    v = w; w = u;
                    fv = fw; fw = fu;
                } else if (fu <= fv || v == x || v == w) {
                    v = u; fv = fu;
                }
            }
        }
        final_result_accumulator += x;
    }
}

void cleanup() {
    free(problems);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%f\n", final_result_accumulator);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
