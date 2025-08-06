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
typedef struct {
    int num_points;
    int num_queries;
    double* x_coords;       // x-coordinates of known data points
    double* y_coords;       // y-coordinates of known data points
    double* weights;        // Precomputed barycentric weights
    double* query_points;   // Points at which to interpolate
    double final_result;    // Accumulated result to prevent dead code elimination
} BenchmarkData;

BenchmarkData g_data;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

// --- Benchmark Implementation ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <num_queries> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_points = atoi(argv[1]);
    g_data.num_queries = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_points <= 1 || g_data.num_queries <= 0) {
        fprintf(stderr, "FATAL: num_points must be > 1 and num_queries > 0.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    g_data.x_coords = (double*)malloc(g_data.num_points * sizeof(double));
    g_data.y_coords = (double*)malloc(g_data.num_points * sizeof(double));
    g_data.weights = (double*)malloc(g_data.num_points * sizeof(double));
    g_data.query_points = (double*)malloc(g_data.num_queries * sizeof(double));

    if (!g_data.x_coords || !g_data.y_coords || !g_data.weights || !g_data.query_points) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate distinct, sorted x_coords to simplify generation and avoid duplicate points.
    g_data.x_coords[0] = mt_rand() / (double)UINT32_MAX;
    for (int i = 1; i < g_data.num_points; ++i) {
        g_data.x_coords[i] = g_data.x_coords[i-1] + (mt_rand() / (double)UINT32_MAX) + 1e-5; // Ensure points are distinct
    }

    // Generate corresponding y_coords.
    for (int i = 0; i < g_data.num_points; ++i) {
        g_data.y_coords[i] = (mt_rand() / (double)UINT32_MAX) * 200.0 - 100.0; // Values between -100 and 100
    }

    // Precompute barycentric weights (O(n^2) complexity).
    for (int i = 0; i < g_data.num_points; ++i) {
        g_data.weights[i] = 1.0;
        for (int j = 0; j < g_data.num_points; ++j) {
            if (i == j) continue;
            g_data.weights[i] /= (g_data.x_coords[i] - g_data.x_coords[j]);
        }
    }

    // Generate query points within the range of x_coords.
    double x_min = g_data.x_coords[0];
    double x_max = g_data.x_coords[g_data.num_points - 1];
    for (int i = 0; i < g_data.num_queries; ++i) {
        g_data.query_points[i] = x_min + (mt_rand() / (double)UINT32_MAX) * (x_max - x_min);
    }
    
    g_data.final_result = 0.0;
}

void run_computation() {
    double total_sum = 0.0;
    const double epsilon = 1e-12;

    for (int i = 0; i < g_data.num_queries; ++i) {
        double t = g_data.query_points[i];
        double numerator_sum = 0.0;
        double denominator_sum = 0.0;
        int exact_match_index = -1;

        // Special case: if query point t is very close to a known point x_j, the interpolated value is y_j.
        // This avoids division by zero.
        for (int j = 0; j < g_data.num_points; ++j) {
            if (fabs(t - g_data.x_coords[j]) < epsilon) {
                exact_match_index = j;
                break;
            }
        }

        if (exact_match_index != -1) {
            total_sum += g_data.y_coords[exact_match_index];
        } else {
            // Use the barycentric interpolation formula.
            for (int j = 0; j < g_data.num_points; ++j) {
                double term = g_data.weights[j] / (t - g_data.x_coords[j]);
                numerator_sum += term * g_data.y_coords[j];
                denominator_sum += term;
            }
            if (fabs(denominator_sum) > epsilon) {
                 total_sum += numerator_sum / denominator_sum;
            }
        }
    }
    g_data.final_result = total_sum;
}

void cleanup() {
    free(g_data.x_coords);
    free(g_data.y_coords);
    free(g_data.weights);
    free(g_data.query_points);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout.
    printf("%f\n", g_data.final_result);

    // Print timing information to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
