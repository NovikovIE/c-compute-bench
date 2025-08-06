#include <stdio.h>
#include <stdlib.h>
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

// --- BENCHMARK DATA AND PARAMETERS ---
int num_points;
int num_eval_points;
double *x_pts;       // x-coordinates of data points
double *y_pts;       // y-coordinates of data points
double *eval_x;     // x-coordinates for evaluation
double *eval_y;     // Storage for interpolated results
double final_result; // Accumulated result

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <num_eval_points> <seed>\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    num_eval_points = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    x_pts = (double *)malloc(num_points * sizeof(double));
    y_pts = (double *)malloc(num_points * sizeof(double));
    eval_x = (double *)malloc(num_eval_points * sizeof(double));
    eval_y = (double *)malloc(num_eval_points * sizeof(double));

    if (!x_pts || !y_pts || !eval_x || !eval_y) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate sorted data points with uniform spacing in x
    for (int i = 0; i < num_points; ++i) {
        x_pts[i] = (double)i;
        y_pts[i] = ((double)mt_rand() / (double)UINT32_MAX) * 100.0;
    }

    // Generate random points to evaluate the interpolation
    double x_max = (double)(num_points - 1);
    for (int i = 0; i < num_eval_points; ++i) {
        eval_x[i] = ((double)mt_rand() / (double)UINT32_MAX) * x_max;
    }
    
    final_result = 0.0;
}

void run_computation() {
    // Akima spline requires calculating slopes and then tangents at each point.
    // We allocate temporary arrays for these intermediate values.
    int n = num_points;
    if (n < 5) {
        // Akima spline is not well-defined for fewer than 5 points due to the 'window' needed.
        // A simpler interpolation or error could be handled here.
        // For this benchmark, we assume n is large enough.
        return;
    }

    // m_padded has space for extrapolated slopes: m[-2], m[-1], m[0]..m[n-2], m[n-1], m[n]
    // Total size: (n-1) real slopes + 4 extrapolated = n+3
    double* m_padded = (double*)malloc((n + 3) * sizeof(double));
    double* m = m_padded + 2; // Pointer to the logical start of the slope array (m[0])

    double* t = (double*)malloc(n * sizeof(double)); // Tangents at each point

    // 1. Calculate slopes between points
    for (int i = 0; i < n - 1; ++i) {
        m[i] = y_pts[i+1] - y_pts[i]; // Assumes x[i+1]-x[i] == 1
    }

    // 2. Extrapolate slopes at boundaries
    m[-1] = 2.0 * m[0] - m[1];
    m[-2] = 2.0 * m[-1] - m[0];
    m[n-1] = 2.0 * m[n-2] - m[n-3];
    m[n]   = 2.0 * m[n-1] - m[n-2];

    // 3. Calculate Akima tangents
    for (int i = 0; i < n; ++i) {
        double w1 = fabs(m[i+1] - m[i]);
        double w2 = fabs(m[i-1] - m[i-2]);
        if ((w1 + w2) > 1e-9) { // Use a small epsilon to avoid division by zero
            t[i] = (w1 * m[i-1] + w2 * m[i]) / (w1 + w2);
        } else {
            t[i] = 0.5 * (m[i-1] + m[i]);
        }
    }

    // 4. Interpolate values for each evaluation point
    for (int j = 0; j < num_eval_points; ++j) {
        double p = eval_x[j];
        int i = (int)floor(p);

        // Clamp index to a valid interval [0, n-2]
        if (i < 0) i = 0;
        if (i >= n - 1) i = n - 2;

        // Use cubic Hermite interpolation basis functions
        double h = 1.0; // x_pts[i+1] - x_pts[i]
        double t_norm = (p - x_pts[i]) / h;

        double t2 = t_norm * t_norm;
        double t3 = t_norm * t2;

        double h00 = 2*t3 - 3*t2 + 1;
        double h10 = t3 - 2*t2 + t_norm;
        double h01 = -2*t3 + 3*t2;
        double h11 = t3 - t2;

        eval_y[j] = h00*y_pts[i] + h10*h*t[i] + h01*y_pts[i+1] + h11*h*t[i+1];
        final_result += eval_y[j];
    }

    free(m_padded);
    free(t);
}

void cleanup() {
    free(x_pts);
    free(y_pts);
    free(eval_x);
    free(eval_y);
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
    printf("%f\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
