#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Benchmark Globals ---
static int num_points;
static int num_queries;
static double *x_values;     // Known x-coordinates (sorted)
static double *y_values;     // Known y-coordinates
static double *query_points; // x-coordinates to interpolate
static double total_sum;      // Accumulator for the result

// --- Helper Functions ---

// Generates a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <num_queries> <seed>\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    num_queries = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_points <= 1 || num_queries <= 0) {
        fprintf(stderr, "FATAL: num_points must be > 1 and num_queries must be > 0.\n");
        exit(1);
    }

    mt_seed(seed);

    x_values = (double *)malloc(num_points * sizeof(double));
    y_values = (double *)malloc(num_points * sizeof(double));
    query_points = (double *)malloc(num_queries * sizeof(double));

    if (!x_values || !y_values || !query_points) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate sorted x_values and corresponding random y_values
    x_values[0] = rand_double(); // Start with a small random value
    y_values[0] = rand_double() * 100.0;
    for (int i = 1; i < num_points; ++i) {
        x_values[i] = x_values[i - 1] + 1.0 + rand_double(); // Ensure they are sorted and spaced
        y_values[i] = rand_double() * 100.0;
    }

    // Generate query points within the range of x_values
    double x_min = x_values[0];
    double x_max = x_values[num_points - 1];
    for (int i = 0; i < num_queries; ++i) {
        query_points[i] = x_min + rand_double() * (x_max - x_min);
    }
}

void run_computation() {
    total_sum = 0.0;
    for (int i = 0; i < num_queries; ++i) {
        double x_q = query_points[i];

        // Binary search to find the index 'j' such that x_values[j] <= x_q < x_values[j+1]
        int low = 0, high = num_points - 1;
        int j = 0; // lower bound index

        while (low <= high) {
            int mid = low + (high - low) / 2;
            if (x_values[mid] <= x_q) {
                j = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }

        // Ensure the upper bound is valid
        if (j >= num_points - 1) {
            j = num_points - 2;
        }

        double x0 = x_values[j];
        double y0 = y_values[j];
        double x1 = x_values[j + 1];
        double y1 = y_values[j + 1];

        // Linear Interpolation: y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
        double interpolated_y = y0 + (x_q - x0) * (y1 - y0) / (x1 - x0);

        total_sum += interpolated_y;
    }
}

void cleanup() {
    free(x_values);
    free(y_values);
    free(query_points);
}

// --- Main Driver ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout to prevent dead code elimination
    printf("%f\n", total_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
