#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Verbatim Mersenne Twister ---
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
// --- End Verbatim Mersenne Twister ---

// Benchmark-specific data structures
typedef struct {
    double x;
    double y;
} Point;

// Global parameters and data pointers
static int num_points;
static int num_queries;
static double power_parameter;

static Point *data_points;
static double *data_values;
static Point *query_points;
static double *query_results;

static double final_result_accumulator;

// Generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_points> <num_queries> <power_parameter> <seed>\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    num_queries = atoi(argv[2]);
    power_parameter = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    data_points = (Point*)malloc(num_points * sizeof(Point));
    data_values = (double*)malloc(num_points * sizeof(double));
    query_points = (Point*)malloc(num_queries * sizeof(Point));
    query_results = (double*)malloc(num_queries * sizeof(double));

    if (!data_points || !data_values || !query_points || !query_results) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate known data points and their associated values
    for (int i = 0; i < num_points; i++) {
        data_points[i].x = rand_double();
        data_points[i].y = rand_double();
        data_values[i] = rand_double();
    }

    // Generate query points for which we will interpolate values
    for (int i = 0; i < num_queries; i++) {
        query_points[i].x = rand_double();
        query_points[i].y = rand_double();
    }
}

void run_computation() {
    final_result_accumulator = 0.0;

    for (int i = 0; i < num_queries; i++) {
        double numerator = 0.0;
        double denominator = 0.0;
        Point q = query_points[i];
        int exact_match = 0;

        for (int j = 0; j < num_points; j++) {
            Point p = data_points[j];
            double dx = p.x - q.x;
            double dy = p.y - q.y;
            double dist_sq = dx * dx + dy * dy;

            // Handle case where query point is identical to a data point
            if (dist_sq == 0.0) {
                query_results[i] = data_values[j];
                exact_match = 1;
                break;
            }

            double dist = sqrt(dist_sq);
            double weight = 1.0 / pow(dist, power_parameter);
            
            numerator += weight * data_values[j];
            denominator += weight;
        }

        if (!exact_match) {
            if (denominator == 0.0) {
                query_results[i] = 0.0; 
            } else {
                query_results[i] = numerator / denominator;
            }
        }
        
        final_result_accumulator += query_results[i];
    }
}

void cleanup() {
    free(data_points);
    free(data_values);
    free(query_points);
    free(query_results);
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
    printf("%f\n", final_result_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
