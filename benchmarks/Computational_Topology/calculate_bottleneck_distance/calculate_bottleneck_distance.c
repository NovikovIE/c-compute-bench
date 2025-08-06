#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// A point in a persistence diagram is defined by its birth and death times.
typedef struct {
    double birth;
    double death;
} Point;

// --- Global Variables for Benchmark Data ---
static int NUM_POINTS_A;
static int NUM_POINTS_B;
static Point* diagram_a;
static Point* diagram_b;
static double final_bottleneck_distance;

// Generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points_diagram_a> <num_points_diagram_b> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_POINTS_A = atoi(argv[1]);
    NUM_POINTS_B = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    diagram_a = (Point*)malloc(NUM_POINTS_A * sizeof(Point));
    diagram_b = (Point*)malloc(NUM_POINTS_B * sizeof(Point));

    if (!diagram_a || !diagram_b) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Generate points for the first persistence diagram
    for (int i = 0; i < NUM_POINTS_A; ++i) {
        double birth = rand_double();
        double persistence = rand_double() * 0.5; // Keep persistence relatively small
        diagram_a[i].birth = birth;
        diagram_a[i].death = birth + persistence; 
    }

    // Generate points for the second persistence diagram
    for (int i = 0; i < NUM_POINTS_B; ++i) {
        double birth = rand_double();
        double persistence = rand_double() * 0.5;
        diagram_b[i].birth = birth;
        diagram_b[i].death = birth + persistence;
    }
}

void run_computation() {
    // This benchmark computes the Hausdorff distance between the two persistence diagrams,
    // which is a core component of many bottleneck distance algorithms.
    // The bottleneck distance is inf_eta sup_x d(x, eta(x)), where eta is a partial bijection.
    // The Hausdorff distance is max(sup_x inf_y d(x,y), sup_y inf_x d(x,y)).
    // The distance between points is the L-infinity norm.
    // The cost of matching a point to the diagonal is (death - birth) / 2.

    double max_of_min_dists_a = 0.0;
    // For each point in A, find the minimum cost to match it to a point in B or the diagonal.
    for (int i = 0; i < NUM_POINTS_A; ++i) {
        double cost_to_diagonal = (diagram_a[i].death - diagram_a[i].birth) / 2.0;
        double min_dist_for_a_i = cost_to_diagonal;

        for (int j = 0; j < NUM_POINTS_B; ++j) {
            double dist = fmax(fabs(diagram_a[i].birth - diagram_b[j].birth),
                               fabs(diagram_a[i].death - diagram_b[j].death));
            min_dist_for_a_i = fmin(min_dist_for_a_i, dist);
        }
        max_of_min_dists_a = fmax(max_of_min_dists_a, min_dist_for_a_i);
    }

    double max_of_min_dists_b = 0.0;
    // For each point in B, find the minimum cost to match it to a point in A or the diagonal.
    for (int j = 0; j < NUM_POINTS_B; ++j) {
        double cost_to_diagonal = (diagram_b[j].death - diagram_b[j].birth) / 2.0;
        double min_dist_for_b_j = cost_to_diagonal;

        for (int i = 0; i < NUM_POINTS_A; ++i) {
            double dist = fmax(fabs(diagram_b[j].birth - diagram_a[i].birth),
                               fabs(diagram_b[j].death - diagram_a[i].death));
            min_dist_for_b_j = fmin(min_dist_for_b_j, dist);
        }
        max_of_min_dists_b = fmax(max_of_min_dists_b, min_dist_for_b_j);
    }

    // The Hausdorff distance is the maximum of these two values.
    final_bottleneck_distance = fmax(max_of_min_dists_a, max_of_min_dists_b);
}

void cleanup() {
    free(diagram_a);
    free(diagram_b);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%f\n", final_bottleneck_distance);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
