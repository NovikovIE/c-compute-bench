#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

// --- START MERSENNE TWISTER ---
// Note: stdio.h, stdlib.h, stdint.h, time.h should be included at the top of the file
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

// --- BENCHMARK DATA AND GLOBALS ---
typedef struct {
    double birth;
    double death;
} Point;

int num_points_diagram_a;
int num_points_diagram_b;
Point *diagram_a = NULL;
Point *diagram_b = NULL;
double final_result = 0.0;

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points_diagram_a> <num_points_diagram_b> <seed>\n", argv[0]);
        exit(1);
    }

    num_points_diagram_a = atoi(argv[1]);
    num_points_diagram_b = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_points_diagram_a <= 0 || num_points_diagram_b <= 0) {
        fprintf(stderr, "Number of points must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    diagram_a = (Point*)malloc(num_points_diagram_a * sizeof(Point));
    diagram_b = (Point*)malloc(num_points_diagram_b * sizeof(Point));
    if (!diagram_a || !diagram_b) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate points for diagram A
    for (int i = 0; i < num_points_diagram_a; ++i) {
        double birth = ((double)mt_rand() / (double)UINT32_MAX) * 1000.0;
        double persistence = ((double)mt_rand() / (double)UINT32_MAX) * 100.0 + 0.1; // Ensure persistence > 0
        diagram_a[i].birth = birth;
        diagram_a[i].death = birth + persistence;
    }

    // Generate points for diagram B
    for (int i = 0; i < num_points_diagram_b; ++i) {
        double birth = ((double)mt_rand() / (double)UINT32_MAX) * 1000.0;
        double persistence = ((double)mt_rand() / (double)UINT32_MAX) * 100.0 + 0.1;
        diagram_b[i].birth = birth;
        diagram_b[i].death = birth + persistence;
    }
}

void run_computation() {
    double total_cost = 0.0;
    
    // For each point in diagram A, find its closest neighbor in diagram B.
    // This is an O(N*M) calculation representative of constructing a cost matrix
    // in Wasserstein distance algorithms. It computes a one-sided Chamfer-like distance.
    for (int i = 0; i < num_points_diagram_a; ++i) {
        double min_dist = DBL_MAX;
        
        Point pa = diagram_a[i];

        for (int j = 0; j < num_points_diagram_b; ++j) {
            Point pb = diagram_b[j];
            
            // L-infinity distance is common for persistence diagrams
            double dist = fmax(fabs(pa.birth - pb.birth), fabs(pa.death - pb.death));
            
            if (dist < min_dist) {
                min_dist = dist;
            }
        }
        total_cost += min_dist;
    }
    
    final_result = total_cost;
}

void cleanup() {
    if (diagram_a) free(diagram_a);
    if (diagram_b) free(diagram_b);
}

// --- MAIN FUNCTION ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
