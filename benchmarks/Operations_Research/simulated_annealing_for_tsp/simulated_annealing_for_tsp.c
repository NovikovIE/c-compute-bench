#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Mersenne Twister (MT19937) generator - DO NOT MODIFY
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
// End of Mersenne Twister

// Benchmark data structures
typedef struct {
    double x;
    double y;
} Point;

typedef struct {
    // Parameters
    int num_cities;
    long long num_iterations;
    double initial_temperature;
    uint32_t seed;
    
    // Data
    Point* cities;
    int* tour;
    double** dist_matrix;
    
    // Result
    double final_cost;
} BenchmarkData;

// Global data handle
BenchmarkData g_data;

// Utility to get random double [0, 1)
double random_double() {
    return (double)mt_rand() / (double)(UINT32_MAX);
}

// Utility to calculate initial total distance of the tour
double calculate_tour_cost(const int* tour, int num_cities, double** dist_matrix) {
    double cost = 0.0;
    for (int i = 0; i < num_cities - 1; ++i) {
        cost += dist_matrix[tour[i]][tour[i + 1]];
    }
    cost += dist_matrix[tour[num_cities - 1]][tour[0]]; // Return to start
    return cost;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_cities> <num_iterations> <initial_temperature> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_cities = atoi(argv[1]);
    g_data.num_iterations = atoll(argv[2]);
    g_data.initial_temperature = atof(argv[3]);
    g_data.seed = (uint32_t)atoi(argv[4]);

    if (g_data.num_cities <= 1 || g_data.num_iterations <= 0 || g_data.initial_temperature <= 0) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }
    
    mt_seed(g_data.seed);

    // Allocate memory
    g_data.cities = (Point*)malloc(g_data.num_cities * sizeof(Point));
    g_data.tour = (int*)malloc(g_data.num_cities * sizeof(int));
    g_data.dist_matrix = (double**)malloc(g_data.num_cities * sizeof(double*));
    if (!g_data.cities || !g_data.tour || !g_data.dist_matrix) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }
    for (int i = 0; i < g_data.num_cities; ++i) {
        g_data.dist_matrix[i] = (double*)malloc(g_data.num_cities * sizeof(double));
        if(!g_data.dist_matrix[i]) {
            fprintf(stderr, "Memory allocation failed.\n");
            exit(1);
        }
    }

    // Generate city coordinates
    for (int i = 0; i < g_data.num_cities; ++i) {
        g_data.cities[i].x = random_double() * 1000.0;
        g_data.cities[i].y = random_double() * 1000.0;
    }

    // Pre-calculate distance matrix
    for (int i = 0; i < g_data.num_cities; ++i) {
        for (int j = 0; j < g_data.num_cities; ++j) {
            if (i == j) {
                g_data.dist_matrix[i][j] = 0.0;
            } else {
                double dx = g_data.cities[i].x - g_data.cities[j].x;
                double dy = g_data.cities[i].y - g_data.cities[j].y;
                g_data.dist_matrix[i][j] = hypot(dx, dy);
            }
        }
    }

    // Initialize tour as a simple sequence
    for (int i = 0; i < g_data.num_cities; ++i) {
        g_data.tour[i] = i;
    }
    
    g_data.final_cost = 0.0;
}

void run_computation() {
    int n_cities = g_data.num_cities;
    double temp = g_data.initial_temperature;
    double current_cost = calculate_tour_cost(g_data.tour, n_cities, g_data.dist_matrix);
    
    // A geometric cooling schedule. The factor is calculated once.
    double cooling_factor = pow(0.001, 1.0 / g_data.num_iterations);

    for (long long k = 0; k < g_data.num_iterations; ++k) {
        // Generate two distinct random indices for 2-opt swap
        int i, j;
        do {
            i = mt_rand() % n_cities;
            j = mt_rand() % n_cities;
        } while (i == j);

        if (i > j) { int temp_idx = i; i = j; j = temp_idx; }

        // Calculate the change in tour cost if we reverse the segment [i, j]
        int i_prev = (i == 0) ? n_cities - 1 : i - 1;
        int j_next = (j == n_cities - 1) ? 0 : j + 1;

        int city_i_prev = g_data.tour[i_prev];
        int city_i = g_data.tour[i];
        int city_j = g_data.tour[j];
        int city_j_next = g_data.tour[j_next];

        double removed_cost = g_data.dist_matrix[city_i_prev][city_i] + g_data.dist_matrix[city_j][city_j_next];
        double added_cost = g_data.dist_matrix[city_i_prev][city_j] + g_data.dist_matrix[city_i][city_j_next];
        double delta_cost = added_cost - removed_cost;
       
        // Acceptance probability
        if (delta_cost < 0 || random_double() < exp(-delta_cost / temp)) {
            // Accept the new tour: reverse the segment from i to j
            int low = i;
            int high = j;
            while (low < high) {
                int temp_city = g_data.tour[low];
                g_data.tour[low] = g_data.tour[high];
                g_data.tour[high] = temp_city;
                low++;
                high--;
            }
            current_cost += delta_cost;
        }
        
        // Cool down the temperature
        temp *= cooling_factor;
    }

    g_data.final_cost = current_cost;
}

void cleanup() {
    for (int i = 0; i < g_data.num_cities; ++i) {
        free(g_data.dist_matrix[i]);
    }
    free(g_data.dist_matrix);
    free(g_data.tour);
    free(g_data.cities);
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
    printf("%f\n", g_data.final_cost);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
