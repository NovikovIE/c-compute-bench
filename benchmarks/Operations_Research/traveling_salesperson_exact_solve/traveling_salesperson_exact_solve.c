#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA AND PARAMETERS ---
int num_cities;
int** dist_matrix; // Adjacency matrix for distances
int min_path_cost; // The final result

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_cities> <seed>\n", argv[0]);
        exit(1);
    }

    num_cities = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);
    mt_seed(seed);

    // Allocate the distance matrix
    dist_matrix = (int**)malloc(num_cities * sizeof(int*));
    if (dist_matrix == NULL) {
        perror("Failed to allocate memory for distance matrix rows");
        exit(1);
    }
    for (int i = 0; i < num_cities; ++i) {
        dist_matrix[i] = (int*)malloc(num_cities * sizeof(int));
        if (dist_matrix[i] == NULL) {
            perror("Failed to allocate memory for distance matrix columns");
            exit(1);
        }
    }

    // Populate the distance matrix with random, symmetric distances
    for (int i = 0; i < num_cities; ++i) {
        for (int j = i; j < num_cities; ++j) {
            if (i == j) {
                dist_matrix[i][j] = 0;
            } else {
                int dist = (mt_rand() % 100) + 1; // Random distance between 1 and 100
                dist_matrix[i][j] = dist;
                dist_matrix[j][i] = dist;
            }
        }
    }
}

// Helper function for the recursive exact solver
static void solve_recursive(int current_city, int count, int current_cost, int* visited) {
    // Pruning step (Branch and Bound): if current path is already worse than best found, stop.
    if (current_cost >= min_path_cost) {
        return;
    }

    // Base case: all cities have been visited
    if (count == num_cities) {
        // Add the cost of returning to the starting city (city 0)
        int total_cost = current_cost + dist_matrix[current_city][0];
        if (total_cost < min_path_cost) {
            min_path_cost = total_cost;
        }
        return;
    }

    // Recursive step: explore unvisited neighbors
    for (int next_city = 0; next_city < num_cities; ++next_city) {
        if (!visited[next_city]) {
            visited[next_city] = 1; // Mark as visited
            solve_recursive(next_city, count + 1, current_cost + dist_matrix[current_city][next_city], visited);
            visited[next_city] = 0; // Backtrack
        }
    }
}

void run_computation() {
    min_path_cost = INT_MAX;

    // Keep track of visited cities for the current path
    int* visited = (int*)calloc(num_cities, sizeof(int));
    if (!visited) {
        perror("Failed to allocate memory for visited array");
        exit(1);
    }

    // Start the tour from city 0
    visited[0] = 1;
    
    solve_recursive(0, 1, 0, visited);

    free(visited);
}

void cleanup() {
    for (int i = 0; i < num_cities; ++i) {
        free(dist_matrix[i]);
    }
    free(dist_matrix);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (minimum path cost) to stdout
    printf("%d\n", min_path_cost);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
