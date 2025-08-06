#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

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
int num_cities;
int **distance_matrix;
long long min_distance; // Use long long to avoid overflow on total path distance

// --- HELPER FUNCTIONS FOR COMPUTATION ---
static void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Recursive helper for run_computation. Finds the shortest path via brute-force permutation.
static void solve_tsp_recursive(int* path, int start_index, int n_permute) {
    // Base case: a full permutation is generated (all cities from 1 to N-1 have been placed)
    if (start_index == n_permute) {
        long long current_distance = 0;
        
        // Distance from starting city (0) to the first city in the permutation
        current_distance += distance_matrix[0][path[0]];
        
        // Distance between cities in the permutation
        for (int i = 0; i < n_permute - 1; i++) {
            current_distance += distance_matrix[path[i]][path[i+1]];
        }
        
        // Distance from the last city in the permutation back to the starting city (0)
        current_distance += distance_matrix[path[n_permute - 1]][0];
        
        if (current_distance < min_distance) {
            min_distance = current_distance;
        }
        return;
    }

    // Recursive step: generate permutations
    for (int i = start_index; i < n_permute; i++) {
        swap(&path[start_index], &path[i]);
        solve_tsp_recursive(path, start_index + 1, n_permute);
        swap(&path[start_index], &path[i]); // Backtrack
    }
}

// --- BENCHMARK CORE FUNCTIONS ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_cities> <seed>\n", argv[0]);
        exit(1);
    }

    num_cities = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_cities < 2) {
        fprintf(stderr, "Error: num_cities must be at least 2 for a valid tour.\n");
        exit(1);
    }
    
    if (num_cities > 15) {
        fprintf(stderr, "Warning: num_cities > 15 is computationally very expensive.\n");
    }

    mt_seed(seed);

    int *x_coords = (int*)malloc(num_cities * sizeof(int));
    int *y_coords = (int*)malloc(num_cities * sizeof(int));
    if (!x_coords || !y_coords) {
        fprintf(stderr, "FATAL: Failed to allocate memory for coordinates.\n");
        exit(1);
    }
    for (int i = 0; i < num_cities; i++) {
        x_coords[i] = mt_rand() % 10000;
        y_coords[i] = mt_rand() % 10000;
    }

    distance_matrix = (int**)malloc(num_cities * sizeof(int*));
    if (!distance_matrix) {
        fprintf(stderr, "FATAL: Failed to allocate memory for distance matrix.\n");
        exit(1);
    }
    for (int i = 0; i < num_cities; i++) {
        distance_matrix[i] = (int*)malloc(num_cities * sizeof(int));
        if (!distance_matrix[i]) {
            fprintf(stderr, "FATAL: Failed to allocate memory for distance matrix column.\n");
            exit(1);
        }
    }

    for (int i = 0; i < num_cities; i++) {
        for (int j = 0; j < num_cities; j++) {
            if (i == j) {
                distance_matrix[i][j] = 0;
            } else {
                long long dx = (long long)x_coords[i] - x_coords[j];
                long long dy = (long long)y_coords[i] - y_coords[j];
                distance_matrix[i][j] = (int)(dx*dx + dy*dy);
            }
        }
    }

    free(x_coords);
    free(y_coords);
}

void run_computation() {
    int n_permute = num_cities - 1;
    if (n_permute <= 0) {
        min_distance = 0;
        return;
    }

    int* path = (int*)malloc(n_permute * sizeof(int));
    if (!path) {
        fprintf(stderr, "FATAL: Failed to allocate path array in computation.\n");
        exit(1);
    }
    for (int i = 0; i < n_permute; i++) {
        path[i] = i + 1;
    }

    min_distance = LLONG_MAX;
    solve_tsp_recursive(path, 0, n_permute);
    
    free(path);
}

void cleanup() {
    if (distance_matrix) {
        for (int i = 0; i < num_cities; i++) {
            if (distance_matrix[i]) {
                free(distance_matrix[i]);
            }
        }
        free(distance_matrix);
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

    printf("%lld\n", min_distance); 
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
