/*
 * BENCHMARK: Ant Colony Optimization for the Traveling Salesperson Problem (TSP)
 * 
 * This program implements the Ant Colony Optimization (ACO) metaheuristic to find an
 * approximate solution to the Traveling Salesperson Problem. The goal of the TSP
 * is to find the shortest possible route that visits each city in a given list
 * exactly once and returns to the origin city.
 * 
 * The benchmark is structured as follows:
 * 1. setup_benchmark: Reads parameters, allocates memory, and generates the problem
 *    instance. This includes creating a set of random 2D city coordinates and 
 *    pre-calculating a distance matrix between all pairs of cities.
 * 2. run_computation: Executes the core ACO algorithm. It simulates multiple
 *    generations (iterations) of ants building tours. In each iteration:
 *    a. Each ant constructs a tour by probabilistically choosing the next city
 *       based on pheromone levels and distance (heuristic information).
 *    b. The pheromone trails are updated: they evaporate over time and are
 *       reinforced by ants based on the quality of their tours.
 *    c. The all-time best tour found is recorded.
 * 3. cleanup: Frees all allocated memory.
 * 
 * The final result is the length of the best tour found.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>

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

// ACO algorithm constants
#define ALPHA 1.0  // Pheromone influence
#define BETA 2.0   // Heuristic (distance) influence
#define INITIAL_PHEROMONE 0.1

// Type definitions and global data structure
typedef struct {
    int x, y;
} City;

typedef struct {
    // Parameters
    int num_ants;
    int num_cities;
    int num_iterations;
    double pheromone_evaporation_rate;

    // Problem Data
    City* cities;
    double** dist_matrix;
    double** pheromone_matrix;
    double** next_city_probabilities; // Pre-allocated for ants

    // Ant-specific Data
    int** ant_tours;
    int* ant_visited_flags;

    // Result
    double best_tour_length;
} BenchmarkData;

static BenchmarkData g_data;

// Helper to get a random double in [0, 1)
double mt_rand_double() {
    return (double)mt_rand() / ((double)UINT32_MAX + 1.0);
}

// Helper to calculate Euclidean distance
double calculate_distance(const City* c1, const City* c2) {
    double dx = c1->x - c2->x;
    double dy = c1->y - c2->y;
    return sqrt(dx * dx + dy * dy);
}

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_ants> <num_cities> <num_iterations> <evaporation_rate> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_ants = atoi(argv[1]);
    g_data.num_cities = atoi(argv[2]);
    g_data.num_iterations = atoi(argv[3]);
    g_data.pheromone_evaporation_rate = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);
    mt_seed(seed);

    // Allocate memory
    g_data.cities = (City*)malloc(g_data.num_cities * sizeof(City));
    
    g_data.dist_matrix = (double**)malloc(g_data.num_cities * sizeof(double*));
    g_data.pheromone_matrix = (double**)malloc(g_data.num_cities * sizeof(double*));
    for (int i = 0; i < g_data.num_cities; ++i) {
        g_data.dist_matrix[i] = (double*)malloc(g_data.num_cities * sizeof(double));
        g_data.pheromone_matrix[i] = (double*)malloc(g_data.num_cities * sizeof(double));
    }

    g_data.ant_tours = (int**)malloc(g_data.num_ants * sizeof(int*));
    for (int i = 0; i < g_data.num_ants; ++i) {
        g_data.ant_tours[i] = (int*)malloc((g_data.num_cities + 1) * sizeof(int));
    }

    g_data.ant_visited_flags = (int*)malloc(g_data.num_cities * sizeof(int));
    g_data.next_city_probabilities = (double**)malloc(g_data.num_ants * sizeof(double*));
    for (int i = 0; i < g_data.num_ants; ++i) {
        g_data.next_city_probabilities[i] = (double*)malloc(g_data.num_cities * sizeof(double));
    }

    // Initialize data
    for (int i = 0; i < g_data.num_cities; ++i) {
        g_data.cities[i].x = mt_rand() % 1000;
        g_data.cities[i].y = mt_rand() % 1000;
    }

    for (int i = 0; i < g_data.num_cities; ++i) {
        for (int j = 0; j < g_data.num_cities; ++j) {
            if (i == j) {
                g_data.dist_matrix[i][j] = 0.0;
            } else {
                g_data.dist_matrix[i][j] = calculate_distance(&g_data.cities[i], &g_data.cities[j]);
            }
            g_data.pheromone_matrix[i][j] = INITIAL_PHEROMONE;
        }
    }

    g_data.best_tour_length = DBL_MAX;
}

void run_computation() {
    double* ant_tour_lengths = (double*)malloc(g_data.num_ants * sizeof(double));

    for (int iter = 0; iter < g_data.num_iterations; ++iter) {
        // For each ant, construct a tour
        for (int ant = 0; ant < g_data.num_ants; ++ant) {
            memset(g_data.ant_visited_flags, 0, g_data.num_cities * sizeof(int));
            int start_city = mt_rand() % g_data.num_cities;
            g_data.ant_tours[ant][0] = start_city;
            g_data.ant_visited_flags[start_city] = 1;
            int current_city = start_city;

            for (int step = 1; step < g_data.num_cities; ++step) {
                double total_attractiveness = 0.0;

                // Calculate probability for moving to each unvisited city
                for (int next_city = 0; next_city < g_data.num_cities; ++next_city) {
                    if (!g_data.ant_visited_flags[next_city]) {
                        double pheromone = g_data.pheromone_matrix[current_city][next_city];
                        double dist = g_data.dist_matrix[current_city][next_city];
                        double attractiveness = pow(pheromone, ALPHA) * pow(1.0 / dist, BETA);
                        g_data.next_city_probabilities[ant][next_city] = attractiveness;
                        total_attractiveness += attractiveness;
                    } else {
                        g_data.next_city_probabilities[ant][next_city] = 0.0;
                    }
                }

                // Select next city using roulette wheel selection
                double rand_val = mt_rand_double() * total_attractiveness;
                double cumulative_prob = 0.0;
                int chosen_city = -1;
                for (int next_city = 0; next_city < g_data.num_cities; ++next_city) {
                    if (!g_data.ant_visited_flags[next_city]) {
                        cumulative_prob += g_data.next_city_probabilities[ant][next_city];
                        if (rand_val <= cumulative_prob) {
                            chosen_city = next_city;
                            break;
                        }
                    }
                }

                if (chosen_city == -1) { // Failsafe for floating point inaccuracies
                    for (int k = g_data.num_cities - 1; k >= 0; --k) {
                        if (!g_data.ant_visited_flags[k]) {
                            chosen_city = k;
                            break;
                        }
                    }
                }

                g_data.ant_tours[ant][step] = chosen_city;
                g_data.ant_visited_flags[chosen_city] = 1;
                current_city = chosen_city;
            }
            g_data.ant_tours[ant][g_data.num_cities] = start_city; // Return to start

            // Calculate tour length
            double tour_length = 0.0;
            for (int i = 0; i < g_data.num_cities; ++i) {
                tour_length += g_data.dist_matrix[g_data.ant_tours[ant][i]][g_data.ant_tours[ant][i + 1]];
            }
            ant_tour_lengths[ant] = tour_length;

            if (tour_length < g_data.best_tour_length) {
                g_data.best_tour_length = tour_length;
            }
        }

        // Pheromone evaporation
        for (int i = 0; i < g_data.num_cities; ++i) {
            for (int j = 0; j < g_data.num_cities; ++j) {
                g_data.pheromone_matrix[i][j] *= (1.0 - g_data.pheromone_evaporation_rate);
            }
        }

        // Pheromone deposition
        for (int ant = 0; ant < g_data.num_ants; ++ant) {
            double pheromone_deposit = 1.0 / ant_tour_lengths[ant];
            for (int i = 0; i < g_data.num_cities; ++i) {
                int from = g_data.ant_tours[ant][i];
                int to = g_data.ant_tours[ant][i + 1];
                g_data.pheromone_matrix[from][to] += pheromone_deposit;
                g_data.pheromone_matrix[to][from] += pheromone_deposit; // Symmetric TSP
            }
        }
    }

    free(ant_tour_lengths);
}

void cleanup() {
    free(g_data.cities);

    for (int i = 0; i < g_data.num_cities; ++i) {
        free(g_data.dist_matrix[i]);
        free(g_data.pheromone_matrix[i]);
    }
    free(g_data.dist_matrix);
    free(g_data.pheromone_matrix);

    for (int i = 0; i < g_data.num_ants; ++i) {
        free(g_data.ant_tours[i]);
        free(g_data.next_city_probabilities[i]);
    }
    free(g_data.ant_tours);
    free(g_data.next_city_probabilities);

    free(g_data.ant_visited_flags);
}


int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%d\n", (int)g_data.best_tour_length);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
