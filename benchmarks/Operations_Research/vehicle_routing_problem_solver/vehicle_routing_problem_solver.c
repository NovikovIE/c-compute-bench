#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

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

// --- BENCHMARK DATA AND STRUCTURES ---

// Location of a customer or depot
typedef struct {
    int x, y;
} Point;

// A saving achieved by merging two routes
typedef struct {
    int customer_i;
    int customer_j;
    double value;
} Saving;

// Global struct to hold all benchmark data
struct {
    int num_customers;
    int num_vehicles;
    int capacity_per_vehicle;
    long final_total_distance;

    Point* locations;       // Size: num_customers + 1 (0 is depot)
    int* demands;           // Size: num_customers + 1
    double** dist_matrix;   // Size: (N+1)x(N+1)

    // Data for computation phase
    Saving* savings_list;
    int* next_node;
    int* prev_node;
    long* route_demands; 

} g_data;

// Comparison function for qsort to sort savings in descending order
int compare_savings(const void* a, const void* b) {
    Saving* s1 = (Saving*)a;
    Saving* s2 = (Saving*)b;
    if (s2->value > s1->value) return 1;
    if (s2->value < s1->value) return -1;
    return 0;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_customers> <num_vehicles> <capacity_per_vehicle> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_customers = atoi(argv[1]);
    g_data.num_vehicles = atoi(argv[2]);
    g_data.capacity_per_vehicle = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    mt_seed(seed);

    int total_locations = g_data.num_customers + 1;
    g_data.final_total_distance = 0;

    // Allocate memory
    g_data.locations = (Point*)malloc(total_locations * sizeof(Point));
    g_data.demands = (int*)malloc(total_locations * sizeof(int));
    g_data.dist_matrix = (double**)malloc(total_locations * sizeof(double*));
    for (int i = 0; i < total_locations; ++i) {
        g_data.dist_matrix[i] = (double*)malloc(total_locations * sizeof(double));
    }

    size_t num_savings = (size_t)g_data.num_customers * (g_data.num_customers - 1) / 2;
    g_data.savings_list = (Saving*) malloc(num_savings * sizeof(Saving));
    g_data.next_node = (int*) malloc(total_locations * sizeof(int));
    g_data.prev_node = (int*) malloc(total_locations * sizeof(int));
    g_data.route_demands = (long*) malloc(total_locations * sizeof(long));

    // Generate random locations and demands
    // Depot at index 0
    g_data.locations[0].x = mt_rand() % 1000;
    g_data.locations[0].y = mt_rand() % 1000;
    g_data.demands[0] = 0;

    for (int i = 1; i < total_locations; ++i) {
        g_data.locations[i].x = mt_rand() % 1000;
        g_data.locations[i].y = mt_rand() % 1000;
        g_data.demands[i] = (mt_rand() % 30) + 1; // Demand from 1 to 30
    }

    // Pre-calculate distance matrix (Euclidean distance)
    for (int i = 0; i < total_locations; ++i) {
        for (int j = 0; j < total_locations; ++j) {
            double dx = g_data.locations[i].x - g_data.locations[j].x;
            double dy = g_data.locations[i].y - g_data.locations[j].y;
            g_data.dist_matrix[i][j] = sqrt(dx * dx + dy * dy);
        }
    }
}

void run_computation() {
    int n_cust = g_data.num_customers;
    int capacity = g_data.capacity_per_vehicle;

    // Phase 1: Calculate savings using the Clarke and Wright savings algorithm.
    // Saving s(i,j) = d(depot,i) + d(depot,j) - d(i,j)
    size_t savings_idx = 0;
    for (int i = 1; i <= n_cust; ++i) {
        for (int j = i + 1; j <= n_cust; ++j) {
            double s_val = g_data.dist_matrix[0][i] + g_data.dist_matrix[0][j] - g_data.dist_matrix[i][j];
            g_data.savings_list[savings_idx].customer_i = i;
            g_data.savings_list[savings_idx].customer_j = j;
            g_data.savings_list[savings_idx].value = s_val;
            savings_idx++;
        }
    }

    // Phase 2: Sort savings in descending order.
    size_t num_savings = (size_t)n_cust * (n_cust - 1) / 2;
    qsort(g_data.savings_list, num_savings, sizeof(Saving), compare_savings);

    // Phase 3: Build routes by merging based on savings.
    // Initially, each customer is in their own route: Depot -> i -> Depot
    for(int i = 0; i <= n_cust; i++) {
        g_data.next_node[i] = 0; // 0 represents the depot
        g_data.prev_node[i] = 0;
        g_data.route_demands[i] = g_data.demands[i]; // Customer i is the head of its own route
    }

    for (size_t i = 0; i < num_savings; ++i) {
        int cust_i = g_data.savings_list[i].customer_i;
        int cust_j = g_data.savings_list[i].customer_j;

        // An customer is an endpoint if it's connected to the depot.
        // If g_data.next_node[cust_i] is 0, cust_i is the end of a route.
        // If g_data.prev_node[cust_j] is 0, cust_j is the start of a route.
        if (g_data.next_node[cust_i] == 0 && g_data.prev_node[cust_j] == 0) {

            // Find the starting customer (route ID) for both routes to check for cycles and demand.
            int head_i = cust_i; while(g_data.prev_node[head_i] != 0) head_i = g_data.prev_node[head_i];
            int head_j = cust_j; while(g_data.prev_node[head_j] != 0) head_j = g_data.prev_node[head_j];

            if (head_i != head_j) { // Not already in the same route
                if (g_data.route_demands[head_i] + g_data.route_demands[head_j] <= capacity) {
                    // Merge routes: connect i -> j
                    g_data.next_node[cust_i] = cust_j;
                    g_data.prev_node[cust_j] = cust_i;

                    // Consolidate demand into head_i's route
                    g_data.route_demands[head_i] += g_data.route_demands[head_j];
                    g_data.route_demands[head_j] = 0; // Mark old route as empty
                }
            }
        }
    }

    // Phase 4: Calculate final total distance of all constructed routes.
    long total_dist = 0;
    int* visited = (int*)calloc(n_cust + 1, sizeof(int));
    for (int i = 1; i <= n_cust; ++i) {
        if (!visited[i] && g_data.prev_node[i] == 0) { // Found the start of a route
            int curr = i;
            total_dist += g_data.dist_matrix[0][curr]; // Depot -> start
            while (curr != 0) {
                visited[curr] = 1;
                int next = g_data.next_node[curr];
                if (next != 0) {
                    total_dist += g_data.dist_matrix[curr][next];
                } else {
                    total_dist += g_data.dist_matrix[curr][0]; // end -> Depot
                }
                curr = next;
            }
        }
    }
    g_data.final_total_distance = total_dist;
    free(visited);
}

void cleanup() {
    int total_locations = g_data.num_customers + 1;
    for (int i = 0; i < total_locations; ++i) {
        free(g_data.dist_matrix[i]);
    }
    free(g_data.dist_matrix);
    free(g_data.locations);
    free(g_data.demands);
    
    free(g_data.savings_list);
    free(g_data.next_node);
    free(g_data.prev_node);
    free(g_data.route_demands);
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
    printf("%ld\n", g_data.final_total_distance);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
