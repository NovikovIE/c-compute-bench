#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>

// DO NOT MODIFY - Mersenne Twister Generator
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
// End of Mersenne Twister Generator

// --- Benchmark Globals ---
// Parameters
int num_agents;

// Data structures
int **cost_matrix;
int *label_x, *label_y;
int *match_x, *match_y;
int *slack, *slack_x;
bool *S_set, *T_set;

// Final Result
int final_result;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "FATAL: requires 2 arguments. Usage: %s <num_agents> <seed>\n", argv[0]);
        exit(1);
    }

    num_agents = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (num_agents <= 0) {
        fprintf(stderr, "FATAL: num_agents must be a positive integer.\n");
        exit(1);
    }

    // Allocate cost matrix
    cost_matrix = (int **)malloc(num_agents * sizeof(int *));
    for (int i = 0; i < num_agents; i++) {
        cost_matrix[i] = (int *)malloc(num_agents * sizeof(int));
    }

    // Fill cost matrix with random values
    for (int i = 0; i < num_agents; i++) {
        for (int j = 0; j < num_agents; j++) {
            cost_matrix[i][j] = (mt_rand() % 1000) + 1; // Costs from 1 to 1000
        }
    }

    // Allocate auxiliary arrays for the Hungarian algorithm
    label_x = (int *)malloc(num_agents * sizeof(int));
    label_y = (int *)malloc(num_agents * sizeof(int));
    match_x = (int *)malloc(num_agents * sizeof(int));
    match_y = (int *)malloc(num_agents * sizeof(int));
    slack = (int *)malloc(num_agents * sizeof(int));
    slack_x = (int *)malloc(num_agents * sizeof(int));
    S_set = (bool *)malloc(num_agents * sizeof(bool));
    T_set = (bool *)malloc(num_agents * sizeof(bool));
}


void run_computation() {
    // Implementation of the Hungarian Algorithm (Kuhn-Munkres)
    // for min-cost perfect matching. Complexity: O(N^3).
    // This variant works on a minimization problem directly.
    
    // Step 1: Initial labeling (l_x[i] + l_y[j] <= cost[i][j])
    memset(label_y, 0, num_agents * sizeof(int));
    for (int i = 0; i < num_agents; i++) {
        label_x[i] = cost_matrix[i][0];
        for (int j = 1; j < num_agents; j++) {
            if (cost_matrix[i][j] < label_x[i]) {
                label_x[i] = cost_matrix[i][j];
            }
        }
    }

    // Step 2: Initial (empty) matching
    memset(match_x, -1, num_agents * sizeof(int));
    memset(match_y, -1, num_agents * sizeof(int));

    // Step 3: Find augmenting paths for each agent/worker x
    for (int x = 0; x < num_agents; x++) {
        
        memset(S_set, false, num_agents * sizeof(bool));
        memset(T_set, false, num_agents * sizeof(bool));
        
        S_set[x] = true; // Start the alternating tree with root x

        // Initialize slack for all tasks y based on the new root x
        // slack[y] = cost[i][y] - l_x[i] - l_y[y]
        for (int y = 0; y < num_agents; y++) {
            slack[y] = cost_matrix[x][y] - label_x[x] - label_y[y];
            slack_x[y] = x;
        }

        int y_leaf;
        while (true) {
            // Find y not in T with minimum slack
            int delta = INT_MAX;
            int y_min_slack = -1;
            for (int y = 0; y < num_agents; y++) {
                if (!T_set[y] && slack[y] < delta) {
                    delta = slack[y];
                    y_min_slack = y;
                }
            }

            // Update labels to introduce new edges into the equality subgraph
            for (int i = 0; i < num_agents; i++) {
                if (S_set[i]) label_x[i] += delta;
                if (T_set[i]) label_y[i] -= delta;
            }
            
            // Due to the label update, slacks for y not in T decrease
            for (int y = 0; y < num_agents; y++) {
                if (!T_set[y]) {
                    slack[y] -= delta;
                }
            }
            
            // y_min_slack is now in the equality subgraph. Add it to the alternating tree T.
            T_set[y_min_slack] = true;

            if (match_y[y_min_slack] == -1) { // Found an augmenting path to an unmatched y
                y_leaf = y_min_slack;
                break; // Exit while loop to augment
            }

            // Path not found yet, extend the alternating tree.
            // The matched partner of y_min_slack is added to S.
            int next_x = match_y[y_min_slack];
            S_set[next_x] = true;
            
            // Update slacks because a new x-vertex is in S
            for (int y = 0; y < num_agents; y++) {
                if (!T_set[y]) {
                    int new_slack = cost_matrix[next_x][y] - label_x[next_x] - label_y[y];
                    if (slack[y] > new_slack) {
                        slack[y] = new_slack;
                        slack_x[y] = next_x; // Track which x provides the best path to y
                    }
                }
            }
        }
        
        // Augment the path, flipping edges.
        // We backtrack from the unmatched y (y_leaf) using the slack_x pointers.
        while (y_leaf != -1) {
            int current_x = slack_x[y_leaf];
            int prev_y = match_x[current_x]; // The y previously matched with this x
            match_y[y_leaf] = current_x;
            match_x[current_x] = y_leaf;
            y_leaf = prev_y; // Move to the previous edge in the path
        }
    }

    // Calculate final cost from the optimal matching
    int total_cost = 0;
    for (int y = 0; y < num_agents; y++) {
        if (match_y[y] != -1) {
            total_cost += cost_matrix[match_y[y]][y];
        }
    }
    final_result = total_cost;
}

void cleanup() {
    for (int i = 0; i < num_agents; i++) {
        free(cost_matrix[i]);
    }
    free(cost_matrix);
    free(label_x);
    free(label_y);
    free(match_x);
    free(match_y);
    free(slack);
    free(slack_x);
    free(S_set);
    free(T_set);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
