#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <limits.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of Mersenne Twister ---


// --- Benchmark Globals ---
int num_keys;
int *freq;
long long *prefix_sum;
long long **cost;
long long final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_keys> <seed>\n", argv[0]);
        exit(1);
    }

    num_keys = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (num_keys <= 0) {
        fprintf(stderr, "Error: num_keys must be a positive integer.\n");
        exit(1);
    }

    // Allocate memory
    freq = (int *)malloc(num_keys * sizeof(int));
    prefix_sum = (long long *)malloc((num_keys + 1) * sizeof(long long));
    cost = (long long **)malloc(num_keys * sizeof(long long *));
    for (int i = 0; i < num_keys; i++) {
        cost[i] = (long long *)malloc(num_keys * sizeof(long long));
    }

    if (!freq || !prefix_sum || !cost) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }
    for (int i = 0; i < num_keys; i++) {
        if (!cost[i]) {
            fprintf(stderr, "Memory allocation failed for cost table row.\n");
            exit(1);
        }
    }
    
    // Generate random frequencies for keys and calculate prefix sums
    prefix_sum[0] = 0;
    for (int i = 0; i < num_keys; i++) {
        freq[i] = (mt_rand() % 100) + 1; // Frequencies from 1 to 100
        prefix_sum[i + 1] = prefix_sum[i] + freq[i];
    }
}

void run_computation() {
    int n = num_keys;

    // Base case: For a single key i, the cost is its frequency
    for (int i = 0; i < n; i++) {
        cost[i][i] = freq[i];
    }

    // Fill the table for lengths from 2 to n
    for (int L = 2; L <= n; L++) {
        // For each possible starting key i
        for (int i = 0; i <= n - L; i++) {
            // Ending key j
            int j = i + L - 1;
            
            // Sum of frequencies for keys from i to j
            long long sum_freq = prefix_sum[j + 1] - prefix_sum[i];
            
            long long min_val = LLONG_MAX;

            // Find the optimal root r for the subtree from i to j
            for (int r = i; r <= j; r++) {
                // Cost is the cost of left subtree + right subtree
                // If a subtree does not exist, its cost is 0
                long long c = ((r > i) ? cost[i][r - 1] : 0) + 
                              ((r < j) ? cost[r + 1][j] : 0);
                if (c < min_val) {
                    min_val = c;
                }
            }
            
            // The cost of the tree is the sum of frequencies plus the minimum cost of subtrees
            cost[i][j] = sum_freq + min_val;
        }
    }

    // The final result is the cost of the OBST for all keys
    final_result = cost[0][n - 1];
}

void cleanup() {
    free(freq);
    free(prefix_sum);
    for (int i = 0; i < num_keys; i++) {
        free(cost[i]);
    }
    free(cost);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
