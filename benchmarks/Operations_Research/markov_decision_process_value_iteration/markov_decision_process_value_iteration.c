#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
// Do Not Modify
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
// --- End of MT19937 ---

// Benchmark parameters
int num_states;
int num_actions;
int horizon;
double discount_factor;

// Data structures
double *transition_probs; // Flattened 3D array: [state][action][next_state]
double *rewards;          // Flattened 2D array: [state][action]
double *values;           // Current value function for each state
double *next_values;      // Next value function being computed
double final_result = 0.0;

// Helper to generate a random double between 0.0 and 1.0
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <num_states> <num_actions> <horizon> <discount_factor> <seed>\n", argv[0]);
        exit(1);
    }
    
    num_states = atoi(argv[1]);
    num_actions = atoi(argv[2]);
    horizon = atoi(argv[3]);
    discount_factor = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    transition_probs = (double *)malloc((size_t)num_states * num_actions * num_states * sizeof(double));
    rewards = (double *)malloc((size_t)num_states * num_actions * sizeof(double));
    values = (double *)malloc((size_t)num_states * sizeof(double));
    next_values = (double *)malloc((size_t)num_states * sizeof(double));

    if (!transition_probs || !rewards || !values || !next_values) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Initialize data
    // Transition probabilities (must sum to 1 for each state-action pair)
    for (int s = 0; s < num_states; ++s) {
        for (int a = 0; a < num_actions; ++a) {
            double sum_prob = 0.0;
            size_t base_idx = (size_t)s * num_actions * num_states + (size_t)a * num_states; 
            for (int s_prime = 0; s_prime < num_states; ++s_prime) {
                double rand_val = rand_double();
                transition_probs[base_idx + s_prime] = rand_val;
                sum_prob += rand_val;
            }
            // Normalize probabilities
            for (int s_prime = 0; s_prime < num_states; ++s_prime) {
                transition_probs[base_idx + s_prime] /= sum_prob;
            }
        }
    }

    // Rewards (random values between 0 and 100)
    for (int s = 0; s < num_states; ++s) {
        for (int a = 0; a < num_actions; ++a) {
            rewards[(size_t)s * num_actions + a] = rand_double() * 100.0;
        }
    }

    // Initial value function (all zeros)
    for (int s = 0; s < num_states; ++s) {
        values[s] = 0.0;
    }
}

void run_computation() {
    double* current_v = values;
    double* next_v = next_values;

    for (int t = 0; t < horizon; ++t) {
        for (int s = 0; s < num_states; ++s) {
            double max_q_value = -INFINITY;
            for (int a = 0; a < num_actions; ++a) {
                size_t reward_idx = (size_t)s * num_actions + a;
                size_t trans_prob_base_idx = reward_idx * num_states;

                double expected_future_value = 0.0;
                for (int s_prime = 0; s_prime < num_states; ++s_prime) {
                    expected_future_value += transition_probs[trans_prob_base_idx + s_prime] * current_v[s_prime];
                }

                double q_value = rewards[reward_idx] + discount_factor * expected_future_value;
                if (q_value > max_q_value) {
                    max_q_value = q_value;
                }
            }
            next_v[s] = max_q_value;
        }
        
        // Swap pointers for the next iteration to avoid memcpy
        double *temp = current_v;
        current_v = next_v;
        next_v = temp;
    }

    // After the loop, current_v points to the final set of values.
    double checksum = 0.0;
    for (int s = 0; s < num_states; ++s) {
        checksum += current_v[s];
    }
    final_result = checksum;
}

void cleanup() {
    free(transition_probs);
    free(rewards);
    free(values);
    free(next_values);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
