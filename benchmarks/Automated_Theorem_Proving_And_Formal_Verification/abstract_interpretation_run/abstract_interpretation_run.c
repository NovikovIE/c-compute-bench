#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
// Parameters
int num_program_variables;
int abstract_domain_size;
int num_loops;

// Data structures
// transition_tables[variable_index][current_value][predecessor_value] -> new_value
int*** transition_tables; 
// The state of abstract variables
int* program_state; 

// Final result to prevent dead code elimination
int final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_program_variables> <abstract_domain_size> <num_loops> <seed>\n", argv[0]);
        exit(1);
    }

    num_program_variables = atoi(argv[1]);
    abstract_domain_size = atoi(argv[2]);
    num_loops = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);

    mt_seed(seed);

    // Allocate initial program state vector
    program_state = (int*)malloc(num_program_variables * sizeof(int));
    if (!program_state) {
        fprintf(stderr, "Failed to allocate memory for program_state.\n");
        exit(1);
    }

    // Initialize program state with random abstract values
    for (int i = 0; i < num_program_variables; ++i) {
        program_state[i] = mt_rand() % abstract_domain_size;
    }

    // Allocate and populate the transition tables
    // These tables represent the abstract transformers for each variable
    transition_tables = (int***)malloc(num_program_variables * sizeof(int**));
    if (!transition_tables) {
        fprintf(stderr, "Failed to allocate memory for transition_tables.\n");
        exit(1);
    }

    for (int i = 0; i < num_program_variables; ++i) {
        transition_tables[i] = (int**)malloc(abstract_domain_size * sizeof(int*));
        if (!transition_tables[i]) {
            fprintf(stderr, "Failed to allocate memory for transition_tables[%d].\n", i);
            exit(1);
        }
        for (int j = 0; j < abstract_domain_size; ++j) {
            transition_tables[i][j] = (int*)malloc(abstract_domain_size * sizeof(int));
            if (!transition_tables[i][j]) {
                fprintf(stderr, "Failed to allocate memory for transition_tables[%d][%d].\n", i, j);
                exit(1);
            }
            for (int k = 0; k < abstract_domain_size; ++k) {
                // The new abstract value is a random function of the current and predecessor values
                transition_tables[i][j][k] = mt_rand() % abstract_domain_size;
            }
        }
    }
}

void run_computation() {
    // A temporary state for the next iteration
    int* next_state = (int*)malloc(num_program_variables * sizeof(int));
    if (!next_state) {
        fprintf(stderr, "Failed to allocate memory for next_state in run_computation.\n");
        exit(1);
    }

    // Simulate the fixpoint iteration process of abstract interpretation
    for (int iter = 0; iter < num_loops; ++iter) {
        // For each variable, compute its next abstract value based on the current state
        for (int i = 0; i < num_program_variables; ++i) {
            // Assume a cyclic dependency: var[i] depends on var[i-1] (and var[0] on var[N-1])
            int pred_idx = (i == 0) ? num_program_variables - 1 : i - 1;
            
            int current_val = program_state[i];
            int pred_val = program_state[pred_idx];

            // Apply the abstract transformer (lookup in the pre-computed table)
            next_state[i] = transition_tables[i][current_val][pred_val];
        }

        // Atomically update the state for the next iteration
        for (int i = 0; i < num_program_variables; ++i) {
            program_state[i] = next_state[i];
        }
    }

    free(next_state);

    // Compute a final result from the final state to prevent optimization
    int accumulator = 0;
    for (int i = 0; i < num_program_variables; ++i) {
        accumulator += program_state[i];
    }
    final_result = accumulator;
}

void cleanup() {
    for (int i = 0; i < num_program_variables; ++i) {
        for (int j = 0; j < abstract_domain_size; ++j) {
            free(transition_tables[i][j]);
        }
        free(transition_tables[i]);
    }
    free(transition_tables);
    free(program_state);
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
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
