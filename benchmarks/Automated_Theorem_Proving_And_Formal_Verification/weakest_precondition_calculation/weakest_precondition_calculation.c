#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator --- DO NOT MODIFY ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Data Structures ---

typedef enum {
    STMT_ASSIGN, // Represents a simple variable assignment
    STMT_LOOP    // Represents a loop construct
} StatementType;

typedef struct {
    StatementType type;
    uint32_t data1; // Coefficients or data used in transformations
    uint32_t data2;
} Statement;

// Global struct to hold all benchmark data and parameters.
// This avoids passing large structs or multiple pointers between functions.
struct {
    int num_statements;
    int loop_nesting_depth;
    int assertion_complexity;
    uint32_t seed;
    
    Statement *program;
    uint32_t *initial_assertion; // The base postcondition
    uint32_t *working_assertion; // The assertion being transformed

    // Volatile to ensure the result is computed and not optimized away
    volatile long long final_result;
} g_data;

// --- Benchmark Functions ---

// Parses arguments, allocates memory, and generates pseudo-random input data.
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_statements loop_nesting_depth assertion_complexity seed\n", argv[0]);
        exit(1);
    }

    g_data.num_statements = atoi(argv[1]);
    g_data.loop_nesting_depth = atoi(argv[2]);
    g_data.assertion_complexity = atoi(argv[3]);
    g_data.seed = (uint32_t)strtoul(argv[4], NULL, 10);

    if (g_data.num_statements <= 0 || g_data.loop_nesting_depth <= 0 || g_data.assertion_complexity <= 0) {
        fprintf(stderr, "FATAL: All parameters must be positive integers.\n");
        exit(1);
    }
    
    mt_seed(g_data.seed);

    g_data.program = (Statement *)malloc(g_data.num_statements * sizeof(Statement));
    g_data.initial_assertion = (uint32_t *)malloc(g_data.assertion_complexity * sizeof(uint32_t));
    g_data.working_assertion = (uint32_t *)malloc(g_data.assertion_complexity * sizeof(uint32_t));

    if (!g_data.program || !g_data.initial_assertion || !g_data.working_assertion) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a program as a sequence of statements
    for (int i = 0; i < g_data.num_statements; i++) {
        // Approx. 5% of statements are complex loops, the rest are simple assignments
        if (mt_rand() % 100 < 5) {
            g_data.program[i].type = STMT_LOOP;
        } else {
            g_data.program[i].type = STMT_ASSIGN;
        }
        g_data.program[i].data1 = mt_rand();
        g_data.program[i].data2 = mt_rand();
    }

    // Generate the initial logical assertion (postcondition)
    for (int i = 0; i < g_data.assertion_complexity; i++) {
        g_data.initial_assertion[i] = mt_rand();
    }
}

// Executes the core computation of the benchmark.
void run_computation() {
    // A large prime for modular arithmetic to keep values within uint32_t range
    const uint32_t PRIME = 2147483647; 

    // Start with the postcondition
    for (int i = 0; i < g_data.assertion_complexity; i++) {
        g_data.working_assertion[i] = g_data.initial_assertion[i];
    }

    // Calculate the weakest precondition by processing the program backwards.
    for (int i = g_data.num_statements - 1; i >= 0; i--) {
        Statement *stmt = &g_data.program[i];

        if (stmt->type == STMT_ASSIGN) {
            // Simulate wp(assignment, Q): transform the assertion based on assignment data.
            for (int j = 0; j < g_data.assertion_complexity; j++) {
                g_data.working_assertion[j] = ((uint64_t)g_data.working_assertion[j] * stmt->data1 + stmt->data2) % PRIME;
            }
        } else { // STMT_LOOP
            // Simulate wp(loop, Q): an iterative process to find a fixed point (loop invariant).
            // The number of iterations is controlled by loop_nesting_depth.
            int iterations = g_data.loop_nesting_depth * 5; // Multiplier simulates complexity of finding the invariant
            for (int k = 0; k < iterations; k++) {
                // Each iteration is a complex transformation on the entire assertion.
                for (int j = 0; j < g_data.assertion_complexity; j++) {
                    // The dependency on the next element helps prevent auto-vectorization, making computation more predictable.
                    uint32_t next_elem = g_data.working_assertion[(j + 1) % g_data.assertion_complexity];
                    g_data.working_assertion[j] = (((uint64_t)g_data.working_assertion[j] ^ stmt->data1) + (next_elem ^ stmt->data2)) % PRIME;
                }
            }
        }
    }
    
    // Accumulate the final computed precondition into a single value.
    // This prevents the compiler from optimizing away the entire computation.
    long long sum = 0;
    for (int i = 0; i < g_data.assertion_complexity; i++) {
        sum += g_data.working_assertion[i];
    }
    g_data.final_result = sum;
}

// Frees all memory allocated in setup_benchmark.
void cleanup() {
    free(g_data.program);
    free(g_data.initial_assertion);
    free(g_data.working_assertion);
    g_data.program = NULL;
    g_data.initial_assertion = NULL;
    g_data.working_assertion = NULL;
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%lld\n", g_data.final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
