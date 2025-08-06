#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- BEGIN Mersenne Twister (Do Not Modify) ---
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
// --- END Mersenne Twister ---

// Benchmark parameters
static int NUM_VARIABLES;
static int NUM_CONSTRAINTS;

// Data structures for difference logic: x_v - x_u <= w
typedef struct {
    int u; // Source variable index
    int v; // Destination variable index
    int w; // Constant weight
} Constraint;

static Constraint* constraints;
static long long final_result = 0;

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_variables> <num_constraints> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_VARIABLES = atoi(argv[1]);
    NUM_CONSTRAINTS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_VARIABLES <= 0 || NUM_CONSTRAINTS <= 0) {
        fprintf(stderr, "Number of variables and constraints must be positive.\n");
        exit(1);
    }
    
    mt_seed(seed);

    constraints = (Constraint*)malloc(NUM_CONSTRAINTS * sizeof(Constraint));
    if (constraints == NULL) {
        fprintf(stderr, "Failed to allocate memory for constraints.\n");
        exit(1);
    }

    for (int i = 0; i < NUM_CONSTRAINTS; i++) {
        int u = mt_rand() % NUM_VARIABLES;
        int v;
        do {
            v = mt_rand() % NUM_VARIABLES;
        } while (u == v);

        constraints[i].u = u;
        constraints[i].v = v;
        // Generate weights in a small range to increase likelihood of interesting cycles
        constraints[i].w = (mt_rand() % 101) - 50;
    }
}

void run_computation() {
    // Bellman-Ford algorithm to solve the system of difference constraints.
    // The system is modeled as a graph where variables are vertices and
    // constraints are weighted, directed edges.
    // A solution exists if and only if there are no negative-weight cycles.

    long* dist = (long*)malloc(NUM_VARIABLES * sizeof(long));
    if (dist == NULL) {
        fprintf(stderr, "Failed to allocate memory for distance array.\n");
        exit(1);
    }

    // Initialize distances. Setting all to 0 is sufficient for cycle detection.
    for (int i = 0; i < NUM_VARIABLES; i++) {
        dist[i] = 0;
    }

    // Relax edges for V-1 iterations.
    for (int i = 0; i < NUM_VARIABLES - 1; i++) {
        for (int j = 0; j < NUM_CONSTRAINTS; j++) {
            Constraint c = constraints[j];
            if (dist[c.u] + c.w < dist[c.v]) {
                dist[c.v] = dist[c.u] + c.w;
            }
        }
    }

    // Check for negative-weight cycles. If a distance can still be reduced,
    // a negative cycle exists, and the system is unsatisfiable (UNSAT).
    int has_negative_cycle = 0;
    for (int j = 0; j < NUM_CONSTRAINTS; j++) {
        Constraint c = constraints[j];
        if (dist[c.u] + c.w < dist[c.v]) {
            has_negative_cycle = 1;
            break;
        }
    }

    if (has_negative_cycle) {
        final_result = -1; // UNSAT
    } else {
        // SAT: Calculate a checksum of distances as the final result.
        // This provides a non-trivial output and prevents dead code elimination.
        // Using XOR is fast and avoids overflow issues.
        long long checksum = 0;
        for (int i = 0; i < NUM_VARIABLES; i++) {
            checksum ^= dist[i];
        }
        final_result = checksum;
    }

    free(dist);
}

void cleanup() {
    free(constraints);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
