#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// --- BENCHMARK DATA AND GLOBALS ---
typedef struct {
    int num_crossings;
    int *crossing_types; // Represents the signs (+1 or -1) of crossings in a knot diagram
    long long result;    // Stores the final computed value
} BenchmarkData;

BenchmarkData *data = NULL;

// --- BENCHMARK-SPECIFIC LOGIC ---

// Recursive function to calculate a value based on the knot structure.
// This simulates the state-sum model for the Jones polynomial by exploring a binary tree.
// The structure of the tree and the operations performed are determined by the crossing types.
long long calculate_jones_recursive(int index) {
    if (index >= data->num_crossings) {
        // Base case: All crossings have been resolved. In the state-sum model,
        // this corresponds to a collection of unknots (a base state).
        return 1;
    }

    // Each crossing has two possible "smoothings" (resolutions), leading to two recursive calls.
    // This creates a binary recursion tree of depth 'num_crossings'.
    long long path_a = calculate_jones_recursive(index + 1);
    long long path_b = calculate_jones_recursive(index + 1);

    // The type of the crossing (pre-generated randomly) determines how the results
    // from the subproblems are combined. This is a simplification of how coefficients
    // (like t or t^-1) would be applied in the actual Jones polynomial calculation.
    if (data->crossing_types[index] > 0) {
        // A positive crossing might correspond to one form of the skein relation
        return path_a + path_b;
    } else {
        // A negative crossing corresponds to another form
        return path_a - path_b;
    }
}

// --- BENCHMARK CORE FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_crossings> <seed>\n", argv[0]);
        exit(1);
    }

    data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!data) {
        perror("Failed to allocate memory for benchmark data");
        exit(1);
    }

    data->num_crossings = atoi(argv[1]);
    uint32_t seed = atol(argv[2]);

    if (data->num_crossings <= 0) {
        fprintf(stderr, "FATAL: num_crossings must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    data->crossing_types = (int*)malloc(data->num_crossings * sizeof(int));
    if (!data->crossing_types) {
        perror("Failed to allocate memory for crossing types");
        free(data);
        exit(1);
    }

    // Generate a random sequence of crossing types for the knot diagram.
    for (int i = 0; i < data->num_crossings; i++) {
        // Assign +1 for a positive crossing, -1 for a negative crossing.
        data->crossing_types[i] = (mt_rand() % 2 == 0) ? 1 : -1;
    }
    
    data->result = 0;
}

void run_computation() {
    // The complexity is exponential, O(2^N), where N is num_crossings.
    // This function call kicks off the recursion.
    data->result = calculate_jones_recursive(0);
}

void cleanup() {
    if (data) {
        if (data->crossing_types) {
            free(data->crossing_types);
        }
        free(data);
        data = NULL;
    }
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start_time, end_time;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start_time);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    
    double time_taken = (end_time.tv_sec - start_time.tv_sec) + 
                        (end_time.tv_nsec - start_time.tv_nsec) / 1e9;

    // Print the final result to stdout to be captured.
    // This value is an integer accumulator that prevents dead code elimination.
    printf("%lld\n", data->result);

    cleanup();

    // Print the time taken to stderr, as required.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
