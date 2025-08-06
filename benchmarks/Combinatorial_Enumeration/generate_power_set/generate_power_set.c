#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Globals ---
int num_elements;
int* element_set;
long long final_result;

// --- Benchmark Functions ---

/**
 * @brief Parses arguments, allocates memory, and generates the input data.
 * 
 * The power set of a set S is the set of all subsets of S. This function
 * sets up the initial set S with random integer values.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_elements> <seed>\n", argv[0]);
        exit(1);
    }

    num_elements = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (num_elements <= 0 || num_elements > 30) { // Keep num_elements in a reasonable range to avoid excessive memory/time
        fprintf(stderr, "FATAL: num_elements must be between 1 and 30.\n");
        exit(1);
    }

    mt_seed(seed);

    element_set = (int*)malloc(num_elements * sizeof(int));
    if (element_set == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_elements; ++i) {
        // Generate small integers to avoid overflow in the sum
        element_set[i] = mt_rand() % 100;
    }
}

/**
 * @brief Executes the core computation: generating the power set and summing elements.
 *
 * This function iterates through all 2^N possible subsets of the `element_set`.
 * It uses a bitmask approach where the i-th bit of a counter corresponds to the i-th
 * element of the set. To prevent the compiler from optimizing away the entire
 * computation, it calculates the sum of all elements across all generated subsets.
 * The final sum is stored in a global variable.
 */
void run_computation() {
    long long total_sum = 0;
    long long num_subsets = 1LL << num_elements; // 2^N

    // Iterate through all possible subsets from 0 to 2^N - 1
    for (long long i = 0; i < num_subsets; ++i) {
        long long current_subset_sum = 0;
        // Check which elements are in the current subset using bitmask 'i'
        for (int j = 0; j < num_elements; ++j) {
            if ((i >> j) & 1) {
                current_subset_sum += element_set[j];
            }
        }
        total_sum += current_subset_sum;
    }
    final_result = total_sum;
}

/**
 * @brief Frees all memory allocated in `setup_benchmark`.
 */
void cleanup() {
    if (element_set != NULL) {
        free(element_set);
        element_set = NULL;
    }
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

    // Print the final accumulated result to stdout
    printf("%lld\n", final_result);

    // Print the measured time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
