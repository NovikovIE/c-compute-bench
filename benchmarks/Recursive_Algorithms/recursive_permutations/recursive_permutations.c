#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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

// --- Benchmark Globals ---
int g_num_elements;
int* g_elements;
int* g_used;
int* g_current_permutation;
unsigned long long g_accumulator = 0;

// Recursive function prototype
void generate_permutations(int depth);

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_elements> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_elements = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);
    
    if (g_num_elements <= 0 || g_num_elements > 15) { // N! grows too fast
         fprintf(stderr, "FATAL: num_elements must be between 1 and 15.\n");
         exit(1);
    }

    mt_seed(seed); // Seed the generator as required even if not used for data generation

    g_elements = (int*)malloc(g_num_elements * sizeof(int));
    g_used = (int*)calloc(g_num_elements, sizeof(int));
    g_current_permutation = (int*)malloc(g_num_elements * sizeof(int));

    if (!g_elements || !g_used || !g_current_permutation) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < g_num_elements; ++i) {
        g_elements[i] = i;
    }
}

void cleanup() {
    free(g_elements);
    free(g_used);
    free(g_current_permutation);
}

void generate_permutations(int depth) {
    if (depth == g_num_elements) {
        // Base case: process the generated permutation
        for (int i = 0; i < g_num_elements; ++i) {
            // Arbitrary computation to use the result and prevent dead code elimination
            g_accumulator += (unsigned long long)g_current_permutation[i] * (i + 1);
            g_accumulator ^= ((unsigned long long)depth << 3);
        }
        return;
    }

    for (int i = 0; i < g_num_elements; ++i) {
        if (!g_used[i]) {
            g_used[i] = 1;
            g_current_permutation[depth] = g_elements[i];

            generate_permutations(depth + 1);

            g_used[i] = 0; // Backtrack
        }
    }
}

void run_computation() {
    generate_permutations(0);
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

    // Print result to stdout
    printf("%llu\n", g_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
