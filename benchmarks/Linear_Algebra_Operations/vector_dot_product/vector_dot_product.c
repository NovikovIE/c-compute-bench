#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Setup ---
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
unsigned long vector_length;
double *vector_a;
double *vector_b;
double dot_product_result;

// --- Utility Function ---
double random_double() {
    // Generate a random double between 0.0 and 1.0
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <vector_length> <seed>\n", argv[0]);
        exit(1);
    }

    vector_length = strtoul(argv[1], NULL, 10);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);

    if (vector_length == 0) {
        fprintf(stderr, "Error: vector_length must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    vector_a = (double*)malloc(vector_length * sizeof(double));
    vector_b = (double*)malloc(vector_length * sizeof(double));

    if (vector_a == NULL || vector_b == NULL) {
        fprintf(stderr, "Failed to allocate memory for vectors.\n");
        exit(1);
    }

    for (unsigned long i = 0; i < vector_length; i++) {
        vector_a[i] = random_double();
        vector_b[i] = random_double();
    }
}

void run_computation() {
    double sum = 0.0;
    for (unsigned long i = 0; i < vector_length; i++) {
        sum += vector_a[i] * vector_b[i];
    }
    dot_product_result = sum;
}

void cleanup() {
    free(vector_a);
    free(vector_b);
    vector_a = NULL;
    vector_b = NULL;
}

// --- Main Execution Logic ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", dot_product_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
