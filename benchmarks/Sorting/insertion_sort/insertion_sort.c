#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- Benchmark Globals ---
int* data;
int num_elements;
long long verification_result = 0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_elements> <initial_data_distribution> <seed>\n", argv[0]);
        fprintf(stderr, "Distributions: random, sorted, reverse_sorted, nearly_sorted\n");
        exit(1);
    }

    num_elements = atoi(argv[1]);
    const char* distribution = argv[2];
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_elements <= 0) {
        fprintf(stderr, "FATAL: num_elements must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    data = (int*)malloc(num_elements * sizeof(int));
    if (data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for data array.\n");
        exit(1);
    }

    if (strcmp(distribution, "random") == 0) {
        for (int i = 0; i < num_elements; ++i) {
            data[i] = mt_rand();
        }
    } else if (strcmp(distribution, "sorted") == 0) {
        for (int i = 0; i < num_elements; ++i) {
            data[i] = i;
        }
    } else if (strcmp(distribution, "reverse_sorted") == 0) {
        for (int i = 0; i < num_elements; ++i) {
            data[i] = num_elements - 1 - i;
        }
    } else if (strcmp(distribution, "nearly_sorted") == 0) {
        for (int i = 0; i < num_elements; ++i) {
            data[i] = i;
        }
        int swaps = num_elements / 100; // 1% of elements swapped
        if (swaps == 0) swaps = 1;
        for(int i = 0; i < swaps; ++i) {
            int idx1 = mt_rand() % num_elements;
            int idx2 = mt_rand() % num_elements;
            int temp = data[idx1];
            data[idx1] = data[idx2];
            data[idx2] = temp;
        }
    } else {
        fprintf(stderr, "FATAL: Unknown data distribution '%s'\n", distribution);
        free(data);
        exit(1);
    }
}

void run_computation() {
    int i, key, j;
    for (i = 1; i < num_elements; i++) {
        key = data[i];
        j = i - 1;
        while (j >= 0 && data[j] > key) {
            data[j + 1] = data[j];
            j = j - 1;
        }
        data[j + 1] = key;
    }

    // Checksum to prevent dead code elimination
    verification_result = data[0] + data[num_elements / 2] + data[num_elements - 1];
}

void cleanup() {
    free(data);
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
    printf("%lld\n", verification_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
