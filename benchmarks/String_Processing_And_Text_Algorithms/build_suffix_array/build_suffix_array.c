#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- BENCHMARK DATA AND PARAMETERS ---
static int text_length;
static char* text;
static int* suffix_array;
static long long final_result; // Use long long for the accumulator

// --- HELPER FUNCTION FOR QSORT ---
// Comparison function for two suffixes, accessed via their start indices.
// The global 'text' pointer is used as context.
int compare_suffixes(const void *a, const void *b) {
    int index1 = *(const int*)a;
    int index2 = *(const int*)b;
    return strcmp(text + index1, text + index2);
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <text_length> <seed>\n", argv[0]);
        exit(1);
    }

    text_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (text_length <= 0) {
        fprintf(stderr, "FATAL: text_length must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory. Add +1 for the null terminator, crucial for strcmp.
    text = (char*)malloc((text_length + 1) * sizeof(char));
    suffix_array = (int*)malloc(text_length * sizeof(int));

    if (!text || !suffix_array) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a random string of lowercase characters.
    for (int i = 0; i < text_length; i++) {
        text[i] = 'a' + (mt_rand() % 26);
        suffix_array[i] = i;
    }
    text[text_length] = '\0'; // Null-terminate the string.

    final_result = 0;
}

void run_computation() {
    // Sort the array of suffix indices using the custom comparison function.
    qsort(suffix_array, text_length, sizeof(int), compare_suffixes);

    // Compute a checksum to prevent dead-code elimination and have a verifiable result.
    // A weighted sum ensures the order of the suffix array matters.
    for (int i = 0; i < text_length; i++) {
        final_result += (long long)suffix_array[i] * (i + 1);
    }
}

void cleanup() {
    free(text);
    free(suffix_array);
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout
    printf("%lld\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
