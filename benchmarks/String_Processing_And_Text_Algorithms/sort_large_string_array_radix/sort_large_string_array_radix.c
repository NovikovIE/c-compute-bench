#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (Do Not Modify) ---
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
static char **string_array;
static long num_strings;
static int max_string_length;
static long final_result;

// --- Helper Functions for Computation ---

// Stable sort using counting sort on a specific character index.
static void counting_sort_by_char(int char_idx, char **temp_array) {
    int count[256] = {0}; // ASCII range

    // 1. Store count of occurrences of each character.
    for (long i = 0; i < num_strings; i++) {
        count[(unsigned char)string_array[i][char_idx]]++;
    }

    // 2. Modify count[i] to store the final position of this character in the output array.
    for (int i = 1; i < 256; i++) {
        count[i] += count[i - 1];
    }

    // 3. Build the output array. Iterate backwards to maintain stability.
    for (long i = num_strings - 1; i >= 0; i--) {
        temp_array[count[(unsigned char)string_array[i][char_idx]] - 1] = string_array[i];
        count[(unsigned char)string_array[i][char_idx]]--;
    }

    // 4. Copy the sorted pointers from the temporary array back to the main array.
    memcpy(string_array, temp_array, num_strings * sizeof(char *));
}

// LSD Radix sort for fixed-length strings.
static void radix_sort() {
    char **temp_array = (char **)malloc(num_strings * sizeof(char *));
    if (!temp_array) {
        fprintf(stderr, "FATAL: Failed to allocate temp array for sorting.\n");
        exit(1);
    }
    
    // Perform counting sort for every character from right to left (LSD).
    for (int i = max_string_length - 1; i >= 0; i--) {
        counting_sort_by_char(i, temp_array);
    }
    
    free(temp_array);
}

// --- Core Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_strings> <max_string_length> <seed>\n", argv[0]);
        exit(1);
    }
    
    num_strings = atol(argv[1]);
    max_string_length = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);
    
    if (num_strings <= 0 || max_string_length <= 0) {
        fprintf(stderr, "FATAL: num_strings and max_string_length must be positive.\n");
        exit(1);
    }
    
    string_array = (char **)malloc(num_strings * sizeof(char *));
    if (!string_array) {
        fprintf(stderr, "FATAL: Failed to allocate string array.\n");
        exit(1);
    }
    
    for (long i = 0; i < num_strings; ++i) {
        string_array[i] = (char *)malloc((max_string_length + 1) * sizeof(char));
        if (!string_array[i]) {
            fprintf(stderr, "FATAL: Failed to allocate string %ld.\n", i);
            // A more robust implementation would free previously allocated memory.
            exit(1);
        }
        for (int j = 0; j < max_string_length; ++j) {
            // Fill with printable ASCII characters ' ' to '~'
            string_array[i][j] = (char)(' ' + (mt_rand() % 95));
        }
        string_array[i][max_string_length] = '\0';
    }
}

void run_computation() {
    radix_sort();

    // Calculate a checksum on the middle string to prevent dead code elimination.
    final_result = 0;
    char* middle_string = string_array[num_strings / 2];
    for (int i = 0; i < max_string_length; ++i) {
        final_result += middle_string[i];
    }
}

void cleanup() {
    for (long i = 0; i < num_strings; i++) {
        free(string_array[i]);
    }
    free(string_array);
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
    printf("%ld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
