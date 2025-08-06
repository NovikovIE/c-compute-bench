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

// Input data: two arrays of strings representing lines of text.
char** file1_data = NULL;
char** file2_data = NULL;
int num_file1_lines;
int num_file2_lines;

// To avoid freeing the same string multiple times, we use a single vocabulary.
char** vocabulary = NULL;
const int VOCAB_SIZE = 1000;

// Output result
int final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s file1_lines file2_lines seed\n", argv[0]);
        exit(1);
    }

    num_file1_lines = atoi(argv[1]);
    num_file2_lines = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    // 1. Create a vocabulary of unique random strings.
    vocabulary = (char**)malloc(VOCAB_SIZE * sizeof(char*));
    if (!vocabulary) {
        perror("malloc failed for vocabulary");
        exit(1);
    }
    for (int i = 0; i < VOCAB_SIZE; i++) {
        int len = 20 + (mt_rand() % 60); // Line length between 20 and 79
        vocabulary[i] = (char*)malloc(len + 1);
        if (!vocabulary[i]) {
            perror("malloc failed for vocabulary string");
            exit(1);
        }
        for (int j = 0; j < len; j++) {
            vocabulary[i][j] = 'a' + (mt_rand() % 26);
        }
        vocabulary[i][len] = '\0';
    }

    // 2. Create file1 by picking random lines from the vocabulary.
    file1_data = (char**)malloc(num_file1_lines * sizeof(char*));
    if (!file1_data) {
        perror("malloc failed for file1_data");
        exit(1);
    }
    for (int i = 0; i < num_file1_lines; i++) {
        file1_data[i] = vocabulary[mt_rand() % VOCAB_SIZE];
    }

    // 3. Create file2 similarly.
    file2_data = (char**)malloc(num_file2_lines * sizeof(char*));
    if (!file2_data) {
        perror("malloc failed for file2_data");
        exit(1);
    }
    for (int i = 0; i < num_file2_lines; i++) {
        file2_data[i] = vocabulary[mt_rand() % VOCAB_SIZE];
    }
}

void run_computation() {
    const int N = num_file1_lines;
    const int M = num_file2_lines;
    const int MAX_D = N + M;
    final_result = -1;

    // The V array stores the furthest reaching x-coordinate for each diagonal k.
    // We use an offset because k can be negative.
    int* V = (int*)malloc((2 * MAX_D + 1) * sizeof(int));
    if (!V) {
      perror("Failed to allocate V array");
      exit(1);
    }
    int offset = MAX_D;

    // Initialize V for the starting point.
    // V[k] = x, means we can reach (x, y) where y = x - k.
    // We start at (0, -1) on diagonal k=1, which after the first step becomes (0,0).
    V[offset + 1] = 0;

    for (int d = 0; d <= MAX_D; d++) {
        for (int k = -d; k <= d; k += 2) {
            int x;
            // Decide whether to come from diagonal k-1 (insertion) or k+1 (deletion).
            // Pick the one that gets us further in the x direction.
            if (k == -d || (k != d && V[offset + k - 1] < V[offset + k + 1])) {
                x = V[offset + k + 1]; // Move down.
            } else {
                x = V[offset + k - 1] + 1; // Move right.
            }
            
            int y = x - k;

            // Follow the diagonal as long as lines match.
            while (x < N && y < M && strcmp(file1_data[x], file2_data[y]) == 0) {
                x++;
                y++;
            }

            V[offset + k] = x;

            // Check if we've reached the end.
            if (x >= N && y >= M) {
                final_result = d;
                free(V);
                return;
            }
        }
    }
    
    // Should not be reached in a typical diff scenario
    final_result = MAX_D;
    free(V);
}

void cleanup() {
    // Free the vocabulary, which holds all unique strings.
    if (vocabulary) {
        for (int i = 0; i < VOCAB_SIZE; i++) {
            free(vocabulary[i]);
        }
        free(vocabulary);
    }
    
    // Free the file data arrays (but not their contents, which point to the vocabulary).
    free(file1_data);
    free(file2_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000000.0;

    // Print the final result (shortest edit distance) to stdout
    printf("%d\n", final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
