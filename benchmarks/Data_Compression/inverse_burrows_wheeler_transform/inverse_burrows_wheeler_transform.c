#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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
char* bwt_output_L;      // The last column of the BWT matrix
char* reconstructed_data; // Buffer for the reconstructed data
int* lf_map;             // The Last-to-First mapping vector
size_t block_size;       // Size of the data block in bytes
int primary_index;       // Original row index from BWT
long long validation_checksum; // To prevent dead code elimination

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s block_size_kb seed\n", argv[0]);
        exit(1);
    }

    int block_size_kb = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);

    block_size = (size_t)block_size_kb * 1024;
    if (block_size == 0) {
        fprintf(stderr, "FATAL: block_size_kb must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for benchmark data structures
    bwt_output_L = (char*)malloc(block_size * sizeof(char));
    reconstructed_data = (char*)malloc(block_size * sizeof(char));
    lf_map = (int*)malloc(block_size * sizeof(int));

    if (!bwt_output_L || !reconstructed_data || !lf_map) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate a random BWT output string (L)
    for (size_t i = 0; i < block_size; ++i) {
        bwt_output_L[i] = (char)(mt_rand() & 0xFF);
    }

    // Generate a random primary index
    primary_index = mt_rand() % block_size;

    validation_checksum = 0;
}

void run_computation() {
    // C-table for character counts and cumulative frequencies
    size_t C[256] = {0};

    // Step 1: Count character frequencies in L
    for (size_t i = 0; i < block_size; ++i) {
        C[(unsigned char)bwt_output_L[i]]++;
    }

    // Step 2: Convert counts to cumulative frequencies to find starting positions
    size_t sum = 0;
    for (int i = 0; i < 256; ++i) {
        size_t count = C[i];
        C[i] = sum;
        sum += count;
    }

    // Step 3: Build the LF-map. C[c] is updated to become C[c] + Occ(c, i)
    for (size_t i = 0; i < block_size; ++i) {
        unsigned char c = (unsigned char)bwt_output_L[i];
        lf_map[i] = C[c];
        C[c]++;
    }

    // Step 4: Reconstruct the original string by repeatedly applying the LF-map
    int current_index = primary_index;
    for (size_t i = 0; i < block_size; ++i) {
        // Reconstruct backwards from the last character
        size_t write_pos = block_size - 1 - i;
        reconstructed_data[write_pos] = bwt_output_L[current_index];
        current_index = lf_map[current_index];
    }

    // Step 5: Calculate a checksum to ensure the result is used
    long long checksum = 0;
    for (size_t i = 0; i < block_size; ++i) {
        checksum += reconstructed_data[i];
    }
    validation_checksum = checksum;
}

void cleanup() {
    free(bwt_output_L);
    free(reconstructed_data);
    free(lf_map);
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
    printf("%lld\n", validation_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
