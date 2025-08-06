// BENCHMARK: Hashing and Checksums - rolling_hash_rabin_karp
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// Mersenne Twister (MT19937) non-cryptographic PRNG
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

// Benchmark specific defines
#define ALPHABET_SIZE 256
#define PRIME_MODULUS 1000000007 // A large prime number

// Global data structure to hold benchmark data
typedef struct {
    size_t window_size;
    size_t text_length;
    char* text;
    uint32_t final_result;
} BenchmarkData;

static BenchmarkData g_data;

// Separated setup, computation, and cleanup functions

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <window_size> <text_length> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.window_size = strtoul(argv[1], NULL, 10);
    g_data.text_length = strtoul(argv[2], NULL, 10);
    uint32_t seed = strtoul(argv[3], NULL, 10);

    if (g_data.window_size == 0 || g_data.window_size > g_data.text_length) {
        fprintf(stderr, "FATAL: Invalid arguments. window_size must be > 0 and <= text_length.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.text = (char*)malloc(g_data.text_length * sizeof(char));
    if (!g_data.text) {
        fprintf(stderr, "FATAL: Memory allocation failed for text.\n");
        exit(1);
    }

    for (size_t i = 0; i < g_data.text_length; i++) {
        g_data.text[i] = (char)(mt_rand() % 256);
    }
    
    g_data.final_result = 0;
}

void run_computation() {
    size_t w_size = g_data.window_size;
    size_t t_len = g_data.text_length;
    char *text = g_data.text;

    long long h = 1;
    // Calculate h = (ALPHABET_SIZE^(w_size-1)) % PRIME_MODULUS
    for (size_t i = 0; i < w_size - 1; i++) {
        h = (h * ALPHABET_SIZE) % PRIME_MODULUS;
    }

    long long current_hash = 0;
    // Calculate hash for the first window
    for (size_t i = 0; i < w_size; i++) {
        current_hash = (current_hash * ALPHABET_SIZE + (unsigned char)text[i]) % PRIME_MODULUS;
    }

    uint64_t hash_accumulator = current_hash;

    // Roll the hash for the rest of the text
    for (size_t i = 1; i <= t_len - w_size; i++) {
        unsigned char prev_char = (unsigned char)text[i - 1];
        unsigned char next_char = (unsigned char)text[i + w_size - 1];
        
        // Rolling hash update formula: hash = (d * (hash - text[old] * h) + text[new])
        long long new_hash_calc = ALPHABET_SIZE * (current_hash - prev_char * h) + next_char;

        // The remainder operator % in C can yield a negative result, so we ensure it's positive
        current_hash = (new_hash_calc % PRIME_MODULUS + PRIME_MODULUS) % PRIME_MODULUS;
        
        hash_accumulator ^= current_hash;
    }

    g_data.final_result = (uint32_t)hash_accumulator;
}

void cleanup() {
    free(g_data.text);
    g_data.text = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%u\n", g_data.final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
