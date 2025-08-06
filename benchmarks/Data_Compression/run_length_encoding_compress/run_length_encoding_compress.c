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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND GLOBALS ---
unsigned char *g_input_data = NULL;
unsigned char *g_compressed_data = NULL;
size_t g_input_size = 0;
size_t g_compressed_size = 0; // Final result

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_size_mb> <seed>\n", argv[0]);
        exit(1);
    }

    int input_size_mb = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    g_input_size = (size_t)input_size_mb * 1024 * 1024;
    if (g_input_size == 0) {
        fprintf(stderr, "FATAL: Input size cannot be zero.\n");
        exit(1);
    }

    // Allocate memory for input data
    g_input_data = (unsigned char*)malloc(g_input_size * sizeof(unsigned char));
    if (g_input_data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for input data.\n");
        exit(1);
    }

    // Allocate memory for compressed data. Worst case is 2x the original size.
    g_compressed_data = (unsigned char*)malloc(g_input_size * 2 * sizeof(unsigned char));
    if (g_compressed_data == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for compressed data.\n");
        free(g_input_data);
        exit(1);
    }

    // Generate pseudo-compressible random data by creating runs of characters
    unsigned char current_char = mt_rand() % 256;
    for (size_t i = 0; i < g_input_size; ++i) {
        // Change character with a small probability to create reasonably long runs
        if ((mt_rand() % 100) < 2) { // 2% chance to change character
            current_char = mt_rand() % 256;
        }
        g_input_data[i] = current_char;
    }
}

void run_computation() {
    size_t input_pos = 0;
    size_t output_pos = 0;
    
    while (input_pos < g_input_size) {
        unsigned char current_value = g_input_data[input_pos];
        unsigned char run_length = 1;
        
        input_pos++;
        
        // Find the length of the run, up to a maximum of 255
        while (input_pos < g_input_size && g_input_data[input_pos] == current_value && run_length < 255) {
            run_length++;
            input_pos++;
        }
        
        // Store the run length and the value
        g_compressed_data[output_pos++] = run_length;
        g_compressed_data[output_pos++] = current_value;
    }
    
    g_compressed_size = output_pos;
}

void cleanup() {
    if (g_input_data) {
        free(g_input_data);
        g_input_data = NULL;
    }
    if (g_compressed_data) {
        free(g_compressed_data);
        g_compressed_data = NULL;
    }
}

//--- MAIN ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print the result (compressed size) to stdout
    printf("%zu\n", g_compressed_size);

    cleanup();

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
