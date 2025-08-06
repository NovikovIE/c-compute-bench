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

// --- BENCHMARK IMPLEMENTATION ---

typedef struct {
    size_t input_size;
    unsigned char *input_data;
    size_t bwt_block_size;
    int num_blocks;
    unsigned long long final_result;
} BenchmarkData;

static BenchmarkData g_data;

// BWT helper data (used by qsort comparison function)
static const unsigned char *bwt_current_block;
static size_t bwt_current_block_size;

// Compare two rotations of the bwt_current_block
// a and b are pointers to integer indices into the block
static int compare_rotations(const void *a, const void *b) {
    int pos1 = *(const int*)a;
    int pos2 = *(const int*)b;
    for (size_t i = 0; i < bwt_current_block_size; ++i) {
        // Perform cyclic access into the block
        unsigned char c1 = bwt_current_block[(pos1 + i) % bwt_current_block_size];
        unsigned char c2 = bwt_current_block[(pos2 + i) % bwt_current_block_size];
        if (c1 != c2) {
            return c1 - c2;
        }
    }
    return 0; // Should not happen for distinct start indices
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <bwt_block_size_kb> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.input_size = (size_t)atoi(argv[1]) * 1024 * 1024;
    g_data.bwt_block_size = (size_t)atoi(argv[2]) * 1024;
    uint32_t seed = atoi(argv[3]);

    if (g_data.bwt_block_size == 0 || g_data.input_size == 0) {
        fprintf(stderr, "FATAL: Input size and block size must be greater than 0.\n");
        exit(1);
    }

    if (g_data.bwt_block_size > g_data.input_size) {
        g_data.bwt_block_size = g_data.input_size;
    }

    g_data.num_blocks = (g_data.input_size + g_data.bwt_block_size - 1) / g_data.bwt_block_size;

    mt_seed(seed);

    g_data.input_data = (unsigned char *)malloc(g_data.input_size);
    if (!g_data.input_data) {
        fprintf(stderr, "FATAL: Memory allocation failed for input_data.\n");
        exit(1);
    }

    // Generate semi-compressible data with runs of characters
    size_t i = 0;
    while (i < g_data.input_size) {
        unsigned char val = mt_rand() % 256;
        size_t run_len = 1 + (mt_rand() % 256);
        for (size_t j = 0; j < run_len && i < g_data.input_size; ++j, ++i) {
            g_data.input_data[i] = val;
        }
    }

    g_data.final_result = 0;
}

void run_computation() {
    g_data.final_result = 0;
    
    // Allocate temporary buffers for a single block processing
    int *rotations = (int *)malloc(g_data.bwt_block_size * sizeof(int));
    unsigned char *bwt_output = (unsigned char *)malloc(g_data.bwt_block_size);
    unsigned char *mtf_output = (unsigned char *)malloc(g_data.bwt_block_size);
    
    if (!rotations || !bwt_output || !mtf_output) {
        fprintf(stderr, "FATAL: Memory allocation failed for temp buffers.\n");
        exit(1);
    }

    for (int i = 0; i < g_data.num_blocks; ++i) {
        const unsigned char* current_block = g_data.input_data + (i * g_data.bwt_block_size);
        size_t current_block_size = (i == g_data.num_blocks - 1 && g_data.input_size % g_data.bwt_block_size != 0) ? 
                                       g_data.input_size % g_data.bwt_block_size : g_data.bwt_block_size;

        // --- BWT Stage ---
        bwt_current_block = current_block;
        bwt_current_block_size = current_block_size;
        for (size_t j = 0; j < current_block_size; ++j) {
            rotations[j] = j;
        }

        qsort(rotations, current_block_size, sizeof(int), compare_rotations);

        int original_index = -1;
        for (size_t j = 0; j < current_block_size; ++j) {
            int prev_char_idx = (rotations[j] + current_block_size - 1) % current_block_size;
            bwt_output[j] = current_block[prev_char_idx];
            if (rotations[j] == 0) {
                original_index = j;
            }
        }

        // --- MTF Stage ---
        unsigned char mtf_symbols[256];
        for (int k = 0; k < 256; ++k) mtf_symbols[k] = k;

        for (size_t j = 0; j < current_block_size; ++j) {
            unsigned char c = bwt_output[j];
            int pos = 0;
            for (pos = 0; pos < 256; ++pos) {
                if (mtf_symbols[pos] == c) break;
            }
            mtf_output[j] = (unsigned char)pos;
            if (pos > 0) {
                memmove(&mtf_symbols[1], &mtf_symbols[0], pos);
                mtf_symbols[0] = c;
            }
        }
        
        // --- Final Checksum Stage ---
        uint32_t checksum = original_index;
        for(size_t j = 0; j < current_block_size; j++) {
            checksum = (checksum * 31) + mtf_output[j];
        }
        g_data.final_result += checksum;
    }

    free(rotations);
    free(bwt_output);
    free(mtf_output);
}

void cleanup() {
    free(g_data.input_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%llu\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
