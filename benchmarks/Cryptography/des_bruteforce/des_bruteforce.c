#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister PRNG (Provided) ---
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

// --- Minimal DES Implementation ---

// Permutation tables
static const uint8_t IP[] = {
    58, 50, 42, 34, 26, 18, 10, 2, 60, 52, 44, 36, 28, 20, 12, 4,
    62, 54, 46, 38, 30, 22, 14, 6, 64, 56, 48, 40, 32, 24, 16, 8,
    57, 49, 41, 33, 25, 17, 9,  1, 59, 51, 43, 35, 27, 19, 11, 3,
    61, 53, 45, 37, 29, 21, 13, 5, 63, 55, 47, 39, 31, 23, 15, 7
};

static const uint8_t FP[] = {
    40, 8, 48, 16, 56, 24, 64, 32, 39, 7, 47, 15, 55, 23, 63, 31,
    38, 6, 46, 14, 54, 22, 62, 30, 37, 5, 45, 13, 53, 21, 61, 29,
    36, 4, 44, 12, 52, 20, 60, 28, 35, 3, 43, 11, 51, 19, 59, 27,
    34, 2, 42, 10, 50, 18, 58, 26, 33, 1, 41,  9, 49, 17, 57, 25
};

static const uint8_t PC1[] = {
    57, 49, 41, 33, 25, 17, 9,  1, 58, 50, 42, 34, 26, 18,
    10, 2, 59, 51, 43, 35, 27, 19, 11, 3, 60, 52, 44, 36,
    63, 55, 47, 39, 31, 23, 15, 7, 62, 54, 46, 38, 30, 22,
    14, 6, 61, 53, 45, 37, 29, 21, 13, 5, 28, 20, 12, 4
};

static const uint8_t PC2[] = {
    14, 17, 11, 24, 1,  5,  3,  28, 15, 6,  21, 10,
    23, 19, 12, 4,  26, 8,  16, 7,  27, 20, 13, 2,
    41, 52, 31, 37, 47, 55, 30, 40, 51, 45, 33, 48,
    44, 49, 39, 56, 34, 53, 46, 42, 50, 36, 29, 32
};

static const uint8_t E[] = {
    32, 1,  2,  3,  4,  5,  4,  5,  6,  7,  8,  9,
    8,  9,  10, 11, 12, 13, 12, 13, 14, 15, 16, 17,
    16, 17, 18, 19, 20, 21, 20, 21, 22, 23, 24, 25,
    24, 25, 26, 27, 28, 29, 28, 29, 30, 31, 32, 1
};

static const uint8_t P[] = {
    16, 7,  20, 21, 29, 12, 28, 17, 1,  15, 23, 26, 5,  18, 31, 10,
    2,  8,  24, 14, 32, 27, 3,  9,  19, 13, 30, 6,  22, 11, 4,  25
};

static const uint8_t S_BOX[8][4][16] = {
    {{14,4,13,1,2,15,11,8,3,10,6,12,5,9,0,7},{0,15,7,4,14,2,13,1,10,6,12,11,9,5,3,8},{4,1,14,8,13,6,2,11,15,12,9,7,3,10,5,0},{15,12,8,2,4,9,1,7,5,11,3,14,10,0,6,13}},
    {{15,1,8,14,6,11,3,4,9,7,2,13,12,0,5,10},{3,13,4,7,15,2,8,14,12,0,1,10,6,9,11,5},{0,14,7,11,10,4,13,1,5,8,12,6,9,3,2,15},{13,8,10,1,3,15,4,2,11,6,7,12,0,5,14,9}},
    {{10,0,9,14,6,3,15,5,1,13,12,7,11,4,2,8},{13,7,0,9,3,4,6,10,2,8,5,14,12,11,15,1},{13,6,4,9,8,15,3,0,11,1,2,12,5,10,14,7},{1,10,13,0,6,9,8,7,4,15,14,3,11,5,2,12}},
    {{7,13,14,3,0,6,9,10,1,2,8,5,11,12,4,15},{13,8,11,5,6,15,0,3,4,7,2,12,1,10,14,9},{10,6,9,0,12,11,7,13,15,1,3,14,5,2,8,4},{3,15,0,6,10,1,13,8,9,4,5,11,12,7,2,14}},
    {{2,12,4,1,7,10,11,6,8,5,3,15,13,0,14,9},{14,11,2,12,4,7,13,1,5,0,15,10,3,9,8,6},{4,2,1,11,10,13,7,8,15,9,12,5,6,3,0,14},{11,8,12,7,1,14,2,13,6,15,0,9,10,4,5,3}},
    {{12,1,10,15,9,2,6,8,0,13,3,4,14,7,5,11},{10,15,4,2,7,12,9,5,6,1,13,14,0,11,3,8},{9,14,15,5,2,8,12,3,7,0,4,10,1,13,11,6},{4,3,2,12,9,5,15,10,11,14,1,7,6,0,8,13}},
    {{4,11,2,14,15,0,8,13,3,12,9,7,5,10,6,1},{13,0,11,7,4,9,1,10,14,3,5,12,2,15,8,6},{1,4,11,13,12,3,7,14,10,15,6,8,0,5,9,2},{6,11,13,8,1,4,10,7,9,5,0,15,14,2,3,12}},
    {{13,2,8,4,6,15,11,1,10,9,3,14,5,0,12,7},{1,15,13,8,10,3,7,4,12,5,6,11,0,14,9,2},{7,11,4,1,9,12,14,2,0,6,10,13,15,3,5,8},{2,1,14,7,4,10,8,13,15,12,9,0,3,5,6,11}}
};

static const uint8_t SHIFTS[] = {1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1};

// Helper to perform a permutation
static uint64_t permute(uint64_t input, const uint8_t *table, int size) {
    uint64_t output = 0;
    for (int i = 0; i < size; ++i) {
        if ((input >> (64 - table[i])) & 1) {
            output |= (1ULL << (63 - i));
        }
    }
    return output;
}

// Generate 16 48-bit subkeys from a 64-bit key
static void generate_subkeys(uint64_t key, uint64_t subkeys[16]) {
    uint64_t permuted_key = 0;
    for (int i = 0; i < 56; ++i) {
        if ((key >> (64 - PC1[i])) & 1) {
            permuted_key |= (1ULL << (55 - i));
        }
    }

    uint32_t c = (uint32_t)(permuted_key >> 28);
    uint32_t d = (uint32_t)(permuted_key & 0x0FFFFFFF);

    for (int i = 0; i < 16; ++i) {
        c = (c << SHIFTS[i]) | (c >> (28 - SHIFTS[i]));
        d = (d << SHIFTS[i]) | (d >> (28 - SHIFTS[i]));
        c &= 0x0FFFFFFF;
        d &= 0x0FFFFFFF;

        uint64_t cd = ((uint64_t)c << 28) | d;
        subkeys[i] = 0;
        for (int j = 0; j < 48; ++j) {
            if ((cd >> (56 - PC2[j])) & 1) {
                subkeys[i] |= (1ULL << (47 - j));
            }
        }
    }
}

// The DES f function
static uint32_t f_function(uint32_t r, uint64_t subkey) {
    uint64_t expanded_r = 0;
    for(int i = 0; i < 48; ++i){
        if((r >> (32 - E[i])) & 1){
            expanded_r |= (1ULL << (47 - i));
        }
    }

    uint64_t xored = expanded_r ^ subkey;
    uint32_t output = 0;

    for (int i = 0; i < 8; ++i) {
        uint8_t six_bits = (xored >> (48 - 6 * (i + 1))) & 0x3F;
        uint8_t row = ((six_bits & 0x20) >> 4) | (six_bits & 0x01);
        uint8_t col = (six_bits >> 1) & 0x0F;
        output |= (uint32_t)S_BOX[i][row][col] << (32 - 4 * (i + 1));
    }

    uint32_t p_permuted = 0;
    for(int i = 0; i < 32; ++i){
        if((output >> (32 - P[i])) & 1) {
            p_permuted |= (1U << (31 - i));
        }
    }
    return p_permuted;
}

// Encrypt a 64-bit block with pre-generated subkeys
static uint64_t des_encrypt_internal(uint64_t block, const uint64_t subkeys[16]) {
    block = permute(block, IP, 64) >> (64-64);

    uint32_t l = (uint32_t)(block >> 32);
    uint32_t r = (uint32_t)(block & 0xFFFFFFFF);

    for (int i = 0; i < 16; ++i) {
        uint32_t temp = r;
        r = l ^ f_function(r, subkeys[i]);
        l = temp;
    }

    uint64_t result = ((uint64_t)r << 32) | l;
    return permute(result, FP, 64) >> (64-64);
}

// Full DES encryption: key schedule + encryption
static uint64_t des_encrypt(uint64_t block, uint64_t key) {
    uint64_t subkeys[16];
    generate_subkeys(key, subkeys);
    return des_encrypt_internal(block, subkeys);
}

// --- Benchmark Globals and Functions ---

typedef struct {
    int key_space_bits_to_search;
    uint64_t plaintext;
    uint64_t target_ciphertext;
    uint64_t num_keys_to_check;
    uint64_t found_key;
} BenchmarkData;

static BenchmarkData* g_data = NULL;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <key_space_bits_to_search> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    g_data->key_space_bits_to_search = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);
    
    if (g_data->key_space_bits_to_search <= 0 || g_data->key_space_bits_to_search >= 56) {
        fprintf(stderr, "Key space bits must be between 1 and 55.\n");
        exit(1);
    }

    g_data->num_keys_to_check = 1ULL << g_data->key_space_bits_to_search;

    uint64_t secret_key = ((uint64_t)mt_rand() << 32) | mt_rand();
    g_data->plaintext = ((uint64_t)mt_rand() << 32) | mt_rand();
    
    // To ensure our search can theoretically succeed, place the key in the search space.
    // Note: this is for benchmark validation; a real attack wouldn't know this.
    uint64_t secret_key_in_search_space = mt_rand() % g_data->num_keys_to_check;

    g_data->target_ciphertext = des_encrypt(g_data->plaintext, secret_key_in_search_space);
    g_data->found_key = (uint64_t)-1; // Sentinel for 'not found'
}

void run_computation() {
    for (uint64_t key_candidate = 0; key_candidate < g_data->num_keys_to_check; ++key_candidate) {
        uint64_t encrypted = des_encrypt(g_data->plaintext, key_candidate);
        if (encrypted == g_data->target_ciphertext) {
            g_data->found_key = key_candidate;
            break;
        }
    }
}

void cleanup() {
    if (g_data) {
        free(g_data);
        g_data = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%llu\n", g_data->found_key);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
