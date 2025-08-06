#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator ---
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

// --- DES Implementation Details ---
// Permutation tables, S-boxes, etc. (Standard FIPS 46-3 constants)

static const uint8_t IP[] = {
    58, 50, 42, 34, 26, 18, 10, 2,
    60, 52, 44, 36, 28, 20, 12, 4,
    62, 54, 46, 38, 30, 22, 14, 6,
    64, 56, 48, 40, 32, 24, 16, 8,
    57, 49, 41, 33, 25, 17, 9,  1,
    59, 51, 43, 35, 27, 19, 11, 3,
    61, 53, 45, 37, 29, 21, 13, 5,
    63, 55, 47, 39, 31, 23, 15, 7
};

static const uint8_t FP[] = {
    40, 8, 48, 16, 56, 24, 64, 32,
    39, 7, 47, 15, 55, 23, 63, 31,
    38, 6, 46, 14, 54, 22, 62, 30,
    37, 5, 45, 13, 53, 21, 61, 29,
    36, 4, 44, 12, 52, 20, 60, 28,
    35, 3, 43, 11, 51, 19, 59, 27,
    34, 2, 42, 10, 50, 18, 58, 26,
    33, 1, 41,  9, 49, 17, 57, 25
};

static const uint8_t E[] = {
    32, 1,  2,  3,  4,  5,
    4,  5,  6,  7,  8,  9,
    8,  9,  10, 11, 12, 13,
    12, 13, 14, 15, 16, 17,
    16, 17, 18, 19, 20, 21,
    20, 21, 22, 23, 24, 25,
    24, 25, 26, 27, 28, 29,
    28, 29, 30, 31, 32, 1
};

static const uint8_t P[] = {
    16, 7,  20, 21, 29, 12, 28, 17,
    1,  15, 23, 26, 5,  18, 31, 10,
    2,  8,  24, 14, 32, 27, 3,  9,
    19, 13, 30, 6,  22, 11, 4,  25
};

static const uint8_t S[8][4][16] = {
    {{14,4,13,1,2,15,11,8,3,10,6,12,5,9,0,7},{0,15,7,4,14,2,13,1,10,6,12,11,9,5,3,8},{4,1,14,8,13,6,2,11,15,12,9,7,3,10,5,0},{15,12,8,2,4,9,1,7,5,11,3,14,10,0,6,13}},
    {{15,1,8,14,6,11,3,4,9,7,2,13,12,0,5,10},{3,13,4,7,15,2,8,14,12,0,1,10,6,9,11,5},{0,14,7,11,10,4,13,1,5,8,12,6,9,3,2,15},{13,8,10,1,3,15,4,2,11,6,7,12,0,5,14,9}},
    {{10,0,9,14,6,3,15,5,1,13,12,7,11,4,2,8},{13,7,0,9,3,4,6,10,2,8,5,14,12,11,15,1},{13,6,4,9,8,15,3,0,11,1,2,12,5,10,14,7},{1,10,13,0,6,9,8,7,4,15,14,3,11,5,2,12}},
    {{7,13,14,3,0,6,9,10,1,2,8,5,11,12,4,15},{13,8,11,5,6,15,0,3,4,7,2,12,1,10,14,9},{10,6,9,0,12,11,7,13,15,1,3,14,5,2,8,4},{3,15,0,6,10,1,13,8,9,4,5,11,12,7,2,14}},
    {{2,12,4,1,7,10,11,6,8,5,3,15,13,0,14,9},{14,11,2,12,4,7,13,1,5,0,15,10,3,9,8,6},{4,2,1,11,10,13,7,8,15,9,12,5,6,3,0,14},{11,8,12,7,1,14,2,13,6,15,0,9,10,4,5,3}},
    {{12,1,10,15,9,2,6,8,0,13,3,4,14,7,5,11},{10,15,4,2,7,12,9,5,6,1,13,14,0,11,3,8},{9,14,15,5,2,8,12,3,7,0,4,10,1,13,11,6},{4,3,2,12,9,5,15,10,11,14,1,7,6,0,8,13}},
    {{4,11,2,14,15,0,8,13,3,12,9,7,5,10,6,1},{13,0,11,7,4,9,1,10,14,3,5,12,2,15,8,6},{1,4,11,13,12,3,7,14,10,15,6,8,0,5,9,2},{6,11,13,8,1,4,10,7,9,5,0,15,14,2,3,12}},
    {{13,2,8,4,6,15,11,1,10,9,3,14,5,0,12,7},{1,15,13,8,10,3,7,4,12,5,6,11,0,14,9,2},{7,11,4,1,9,12,14,2,0,6,10,13,15,3,5,8},{2,1,14,7,4,10,8,13,15,12,9,0,3,5,6,11}}
};

static const uint8_t PC1[] = {
    57, 49, 41, 33, 25, 17, 9,  1, 58, 50, 42, 34, 26, 18,
    10, 2,  59, 51, 43, 35, 27, 19, 11, 3,  60, 52, 44, 36,
    63, 55, 47, 39, 31, 23, 15, 7,  62, 54, 46, 38, 30, 22,
    14, 6,  61, 53, 45, 37, 29, 21, 13, 5,  28, 20, 12, 4
};

static const uint8_t PC2[] = {
    14, 17, 11, 24, 1,  5,  3,  28, 15, 6,  21, 10,
    23, 19, 12, 4,  26, 8,  16, 7,  27, 20, 13, 2,
    41, 52, 31, 37, 47, 55, 30, 40, 51, 45, 33, 48,
    44, 49, 39, 56, 34, 53, 46, 42, 50, 36, 29, 32
};

static const uint8_t key_shifts[] = {1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1};

// --- Benchmark Globals ---
static int KEY_SPACE_BITS_TO_SEARCH;
static uint8_t TARGET_CIPHERTEXT[8];
static const uint8_t KNOWN_PLAINTEXT[8] = {0x12, 0x34, 0x56, 0x78, 0x90, 0xAB, 0xCD, 0xEF};
static uint64_t *found_count; // Using heap allocation as per requirements
static uint64_t num_keys_to_check;

// --- DES Helper Functions ---

static void permute(uint8_t *out, const uint8_t *in, const uint8_t *table, int n) {
    for (int i = 0; i < n; i++) {
        int pos = table[i] - 1;
        out[i / 8] |= ((in[pos / 8] >> (7 - (pos % 8))) & 1) << (7 - (i % 8));
    }
}

static void bit_slice(uint8_t *out, const uint8_t *in, int start, int len) {
    for(int i = 0; i < len; ++i) {
        int pos = start + i;
        if ((in[pos / 8] >> (7 - (pos % 8))) & 1) {
            out[i/8] |= (1 << (7 - (i % 8)));
        }
    }
}

static void left_shift_28(uint8_t *data, int shifts) {
    uint32_t val = (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | (data[3] & 0xF0);
    val >>= 4;
    val = (val << shifts) | (val >> (28 - shifts));
    val <<= 4;
    data[0] = (val >> 24) & 0xFF;
    data[1] = (val >> 16) & 0xFF;
    data[2] = (val >> 8) & 0xFF;
    data[3] = val & 0xF0;
}

static void generate_subkeys(const uint8_t *key64, uint8_t subkeys[16][6]) {
    uint8_t key56[7] = {0};
    uint8_t C[4] = {0}, D[4] = {0};
    
    permute(key56, key64, PC1, 56);
    bit_slice(C, key56, 0, 28);
    bit_slice(D, key56, 28, 28);

    for (int i = 0; i < 16; i++) {
        left_shift_28(C, key_shifts[i]);
        left_shift_28(D, key_shifts[i]);

        uint8_t CD[7] = {0};
        memcpy(CD, C, 3);
        CD[3] = C[3] | (D[0] >> 4); 
        memcpy(CD+4, D+1, 3);

        memset(subkeys[i], 0, 6);
        permute(subkeys[i], CD, PC2, 48);
    }
}

static void des_crypt(uint8_t *out, const uint8_t *in, const uint8_t subkeys[16][6], int is_decrypt) {
    memset(out, 0, 8); // Fix: `permute` requires zeroed output buffer.
    uint8_t permuted_block[8] = {0}, L[4], R[4], temp[4], f_res[4] = {0}; // Fix: `f_res` must be initialized.
    
    permute(permuted_block, in, IP, 64);
    memcpy(L, permuted_block, 4);
    memcpy(R, permuted_block + 4, 4);

    for (int i = 0; i < 16; i++) {
        memcpy(temp, R, 4);
        // Feistel function
        // Fix: `xor_res` was read out-of-bounds. Padded to 8 bytes and zero-initialized.
        uint8_t expanded_R[6] = {0}, xor_res[8] = {0}, sbox_out[4] = {0};
        permute(expanded_R, R, E, 48);
        
        const uint8_t *subkey = subkeys[is_decrypt ? (15 - i) : i];
        for(int j=0; j<6; j++) xor_res[j] = expanded_R[j] ^ subkey[j];

        for(int j=0; j<8; j++) {
            uint8_t b = xor_res[j*6 / 8];
            uint8_t val6 = (b << (j*6 % 8)) | (xor_res[j*6/8+1] >> (8 - (j*6%8)));
            val6 &= 0xFC; // get 6 bits
            
            uint8_t row = ((val6 & 0x80) >> 6) | (val6 & 0x04) >> 2;
            uint8_t col = (val6 & 0x78) >> 3;
            uint8_t s_val = S[j][row][col];
            sbox_out[j/2] |= (j%2 == 0) ? (s_val << 4) : s_val;
        }

        permute(f_res, sbox_out, P, 32);
        for(int j=0; j<4; j++) R[j] = L[j] ^ f_res[j];
        memcpy(L, temp, 4);
    }

    uint8_t final_block[8];
    memcpy(final_block, R, 4);
    memcpy(final_block + 4, L, 4);
    permute(out, final_block, FP, 64);
}

static void triple_des_decrypt_block(const uint8_t *in, uint8_t *out, 
                                       const uint8_t *key1, const uint8_t *key2, const uint8_t *key3) {
    uint8_t subkeys1[16][6], subkeys2[16][6], subkeys3[16][6];
    uint8_t temp1[8], temp2[8];

    generate_subkeys(key1, subkeys1);
    generate_subkeys(key2, subkeys2);
    generate_subkeys(key3, subkeys3);

    des_crypt(temp1, in,    subkeys3, 1); // Decrypt with K3
    des_crypt(temp2, temp1, subkeys2, 0); // Encrypt with K2
    des_crypt(out,   temp2, subkeys1, 1); // Decrypt with K1
}


// --- Benchmark Functions ---

static int hex_char_to_int(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    return -1;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <key_space_bits_to_search> <ciphertext_hex_8_bytes> <seed>\n", argv[0]);
        exit(1);
    }

    KEY_SPACE_BITS_TO_SEARCH = atoi(argv[1]);
    if (KEY_SPACE_BITS_TO_SEARCH <= 0 || KEY_SPACE_BITS_TO_SEARCH > 63) {
        fprintf(stderr, "FATAL: key_space_bits_to_search must be between 1 and 63.\n");
        exit(1);
    }
    
    if (strlen(argv[2]) != 16) {
        fprintf(stderr, "FATAL: ciphertext_hex must be 16 hex characters (8 bytes).\n");
        exit(1);
    }

    for (size_t i = 0; i < 8; i++) {
        int high = hex_char_to_int(argv[2][i*2]);
        int low = hex_char_to_int(argv[2][i*2 + 1]);
        if (high == -1 || low == -1) {
            fprintf(stderr, "FATAL: Invalid hex character in ciphertext.\n");
            exit(1);
        }
        TARGET_CIPHERTEXT[i] = (high << 4) | low;
    }

    uint32_t seed = atoi(argv[3]);
    mt_seed(seed);

    num_keys_to_check = 1ULL << KEY_SPACE_BITS_TO_SEARCH;

    found_count = (uint64_t *)malloc(sizeof(uint64_t));
    if (found_count == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
    *found_count = 0;
}

void run_computation() {
    uint8_t current_k1[8];
    uint8_t current_k2[8];
    uint8_t decrypted_block[8];

    // K2 is fixed for this benchmark run. In a real scenario, it might be part of the search.
    // Here we generate it once to add to setup cost but it's not used by MT rand since seed is fixed.
    uint64_t k2_val = ((uint64_t)mt_rand() << 32) | mt_rand();
    memcpy(current_k2, &k2_val, 8);

    // We search the space of K1. K3=K1 for keying option 2.
    for (uint64_t i = 0; i < num_keys_to_check; i++) {
        memcpy(current_k1, &i, 8);

        // Using 3DES with keying option 2 (K1, K2, K1)
        // We decrypt the target and see if it matches our known plaintext.
        triple_des_decrypt_block(TARGET_CIPHERTEXT, decrypted_block, current_k1, current_k2, current_k1);

        if (memcmp(decrypted_block, KNOWN_PLAINTEXT, 8) == 0) {
            (*found_count)++;
        }
    }
}

void cleanup() {
    if (found_count) {
        free(found_count);
        found_count = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    uint64_t final_result;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    final_result = *found_count;

    cleanup();

    // Print result to stdout
    printf("%llu\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
