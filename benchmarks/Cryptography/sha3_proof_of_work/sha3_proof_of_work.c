#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>

/*******************************************************************************
 *                        MERSENNE TWISTER (VERBATIM)
 ******************************************************************************/
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

/*******************************************************************************
 *                        SHA-3 (KECCAK) IMPLEMENTATION
 ******************************************************************************/
#define SHA3_STATE_SIZE_BYTES 200 // 1600 bits
#define SHA3_256_RATE_BYTES 136   // r = 1600 - 2*256 = 1088 bits = 136 bytes
#define SHA3_DIGEST_SIZE 32     // 256 bits

#define ROTL64(x, y) (((x) << (y)) | ((x) >> (64 - (y))))

typedef struct {
    uint64_t state[25]; // 25 * 64 bits = 1600 bits
    unsigned int buffer_pos;
    uint8_t buffer[SHA3_256_RATE_BYTES];
} keccak_ctx;

static const uint64_t KeccakF_RoundConstants[24] = {
    0x0000000000000001, 0x0000000000008082, 0x800000000000808a, 0x8000000080008000, 
    0x000000000000808b, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009, 
    0x000000000000008a, 0x0000000000000088, 0x0000000080008009, 0x000000008000000a, 
    0x000000008000808b, 0x800000000000008b, 0x8000000000008089, 0x8000000000008003, 
    0x8000000000008002, 0x8000000000000080, 0x000000000000800a, 0x800000008000000a, 
    0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008
};

static void keccak_f1600(uint64_t state[25]) {
    uint64_t C[5], D, B[25];
    int i, j, round;

    for (round = 0; round < 24; round++) {
        // Theta
        for (i = 0; i < 5; i++) C[i] = state[i] ^ state[i+5] ^ state[i+10] ^ state[i+15] ^ state[i+20];
        for (i = 0; i < 5; i++) {
            D = C[(i + 4) % 5] ^ ROTL64(C[(i + 1) % 5], 1);
            for (j = 0; j < 25; j += 5) state[i + j] ^= D;
        }

        // Rho and Pi
        D = state[1];
        B[0] = state[0];
        B[10] = ROTL64(state[1], 1);    B[7]  = ROTL64(state[2], 3);    B[11] = ROTL64(state[3], 6);
        B[17] = ROTL64(state[4], 10);   B[18] = ROTL64(state[5], 15);   B[3]  = ROTL64(state[6], 21);
        B[5]  = ROTL64(state[7], 28);   B[16] = ROTL64(state[8], 36);   B[8]  = ROTL64(state[9], 45);
        B[21] = ROTL64(state[10], 55);  B[24] = ROTL64(state[11], 2);   B[4]  = ROTL64(state[12], 14);
        B[15] = ROTL64(state[13], 27);  B[23] = ROTL64(state[14], 41);  B[19] = ROTL64(state[15], 56);
        B[13] = ROTL64(state[16], 8);   B[12] = ROTL64(state[17], 25);  B[2]  = ROTL64(state[18], 43);
        B[20] = ROTL64(state[19], 62);  B[14] = ROTL64(state[20], 18);  B[22] = ROTL64(state[21], 39);
        B[9]  = ROTL64(state[22], 61);  B[6]  = ROTL64(state[23], 20);  B[1]  = ROTL64(state[24], 44);
        for(i = 0; i < 25; i++) state[i] = B[i];

        // Chi
        for (j = 0; j < 25; j += 5) {
            for (i = 0; i < 5; i++) B[i] = state[j + i];
            for (i = 0; i < 5; i++) state[j+i] = B[i] ^ ((~B[(i + 1) % 5]) & B[(i + 2) % 5]);
        }

        // Iota
        state[0] ^= KeccakF_RoundConstants[round];
    }
}

static void sha3_init(keccak_ctx *ctx) {
    memset(ctx, 0, sizeof(keccak_ctx));
}

static void sha3_update(keccak_ctx *ctx, const uint8_t *data, size_t len) {
    size_t i;
    while (len > 0) {
        i = (SHA3_256_RATE_BYTES - ctx->buffer_pos < len) ? (SHA3_256_RATE_BYTES - ctx->buffer_pos) : len;
        memcpy(ctx->buffer + ctx->buffer_pos, data, i);
        ctx->buffer_pos += i;
        data += i;
        len -= i;

        if (ctx->buffer_pos == SHA3_256_RATE_BYTES) {
            for (i = 0; i < SHA3_256_RATE_BYTES / 8; i++) {
                ctx->state[i] ^= ((uint64_t*)ctx->buffer)[i];
            }
            keccak_f1600(ctx->state);
            ctx->buffer_pos = 0;
        }
    }
}

static void sha3_final(keccak_ctx *ctx, uint8_t *out) {
    // Padding
    ctx->buffer[ctx->buffer_pos++] = 0x06; // SHA-3 padding byte
    memset(ctx->buffer + ctx->buffer_pos, 0, SHA3_256_RATE_BYTES - ctx->buffer_pos);
    ctx->buffer[SHA3_256_RATE_BYTES - 1] |= 0x80;

    for (size_t i = 0; i < SHA3_256_RATE_BYTES / 8; i++) {
        ctx->state[i] ^= ((uint64_t*)ctx->buffer)[i];
    }
    keccak_f1600(ctx->state);
    memcpy(out, ctx->state, SHA3_DIGEST_SIZE);
}

static void compute_sha3_256(const uint8_t *data, size_t len, uint8_t *out) {
    keccak_ctx ctx;
    sha3_init(&ctx);
    sha3_update(&ctx, data, len);
    sha3_final(&ctx, out);
}

/*******************************************************************************
 *                        BENCHMARK SETUP & EXECUTION
 ******************************************************************************/

typedef struct {
    int difficulty_leading_zeros; // Interpreted as number of zero-nibbles
    size_t data_block_size;
    uint8_t *data_block;
    uint64_t final_nonce;
} benchmark_config;

static benchmark_config config;

static int check_difficulty(const uint8_t *hash, int leading_zero_nibbles) {
    int full_zero_bytes = leading_zero_nibbles / 2;
    for (int i = 0; i < full_zero_bytes; i++) {
        if (hash[i] != 0) {
            return 0;
        }
    }
    if (leading_zero_nibbles % 2 != 0) {
        if ((hash[full_zero_bytes] & 0xF0) != 0) {
            return 0;
        }
    }
    return 1;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <difficulty_leading_zeros> <data_block_size> <seed>\n", argv[0]);
        exit(1);
    }

    config.difficulty_leading_zeros = atoi(argv[1]);
    config.data_block_size = strtoul(argv[2], NULL, 10);
    uint32_t seed = strtoul(argv[3], NULL, 10);
    mt_seed(seed);

    config.data_block = (uint8_t*) malloc(config.data_block_size);
    if (!config.data_block) {
        fprintf(stderr, "Failed to allocate memory for data block\n");
        exit(1);
    }

    for (size_t i = 0; i < config.data_block_size; i++) {
        config.data_block[i] = mt_rand() & 0xFF;
    }
    config.final_nonce = 0;
}

void run_computation() {
    size_t buffer_size = config.data_block_size + sizeof(uint64_t);
    uint8_t *work_buffer = (uint8_t*)malloc(buffer_size);
    if (!work_buffer) {
        fprintf(stderr, "Failed to allocate memory for work buffer\n");
        exit(1);
    }
    uint8_t hash[SHA3_DIGEST_SIZE];
    
    memcpy(work_buffer, config.data_block, config.data_block_size);

    uint64_t nonce = 0;
    while(1) {
        memcpy(work_buffer + config.data_block_size, &nonce, sizeof(uint64_t));
        compute_sha3_256(work_buffer, buffer_size, hash);

        if (check_difficulty(hash, config.difficulty_leading_zeros)) {
            config.final_nonce = nonce;
            break;
        }
        nonce++;
        if (nonce == 0) {
            fprintf(stderr, "FATAL: Nonce overflowed without finding solution.\n");
            config.final_nonce = -1; // Indicate failure
            break;
        }
    }
    free(work_buffer);
}

void cleanup() {
    free(config.data_block);
}

/*******************************************************************************
 *                                  MAIN
 ******************************************************************************/
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%" PRIu64 "\n", config.final_nonce);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
