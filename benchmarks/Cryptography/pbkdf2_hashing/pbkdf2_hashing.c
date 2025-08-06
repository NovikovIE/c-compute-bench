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

// --- Cryptography Implementation (SHA-256, HMAC, PBKDF2) ---

#define SHA256_BLOCK_SIZE 64
#define SHA256_DIGEST_LENGTH 32

typedef struct {
    uint8_t data[SHA256_BLOCK_SIZE];
    uint32_t datalen;
    uint64_t bitlen;
    uint32_t state[8];
} SHA256_CTX;

static const uint32_t k[64] = {
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
};

#define ROTR(x, n) (((x) >> (n)) | ((x) << (32 - (n))))
#define CH(x, y, z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define SIGMA0(x) (ROTR(x, 2) ^ ROTR(x, 13) ^ ROTR(x, 22))
#define SIGMA1(x) (ROTR(x, 6) ^ ROTR(x, 11) ^ ROTR(x, 25))
#define sigma0(x) (ROTR(x, 7) ^ ROTR(x, 18) ^ ((x) >> 3))
#define sigma1(x) (ROTR(x, 17) ^ ROTR(x, 19) ^ ((x) >> 10))

static void sha256_transform(SHA256_CTX *ctx, const uint8_t data[]) {
    uint32_t a, b, c, d, e, f, g, h, w[64], t1, t2;
    a = ctx->state[0]; b = ctx->state[1]; c = ctx->state[2]; d = ctx->state[3];
    e = ctx->state[4]; f = ctx->state[5]; g = ctx->state[6]; h = ctx->state[7];

    for (int i = 0, j = 0; i < 16; ++i, j += 4)
        w[i] = (data[j] << 24) | (data[j + 1] << 16) | (data[j + 2] << 8) | (data[j + 3]);

    for (int i = 16; i < 64; ++i)
        w[i] = sigma1(w[i - 2]) + w[i - 7] + sigma0(w[i - 15]) + w[i - 16];

    for (int i = 0; i < 64; ++i) {
        t1 = h + SIGMA1(e) + CH(e, f, g) + k[i] + w[i];
        t2 = SIGMA0(a) + MAJ(a, b, c);
        h = g; g = f; f = e; e = d + t1;
        d = c; c = b; b = a; a = t1 + t2;
    }
    ctx->state[0] += a; ctx->state[1] += b; ctx->state[2] += c; ctx->state[3] += d;
    ctx->state[4] += e; ctx->state[5] += f; ctx->state[6] += g; ctx->state[7] += h;
}

static void sha256_init(SHA256_CTX *ctx) {
    ctx->datalen = 0;
    ctx->bitlen = 0;
    ctx->state[0] = 0x6a09e667; ctx->state[1] = 0xbb67ae85; ctx->state[2] = 0x3c6ef372; ctx->state[3] = 0xa54ff53a;
    ctx->state[4] = 0x510e527f; ctx->state[5] = 0x9b05688c; ctx->state[6] = 0x1f83d9ab; ctx->state[7] = 0x5be0cd19;
}

static void sha256_update(SHA256_CTX *ctx, const uint8_t data[], size_t len) {
    for (size_t i = 0; i < len; ++i) {
        ctx->data[ctx->datalen] = data[i];
        ctx->datalen++;
        if (ctx->datalen == SHA256_BLOCK_SIZE) {
            sha256_transform(ctx, ctx->data);
            ctx->bitlen += 512;
            ctx->datalen = 0;
        }
    }
}

static void sha256_final(SHA256_CTX *ctx, uint8_t hash[]) {
    uint32_t i = ctx->datalen;
    ctx->data[i++] = 0x80;

    if (ctx->datalen >= 56) {
        memset(ctx->data + i, 0, SHA256_BLOCK_SIZE - i);
        sha256_transform(ctx, ctx->data);
        i = 0;
    }

    memset(ctx->data + i, 0, SHA256_BLOCK_SIZE - i - 8);

    ctx->bitlen += ctx->datalen * 8;
    ctx->data[56] = (ctx->bitlen >> 56) & 0xFF; ctx->data[57] = (ctx->bitlen >> 48) & 0xFF;
    ctx->data[58] = (ctx->bitlen >> 40) & 0xFF; ctx->data[59] = (ctx->bitlen >> 32) & 0xFF;
    ctx->data[60] = (ctx->bitlen >> 24) & 0xFF; ctx->data[61] = (ctx->bitlen >> 16) & 0xFF;
    ctx->data[62] = (ctx->bitlen >> 8) & 0xFF;  ctx->data[63] = (ctx->bitlen) & 0xFF;

    sha256_transform(ctx, ctx->data);

    for (int j = 0; j < 8; ++j) {
        hash[j * 4 + 0] = (ctx->state[j] >> 24) & 0xFF;
        hash[j * 4 + 1] = (ctx->state[j] >> 16) & 0xFF;
        hash[j * 4 + 2] = (ctx->state[j] >> 8) & 0xFF;
        hash[j * 4 + 3] = (ctx->state[j]) & 0xFF;
    }
}

static void hmac_sha256(const uint8_t *key, size_t keylen, const uint8_t *data, size_t datalen, uint8_t *out) {
    SHA256_CTX ctx;
    uint8_t k_ipad[SHA256_BLOCK_SIZE];
    uint8_t k_opad[SHA256_BLOCK_SIZE];
    uint8_t temp_key[SHA256_DIGEST_LENGTH];

    if (keylen > SHA256_BLOCK_SIZE) {
        sha256_init(&ctx);
        sha256_update(&ctx, key, keylen);
        sha256_final(&ctx, temp_key);
        key = temp_key;
        keylen = SHA256_DIGEST_LENGTH;
    }

    memset(k_ipad, 0x36, SHA256_BLOCK_SIZE);
    memset(k_opad, 0x5c, SHA256_BLOCK_SIZE);
    memcpy(k_ipad, key, keylen);
    memcpy(k_opad, key, keylen);

    for(int i = 0; i < SHA256_BLOCK_SIZE; ++i) {
        k_ipad[i] ^= 0x36;
        k_opad[i] ^= 0x5c;
    }

    sha256_init(&ctx);
    sha256_update(&ctx, k_ipad, SHA256_BLOCK_SIZE);
    sha256_update(&ctx, data, datalen);
    sha256_final(&ctx, out);

    sha256_init(&ctx);
    sha256_update(&ctx, k_opad, SHA256_BLOCK_SIZE);
    sha256_update(&ctx, out, SHA256_DIGEST_LENGTH);
    sha256_final(&ctx, out);
}

static void pbkdf2_F(const uint8_t *password, size_t password_len, const uint8_t *salt, size_t salt_len,
                       unsigned int iterations, const uint8_t i_be[4], uint8_t *out_T) {
    uint8_t temp_U[SHA256_DIGEST_LENGTH];
    uint8_t salt_plus_i[salt_len + 4];

    memcpy(salt_plus_i, salt, salt_len);
    memcpy(salt_plus_i + salt_len, i_be, 4);

    hmac_sha256(password, password_len, salt_plus_i, salt_len + 4, temp_U);
    memcpy(out_T, temp_U, SHA256_DIGEST_LENGTH);

    for (unsigned int c = 1; c < iterations; ++c) {
        hmac_sha256(password, password_len, temp_U, SHA256_DIGEST_LENGTH, temp_U);
        for (int k = 0; k < SHA256_DIGEST_LENGTH; ++k) {
            out_T[k] ^= temp_U[k];
        }
    }
}

static void pbkdf2_hmac_sha256(const uint8_t *password, size_t password_len, const uint8_t *salt, size_t salt_len,
                               unsigned int iterations, uint8_t *dk, size_t dk_len) {
    size_t l = (dk_len + SHA256_DIGEST_LENGTH - 1) / SHA256_DIGEST_LENGTH;
    size_t r = dk_len - (l - 1) * SHA256_DIGEST_LENGTH;
    uint8_t T[SHA256_DIGEST_LENGTH];

    for (uint32_t i = 1; i <= l; ++i) {
        uint8_t i_be[4] = {(i >> 24) & 0xff, (i >> 16) & 0xff, (i >> 8) & 0xff, i & 0xff};
        pbkdf2_F(password, password_len, salt, salt_len, iterations, i_be, T);

        size_t copy_len = (i == l) ? r : SHA256_DIGEST_LENGTH;
        memcpy(dk + (i - 1) * SHA256_DIGEST_LENGTH, T, copy_len);
    }
}

// --- Benchmark Globals & Functions ---

unsigned int password_length = 0;
unsigned int salt_length = 0;
unsigned int iterations = 0;

unsigned char *password = NULL;
unsigned char *salt = NULL;
unsigned char *derived_key = NULL;

long final_result_sum = 0;

// Fixed output key length for this benchmark
#define DERIVED_KEY_LENGTH 32

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <password_length> <salt_length> <iterations> <seed>\n", argv[0]);
        exit(1);
    }

    password_length = atoi(argv[1]);
    salt_length = atoi(argv[2]);
    iterations = atoi(argv[3]);
    uint32_t seed = atol(argv[4]);
    mt_seed(seed);

    password = (unsigned char *)malloc(password_length);
    salt = (unsigned char *)malloc(salt_length);
    derived_key = (unsigned char *)malloc(DERIVED_KEY_LENGTH);

    if (!password || !salt || !derived_key) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (unsigned int i = 0; i < password_length; ++i) {
        password[i] = mt_rand() & 0xFF;
    }
    for (unsigned int i = 0; i < salt_length; ++i) {
        salt[i] = mt_rand() & 0xFF;
    }
}

void run_computation() {
    pbkdf2_hmac_sha256(password, password_length, salt, salt_length, iterations, derived_key, DERIVED_KEY_LENGTH);

    final_result_sum = 0;
    for (int i = 0; i < DERIVED_KEY_LENGTH; ++i) {
        final_result_sum += derived_key[i];
    }
}

void cleanup() {
    free(password);
    free(salt);
    free(derived_key);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%ld\n", final_result_sum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
