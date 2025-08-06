#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// --- MERSENNE TWISTER (Verbatim) ---
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

// --- CRYPTO HELPERS ---

static inline uint32_t le32dec(const void *pp) {
    const uint8_t *p = (uint8_t const *)pp;
    return ((uint32_t)(p[3]) << 24) | ((uint32_t)(p[2]) << 16) |
           ((uint32_t)(p[1]) << 8) | p[0];
}

// --- SHA256 IMPLEMENTATION ---
#define SHA256_BLOCK_SIZE 64
#define SHA256_DIGEST_SIZE 32
#define ROTR(a,b) (((a) >> (b)) | ((a) << (32-(b))))

typedef struct {
	uint8_t data[SHA256_BLOCK_SIZE];
	uint32_t datalen;
	unsigned long long bitlen;
	uint32_t state[8];
} SHA256_CTX;

static const uint32_t k[64] = {
	0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
	0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
	0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
	0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
	0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
	0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
	0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
	0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

static void sha256_transform(SHA256_CTX *ctx, const uint8_t data[]) {
	uint32_t a,b,c,d,e,f,g,h,i,j,m[64];

	for (i=0,j=0; i < 16; ++i, j += 4)
		m[i] = (data[j] << 24) | (data[j+1] << 16) | (data[j+2] << 8) | (data[j+3]);
	for ( ; i < 64; ++i)
		m[i] = ROTR(m[i-2],17) ^ ROTR(m[i-2],19) ^ (m[i-2] >> 10)   + m[i-7] +   ROTR(m[i-15],7) ^ ROTR(m[i-15],18) ^ (m[i-15] >> 3)   + m[i-16];

	a = ctx->state[0]; b = ctx->state[1]; c = ctx->state[2]; d = ctx->state[3];
	e = ctx->state[4]; f = ctx->state[5]; g = ctx->state[6]; h = ctx->state[7];

	for (i = 0; i < 64; ++i) {
		uint32_t t1 = h + (ROTR(e,6) ^ ROTR(e,11) ^ ROTR(e,25)) + ((e & f) ^ (~e & g)) + k[i] + m[i];
		uint32_t t2 = (ROTR(a,2) ^ ROTR(a,13) ^ ROTR(a,22)) + ((a & b) ^ (a & c) ^ (b & c));
		h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
	}

	ctx->state[0] += a; ctx->state[1] += b; ctx->state[2] += c; ctx->state[3] += d;
	ctx->state[4] += e; ctx->state[5] += f; ctx->state[6] += g; ctx->state[7] += h;
}

static void sha256_init(SHA256_CTX *ctx) {
	ctx->datalen = 0; ctx->bitlen = 0;
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
	if (ctx->datalen < 56) {
		memset(ctx->data + i, 0, 56 - i);
	} else {
		memset(ctx->data + i, 0, SHA256_BLOCK_SIZE - i);
		sha256_transform(ctx, ctx->data);
		memset(ctx->data, 0, 56);
	}

	ctx->bitlen += ctx->datalen * 8;
	ctx->data[63] = ctx->bitlen;
	ctx->data[62] = ctx->bitlen >> 8;
	ctx->data[61] = ctx->bitlen >> 16;
	ctx->data[60] = ctx->bitlen >> 24;
	ctx->data[59] = ctx->bitlen >> 32;
	ctx->data[58] = ctx->bitlen >> 40;
	ctx->data[57] = ctx->bitlen >> 48;
	ctx->data[56] = ctx->bitlen >> 56;
	sha256_transform(ctx, ctx->data);

	for (i = 0; i < 8; ++i) {
		hash[i*4+0] = (ctx->state[i] >> 24) & 0xff;
		hash[i*4+1] = (ctx->state[i] >> 16) & 0xff;
		hash[i*4+2] = (ctx->state[i] >> 8) & 0xff;
		hash[i*4+3] = (ctx->state[i] >> 0) & 0xff;
	}
}

// --- HMAC-SHA256 ---
static void hmac_sha256(const uint8_t *key, size_t keylen, const uint8_t *in, size_t inlen, uint8_t out[SHA256_DIGEST_SIZE]) {
	SHA256_CTX ctx;
	uint8_t k[SHA256_BLOCK_SIZE];
	uint8_t temp[SHA256_DIGEST_SIZE];
	
	if (keylen > SHA256_BLOCK_SIZE) {
		sha256_init(&ctx);
		sha256_update(&ctx, key, keylen);
		sha256_final(&ctx, k);
		memset(k + SHA256_DIGEST_SIZE, 0, SHA256_BLOCK_SIZE - SHA256_DIGEST_SIZE);
	} else {
		memcpy(k, key, keylen);
		memset(k + keylen, 0, SHA256_BLOCK_SIZE - keylen);
	}

	uint8_t k_ipad[SHA256_BLOCK_SIZE];
	for(size_t i = 0; i < SHA256_BLOCK_SIZE; i++) k_ipad[i] = k[i] ^ 0x36;
	
	sha256_init(&ctx);
	sha256_update(&ctx, k_ipad, SHA256_BLOCK_SIZE);
	sha256_update(&ctx, in, inlen);
	sha256_final(&ctx, temp);

	uint8_t k_opad[SHA256_BLOCK_SIZE];
	for(size_t i = 0; i < SHA256_BLOCK_SIZE; i++) k_opad[i] = k[i] ^ 0x5c;

	sha256_init(&ctx);
	sha256_update(&ctx, k_opad, SHA256_BLOCK_SIZE);
	sha256_update(&ctx, temp, SHA256_DIGEST_SIZE);
	sha256_final(&ctx, out);
}

// --- PBKDF2-HMAC-SHA256 ---
static void pbkdf2_hmac_sha256(const uint8_t *passwd, size_t passwdlen, const uint8_t *salt, size_t saltlen, uint32_t c, uint8_t *dk, size_t dklen) {
	uint8_t U[SHA256_DIGEST_SIZE];
	uint8_t T[SHA256_DIGEST_SIZE];
	uint8_t salt_plus[saltlen + 4];
	size_t dk_pos = 0;

	for (uint32_t i = 1; dk_pos < dklen; i++) {
		memcpy(salt_plus, salt, saltlen);
		salt_plus[saltlen] = (i >> 24) & 0xFF;
		salt_plus[saltlen+1] = (i >> 16) & 0xFF;
		salt_plus[saltlen+2] = (i >> 8) & 0xFF;
		salt_plus[saltlen+3] = i & 0xFF;
		
		hmac_sha256(passwd, passwdlen, salt_plus, saltlen + 4, U);
		memcpy(T, U, SHA256_DIGEST_SIZE);

		for (uint32_t j = 1; j < c; j++) {
			hmac_sha256(passwd, passwdlen, U, SHA256_DIGEST_SIZE, U);
			for (int k = 0; k < SHA256_DIGEST_SIZE; k++) {
				T[k] ^= U[k];
			}
		}
		
		size_t copy_len = (dklen - dk_pos) < SHA256_DIGEST_SIZE ? (dklen - dk_pos) : SHA256_DIGEST_SIZE;
		memcpy(dk + dk_pos, T, copy_len);
		dk_pos += copy_len;
	}
}

// --- SCRYPT IMPLEMENTATION ---
#define R_SALSA(a,b) (((a) << (b)) | ((a) >> (32 - (b))))
static void salsa20_8_core(uint32_t B[16]) {
    uint32_t x[16];
    memcpy(x, B, sizeof(x));
    for (int i = 0; i < 8; i += 2) {
        x[ 4] ^= R_SALSA(x[ 0]+x[12], 7);  x[ 8] ^= R_SALSA(x[ 4]+x[ 0], 9);
        x[12] ^= R_SALSA(x[ 8]+x[ 4],13);  x[ 0] ^= R_SALSA(x[12]+x[ 8],18);
        x[ 9] ^= R_SALSA(x[ 5]+x[ 1], 7);  x[13] ^= R_SALSA(x[ 9]+x[ 5], 9);
        x[ 1] ^= R_SALSA(x[13]+x[ 9],13);  x[ 5] ^= R_SALSA(x[ 1]+x[13],18);
        x[14] ^= R_SALSA(x[10]+x[ 6], 7);  x[ 2] ^= R_SALSA(x[14]+x[10], 9);
        x[ 6] ^= R_SALSA(x[ 2]+x[14],13);  x[10] ^= R_SALSA(x[ 6]+x[ 2],18);
        x[ 3] ^= R_SALSA(x[15]+x[11], 7);  x[ 7] ^= R_SALSA(x[ 3]+x[15], 9);
        x[11] ^= R_SALSA(x[ 7]+x[ 3],13);  x[15] ^= R_SALSA(x[11]+x[ 7],18);
        x[ 1] ^= R_SALSA(x[ 0]+x[ 3], 7);  x[ 2] ^= R_SALSA(x[ 1]+x[ 0], 9);
        x[ 3] ^= R_SALSA(x[ 2]+x[ 1],13);  x[ 0] ^= R_SALSA(x[ 3]+x[ 2],18);
        x[ 6] ^= R_SALSA(x[ 5]+x[ 4], 7);  x[ 7] ^= R_SALSA(x[ 6]+x[ 5], 9);
        x[ 4] ^= R_SALSA(x[ 7]+x[ 6],13);  x[ 5] ^= R_SALSA(x[ 4]+x[ 7],18);
        x[11] ^= R_SALSA(x[10]+x[ 9], 7);  x[ 8] ^= R_SALSA(x[11]+x[10], 9);
        x[ 9] ^= R_SALSA(x[ 8]+x[11],13);  x[10] ^= R_SALSA(x[ 9]+x[ 8],18);
        x[12] ^= R_SALSA(x[15]+x[14], 7);  x[13] ^= R_SALSA(x[12]+x[15], 9);
        x[14] ^= R_SALSA(x[13]+x[12],13);  x[15] ^= R_SALSA(x[14]+x[13],18);
    }
    for(int i = 0; i < 16; ++i) B[i] += x[i];
}

static void blockmix_salsa8(const uint32_t *Bin, uint32_t *Bout, size_t r) {
    uint32_t T[16];
    memcpy(T, &Bin[(2 * r - 1) * 16], sizeof(T));

    for (size_t i = 0; i < 2 * r; i++) {
        for(int j=0; j<16; j++) T[j] ^= Bin[i * 16 + j];
        salsa20_8_core(T);
        if (i % 2 == 0) {
            memcpy(&Bout[(i/2) * 16], T, sizeof(T));
        } else {
            memcpy(&Bout[(r + i/2) * 16], T, sizeof(T));
        }
    }
}

static void romix(uint8_t *B, uint32_t N, uint32_t r, uint8_t *V) {
    uint8_t *X = V;
    uint8_t *T = V + (128*r);
    
    memcpy(X, B, 128 * r);

    for (uint32_t i = 0; i < N; i++) {
        memcpy(&V[i * 128 * r], X, 128 * r);
        blockmix_salsa8((uint32_t*)X, (uint32_t*)X, r);
    }

    for (uint32_t i = 0; i < N; i++) {
        uint32_t j = le32dec(&X[(2 * r - 1) * 64]) & (N - 1);
        uint8_t *Vj = &V[j * 128 * r];
        for(uint32_t k = 0; k < 128 * r; k++) T[k] = X[k] ^ Vj[k];
        blockmix_salsa8((uint32_t*)T, (uint32_t*)X, r);
    }
    memcpy(B, X, 128 * r);
}

static void scrypt_prealloc(const uint8_t *passwd, size_t passwdlen, const uint8_t *salt, size_t saltlen,
                uint32_t N, uint32_t r, uint32_t p,
                uint8_t *dk, size_t dklen, uint8_t *B_buf, uint8_t *V_buf) {

    size_t B_len = (size_t)p * 128 * r;
    pbkdf2_hmac_sha256(passwd, passwdlen, salt, saltlen, 1, B_buf, B_len);

    for (uint32_t i = 0; i < p; i++) {
        romix(&B_buf[i * 128 * r], N, r, V_buf);
    }

    pbkdf2_hmac_sha256(passwd, passwdlen, B_buf, B_len, 1, dk, dklen);
}


// --- BENCHMARK STRUCTURE ---

typedef struct {
    uint32_t N, r, p;
    uint8_t *password;
    size_t password_len;
    uint8_t *salt;
    size_t salt_len;
    uint8_t *derived_key;
    size_t derived_key_len;
    uint8_t *B_buf; // Workspace for PBKDF2 intermediate
    uint8_t *V_buf; // Workspace for ROMix
    long long final_result;
} BenchmarkData;

static BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <cpu_memory_cost_n> <block_size_r> <parallelization_p> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.N = atoi(argv[1]);
    g_data.r = atoi(argv[2]);
    g_data.p = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);
    mt_seed(seed);

    if ((g_data.N & (g_data.N - 1)) != 0 || g_data.N == 0) {
        fprintf(stderr, "ERROR: N must be a power of 2\n");
        exit(1);
    }
    if ((uint64_t)g_data.r * g_data.p >= (1 << 30)) {
        fprintf(stderr, "ERROR: r * p is too large\n");
        exit(1);
    }
    
    g_data.password_len = 64;
    g_data.salt_len = 32;
    g_data.derived_key_len = 64;

    g_data.password = (uint8_t*)malloc(g_data.password_len);
    g_data.salt = (uint8_t*)malloc(g_data.salt_len);
    g_data.derived_key = (uint8_t*)malloc(g_data.derived_key_len);

    if (!g_data.password || !g_data.salt || !g_data.derived_key) {
        fprintf(stderr, "ERROR: Failed to allocate password/salt/dk buffers\n");
        exit(1);
    }

    for (size_t i = 0; i < g_data.password_len; i++) g_data.password[i] = mt_rand() & 0xFF;
    for (size_t i = 0; i < g_data.salt_len; i++) g_data.salt[i] = mt_rand() & 0xFF;

    size_t B_size = (size_t)128 * g_data.r * g_data.p;
    // ROMix needs workspace for V, plus a copy of X and T. V is N * 128 * r.
    size_t V_size = (size_t)128 * g_data.r * g_data.N + 2 * 128 * g_data.r;
    g_data.B_buf = (uint8_t*)malloc(B_size);
    g_data.V_buf = (uint8_t*)malloc(V_size);
    if (!g_data.B_buf || !g_data.V_buf) {
        fprintf(stderr, "ERROR: Failed to allocate scrypt workspace\n");
        exit(1);
    }

    g_data.final_result = 0;
}

void run_computation() {
    scrypt_prealloc(g_data.password, g_data.password_len,
                    g_data.salt, g_data.salt_len,
                    g_data.N, g_data.r, g_data.p,
                    g_data.derived_key, g_data.derived_key_len,
                    g_data.B_buf, g_data.V_buf);
    
    // Accumulate result to prevent dead code elimination
    long long sum = 0;
    for (size_t i = 0; i < g_data.derived_key_len; i++) {
        sum += g_data.derived_key[i];
    }
    g_data.final_result = sum;
}

void cleanup() {
    free(g_data.password);
    free(g_data.salt);
    free(g_data.derived_key);
    free(g_data.B_buf);
    free(g_data.V_buf);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
