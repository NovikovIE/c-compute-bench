#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// MERSENNE TWISTER (Original C code by Makoto Matsumoto and Takuji Nishimura)
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

// BENCHMARK PARAMETERS & DATA
#define NUM_VERIFICATIONS 250 // Tuned for ~1s runtime with 256-bit curves

int CURVE_BIT_LENGTH;
int MESSAGE_SIZE_BYTES;
int BIGNUM_WORDS;

// ECDSA Data
uint8_t **messages;
uint32_t **pub_keys_x, **pub_keys_y;
uint32_t **signatures_r, **signatures_s;
uint32_t *curve_n; // Group order
uint32_t *gen_point_x, *gen_point_y; // Generator Point G

int final_result;

// --- SIMULATED CRYPTOGRAPHIC HELPER FUNCTIONS ---
// These functions model the computational complexity profile of real
// cryptographic operations, but are not cryptographically correct.

// res = (a + b)
void bignum_add_sim(uint32_t* res, const uint32_t* a, const uint32_t* b) {
    uint64_t carry = 0;
    for (int i = 0; i < BIGNUM_WORDS; i++) {
        uint64_t sum = (uint64_t)a[i] + b[i] + carry;
        res[i] = (uint32_t)sum;
        carry = sum >> 32;
    }
}

// res = (a * b)
void bignum_mul_sim(uint32_t* res, const uint32_t* a, const uint32_t* b) {
    uint64_t* temp = (uint64_t*)calloc(BIGNUM_WORDS * 2, sizeof(uint64_t));
    if (!temp) exit(1);

    for (int i = 0; i < BIGNUM_WORDS; i++) {
        for (int j = 0; j < BIGNUM_WORDS; j++) {
            temp[i + j] += (uint64_t)a[i] * b[j];
        }
    }

    for (int i = 0; i < BIGNUM_WORDS * 2 - 1; i++) {
        temp[i + 1] += temp[i] >> 32;
    }

    for (int i = 0; i < BIGNUM_WORDS; i++) {
        res[i] = (uint32_t)temp[i];
    }

    free(temp);
}

// res = a^-1 mod m (simulated with O(N^2) complexity)
void bignum_mod_inv_sim(uint32_t* res, const uint32_t* a, const uint32_t* m) {
    uint32_t* temp1 = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t)); 
    uint32_t* temp2 = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    if (!temp1 || !temp2) exit(1);

    memcpy(temp1, a, BIGNUM_WORDS * sizeof(uint32_t));

    for (int i = 0; i < BIGNUM_WORDS; i++) {
        bignum_mul_sim(temp2, temp1, temp1);
        bignum_add_sim(temp1, temp2, m);
    }
    memcpy(res, temp1, BIGNUM_WORDS * sizeof(uint32_t));

    free(temp1);
    free(temp2);
}

// point_r = point_p + point_q (simulated)
void ec_point_add_sim(uint32_t* rx, uint32_t* ry, const uint32_t* px, const uint32_t* py, const uint32_t* qx, const uint32_t* qy) {
    uint32_t *t1 = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t *t2 = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    if (!t1 || !t2) exit(1);

    bignum_add_sim(t1, py, qy); // Simulate dy
    bignum_add_sim(t2, px, qx); // Simulate dx
    bignum_mod_inv_sim(t2, t2, curve_n); // Simulate 1/dx
    bignum_mul_sim(t1, t1, t2); // s = dy/dx

    bignum_mul_sim(t2, t1, t1); // s^2
    bignum_add_sim(t2, t2, px); // s^2 - px
    bignum_add_sim(rx, t2, qx); // s^2 - px - qx

    bignum_add_sim(t2, px, rx); // px - rx
    bignum_mul_sim(t1, t1, t2); // s(px-rx)
    bignum_add_sim(ry, t1, py); // s(px-rx) - py
    
    free(t1); free(t2);
}

// point_r = k * point_p (simulated double-and-add)
void ec_scalar_mul_sim(uint32_t* rx, uint32_t* ry, const uint32_t* k, const uint32_t* px, const uint32_t* py) {
    uint32_t* temp_px = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* temp_py = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    if (!temp_px || !temp_py) exit(1);
    
    memcpy(temp_px, px, BIGNUM_WORDS * sizeof(uint32_t));
    memcpy(temp_py, py, BIGNUM_WORDS * sizeof(uint32_t));
    memset(rx, 0, BIGNUM_WORDS * sizeof(uint32_t)); // Point at infinity
    memset(ry, 0, BIGNUM_WORDS * sizeof(uint32_t));

    for (int i = 0; i < CURVE_BIT_LENGTH; i++) {
        int word_idx = i / 32;
        int bit_idx = i % 32;
        if (((k[word_idx] >> bit_idx) & 1)) {
            ec_point_add_sim(rx, ry, rx, ry, temp_px, temp_py);
        }
        ec_point_add_sim(temp_px, temp_py, temp_px, temp_py, temp_px, temp_py);
    }

    free(temp_px); free(temp_py);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <curve_bit_length> <message_size_bytes> <seed>\n", argv[0]);
        exit(1);
    }

    CURVE_BIT_LENGTH = atoi(argv[1]);
    MESSAGE_SIZE_BYTES = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    mt_seed(seed);

    BIGNUM_WORDS = (CURVE_BIT_LENGTH + 31) / 32;
    if (BIGNUM_WORDS <= 0) exit(1);

    // Allocate arrays of pointers
    messages = (uint8_t**)malloc(NUM_VERIFICATIONS * sizeof(uint8_t*));
    pub_keys_x = (uint32_t**)malloc(NUM_VERIFICATIONS * sizeof(uint32_t*));
    pub_keys_y = (uint32_t**)malloc(NUM_VERIFICATIONS * sizeof(uint32_t*));
    signatures_r = (uint32_t**)malloc(NUM_VERIFICATIONS * sizeof(uint32_t*));
    signatures_s = (uint32_t**)malloc(NUM_VERIFICATIONS * sizeof(uint32_t*));
    if (!messages || !pub_keys_x || !pub_keys_y || !signatures_r || !signatures_s) exit(1);

    for (int i = 0; i < NUM_VERIFICATIONS; i++) {
        messages[i] = (uint8_t*)malloc(MESSAGE_SIZE_BYTES * sizeof(uint8_t));
        pub_keys_x[i] = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
        pub_keys_y[i] = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
        signatures_r[i] = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
        signatures_s[i] = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
        if (!messages[i] || !pub_keys_x[i] || !pub_keys_y[i] || !signatures_r[i] || !signatures_s[i]) exit(1);

        for (int j = 0; j < MESSAGE_SIZE_BYTES; j++) messages[i][j] = mt_rand() & 0xFF;
        for (int j = 0; j < BIGNUM_WORDS; j++) pub_keys_x[i][j] = mt_rand();
        for (int j = 0; j < BIGNUM_WORDS; j++) pub_keys_y[i][j] = mt_rand();
        for (int j = 0; j < BIGNUM_WORDS; j++) signatures_r[i][j] = mt_rand();
        for (int j = 0; j < BIGNUM_WORDS; j++) signatures_s[i][j] = mt_rand();
    }

    curve_n = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    gen_point_x = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    gen_point_y = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    if (!curve_n || !gen_point_x || !gen_point_y) exit(1);

    for (int j = 0; j < BIGNUM_WORDS; j++) curve_n[j] = mt_rand();
    for (int j = 0; j < BIGNUM_WORDS; j++) gen_point_x[j] = mt_rand();
    for (int j = 0; j < BIGNUM_WORDS; j++) gen_point_y[j] = mt_rand();
}

void run_computation() {
    int valid_signatures = 0;
    uint32_t* h = (uint32_t*)calloc(BIGNUM_WORDS, sizeof(uint32_t));
    uint32_t* s_inv = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* u1 = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* u2 = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* p1x = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* p1y = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* p2x = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* p2y = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* final_x = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    uint32_t* final_y = (uint32_t*)malloc(BIGNUM_WORDS * sizeof(uint32_t));
    if (!h || !s_inv || !u1 || !u2 || !p1x || !p1y || !p2x || !p2y || !final_x || !final_y) exit(1);

    for (int i = 0; i < NUM_VERIFICATIONS; i++) {
        // 1. Hash message (simplified)
        uint32_t hash_val = 0;
        for (int j = 0; j < MESSAGE_SIZE_BYTES; j++) hash_val = (hash_val * 31) + messages[i][j];
        memset(h, 0, BIGNUM_WORDS * sizeof(uint32_t));
        h[0] = hash_val;

        // 2. s_inv = s^-1 mod n
        bignum_mod_inv_sim(s_inv, signatures_s[i], curve_n);

        // 3. u1 = h * s_inv mod n; u2 = r * s_inv mod n
        bignum_mul_sim(u1, h, s_inv); 
        bignum_mul_sim(u2, signatures_r[i], s_inv);

        // 4. P1 = u1*G
        ec_scalar_mul_sim(p1x, p1y, u1, gen_point_x, gen_point_y);
        
        // 5. P2 = u2*Q
        ec_scalar_mul_sim(p2x, p2y, u2, pub_keys_x[i], pub_keys_y[i]);

        // 6. R = P1 + P2
        ec_point_add_sim(final_x, final_y, p1x, p1y, p2x, p2y);

        // 7. verification: final_x mod n == r mod n (simplified comparison)
        if ((final_x[0] & 0xFF) == (signatures_r[i][0] & 0xFF)) {
            valid_signatures++;
        }
    }

    final_result = valid_signatures;

    free(h); free(s_inv); free(u1); free(u2);
    free(p1x); free(p1y); free(p2x); free(p2y);
    free(final_x); free(final_y);
}

void cleanup() {
    free(curve_n);
    free(gen_point_x);
    free(gen_point_y);

    for (int i = 0; i < NUM_VERIFICATIONS; i++) {
        free(messages[i]);
        free(pub_keys_x[i]);
        free(pub_keys_y[i]);
        free(signatures_r[i]);
        free(signatures_s[i]);
    }
    free(messages);
    free(pub_keys_x);
    free(pub_keys_y);
    free(signatures_r);
    free(signatures_s);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
