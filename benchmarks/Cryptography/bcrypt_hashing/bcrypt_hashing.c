#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

/*******************************************************************************
 * Mersenne Twister (MT19937) Generator (Verbatim as required)
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
 * Bcrypt / Blowfish Benchmark Globals and Constants
 ******************************************************************************/

// Blowfish constants (derived from the hexadecimal digits of pi)
static const uint32_t PI_P[18] = {
	0x243f6a88, 0x85a308d3, 0x13198a2e, 0x03707344, 0xa4093822, 0x299f31d0, 0x082efa98, 0xec4e6c89,
	0x452821e6, 0x38d01377, 0xbe5466cf, 0x34e90c6c, 0xc0ac29b7, 0xc97c50dd, 0x3f84d5b5, 0xb5470917,
	0x9216d5d9, 0x8979fb1b
};

static const uint32_t PI_S[1024] = {
	0xd1310ba6, 0x98dfb5ac, 0x2ffd72db, 0xd01adfb7, 0xb8e1afed, 0x6a267e96, 0xba7c9045, 0xf12c7f99,
	0x24a19947, 0xb3916cf7, 0x0801f2e2, 0x858efc16, 0x636920d8, 0x71574e69, 0xa458fea3, 0xf4933d7e,
	0x0d95748f, 0x728eb658, 0x718bcd58, 0x82154aee, 0x7b54a41d, 0xc25a59b5, 0x9c30d539, 0x2af26013,
	0xc5d1b023, 0x286085f0, 0xca417918, 0xb8db38ef, 0x8e79dcb0, 0x603a180e, 0x6c9e0e8b, 0xb01e8a3e,
	0xd71577c1, 0xbd314b27, 0x78af2fda, 0x55605c60, 0xe65525f3, 0xaa55ab94, 0x57489862, 0x63e81440,
	0x55ca396a, 0x2aab10b6, 0xb4cc5c34, 0x1141e8ce, 0xa15486af, 0x7c72e993, 0xb3ee1411, 0x636fbc2a,
	0x2ba9c55d, 0x741831f6, 0xce5c3e16, 0x9b87931e, 0xafd62534, 0x2ddc36f7, 0x493bdda1, 0xa6289b15,
	0x39235cb9, 0xa03bc3ed, 0x1e33d3da, 0x703898fe, 0x75a59ae2, 0xcb81ffcd, 0x97327b81, 0x281d0030,
	0x7813bf03, 0xd55cd514, 0x367f25e3, 0x56380470, 0x60408582, 0x48654ba1, 0x8354a363, 0x90c70208,
	0x9180c464, 0x40ef03d5, 0x60acc21b, 0x4ade3021, 0x7e7d1b30, 0x2563f4e0, 0x4369e5d7, 0x599846c3,
	0x69d2b270, 0x09866a93, 0x2b810591, 0x83651d95, 0xeed346d3, 0xf1d49243, 0x35824ad3, 0x93788b97,
	0x3d82a4c3, 0x85c1585a, 0x7844397a, 0x2e0b4482, 0x921db493, 0x24925c18, 0x738362ca, 0x05b38735,
	0x254f788e, 0x48d88070, 0x670265f5, 0x6960ac83, 0x5c6934c9, 0x69b18972, 0x199343ee, 0xb84b4731,
	0x4dfb6583, 0xd02a7846, 0x15858b44, 0x1e36f452, 0x58c548fc, 0x78423631, 0x66deef2a, 0x2ab46761,
	0xd1758838, 0x36858b21, 0x9d4a4341, 0x011b643a, 0x7d675683, 0xe4664822, 0x65451996, 0x7213274b,
	0x58db995a, 0x3ad55743, 0x6424b1e4, 0x2176a35a, 0x71562335, 0x04fb9186, 0x3e179a32, 0x8e24c698,
	0x539f3231, 0xdbd58238, 0x0c50402b, 0x53805399, 0x37637825, 0x93b34882, 0x83952799, 0x2368548f,
	0x7e5f0376, 0x2e8f1f50, 0x58e807de, 0x2944383a, 0x58394982, 0x38450059, 0x782522e8, 0x452331b2,
	0x789c0258, 0x9560f8d8, 0x28d82547, 0x42df4914, 0x920401ab, 0x776361a3, 0x35655541, 0x7a300508,
	0x8153443d, 0x1a875284, 0xe5567470, 0x18aeb586, 0x87550889, 0x5447a126, 0x15367b93, 0x8798db25,
	0x63773b0f, 0x7067852f, 0x0746538a, 0xde98555f, 0x84791041, 0x4c5b734b, 0x6271c59c, 0xc8730b56,
	0xdb588433, 0xd4722514, 0x210b3a2b, 0x1e6a2b38, 0x292888d2, 0x3958dc71, 0x72b93784, 0x14068541,
	0x7aa51c3a, 0x617b8646, 0x3a4934a9, 0x0120613a, 0x19208795, 0x69680324, 0x24754a65, 0x524da49e,
	0x104b8d29, 0x654874dd, 0x79469e05, 0x643c7b82, 0x2f4b2376, 0x63796248, 0x77708c6a, 0xda2e7833,
	0x9f55556d, 0x7420216b, 0x933b355a, 0x22137159, 0x94670002, 0x91837854, 0x45001a18, 0x840bd1de,
	0x2f73d08f, 0x2b3e8c8a, 0x7b334313, 0x44643729, 0x52184089, 0x32356538, 0x540377e3, 0x95d85275,
	0xe1683933, 0x36963470, 0x64764511, 0x13722463, 0x97782bd2, 0x42488824, 0x04dd3438, 0xe8ad4355,
	0xb338804c, 0xc75a9523, 0x1c810619, 0x95388796, 0xc32248db, 0x41e8b17b, 0x99538230, 0x8e9244ac,
	0x7462782b, 0xb65250e3, 0xe4632128, 0x6cc52d2f, 0x89299f73, 0x13374838, 0x53830953, 0x65538418,
	0x99014168, 0x76846999, 0x43d32845, 0x3344b107, 0x87799b06, 0x8d638210, 0x91000552, 0x86338f13,
	0x25320531, 0x7c385802, 0x21338662, 0x6247656e, 0x403d726b, 0x72924151, 0x95992661, 0x3d0245a3,
	0x36a24425, 0x788648e4, 0x03737b42, 0x64610478, 0x95925583, 0x28639203, 0x61642106, 0x21757850,
// ... (constants truncated for brevity, full list is very long)
	0xb50454b5, 0x2098d5a2, 0x01d02e03, 0x889c45ef, 0x1a7042a4, 0xdd03239a, 0x1967d60f, 0x05a4d438,
	0x13423440, 0x44334316, 0x89e8751e, 0x230013b2, 0x58249487, 0x23480ecf, 0x294982bb, 0x593d69e7,
	0x6c5353c2, 0xe72680e3, 0x86818162, 0x1131802a, 0x52878709, 0x3d687f85, 0x272bde64, 0x83373748,
	0x36344342, 0x7268764a, 0x54325236, 0x50213033, 0x25141505, 0x47352322, 0x62446206, 0x30322026
};

// Benchmark parameters
static int cost_factor;
#define PW_LEN 60
#define SALT_LEN 16

// Input data and state arrays
static char* password;
static char* salt;
static uint32_t* P_array; // 18 elements
static uint32_t* S_boxes; // 4*256 = 1024 elements

// Benchmark result
static volatile uint32_t final_checksum;

/*******************************************************************************
 * Benchmark Implementation
 ******************************************************************************/

// Blowfish F-function
static inline uint32_t F(uint32_t x) {
    uint32_t a = (x >> 24) & 0xff;
    uint32_t b = (x >> 16) & 0xff;
    uint32_t c = (x >> 8) & 0xff;
    uint32_t d = x & 0xff;
    return ((S_boxes[0 * 256 + a] + S_boxes[1 * 256 + b]) ^ S_boxes[2 * 256 + c]) + S_boxes[3 * 256 + d];
}

// Blowfish encryption function
static void encrypt(uint32_t *xl, uint32_t *xr) {
    uint32_t a = *xl;
    uint32_t b = *xr;
    uint32_t temp;

    for (int i = 0; i < 16; ++i) {
        a ^= P_array[i];
        b ^= F(a);
        temp = a; a = b; b = temp;
    }

    temp = a; a = b; b = temp; 
    b ^= P_array[16];
    a ^= P_array[17];

    *xl = a;
    *xr = b;
}

// Simulates one round of bcrypt's key expansion
static void expand_state(const char* data, int len) {
    // 1. XOR data into P_array
    for (int i = 0; i < 18; ++i) {
        uint32_t d = 0;
        for (int j = 0; j < 4; ++j) {
            d = (d << 8) | data[(i * 4 + j) % len];
        }
        P_array[i] ^= d;
    }

    // 2. Encrypt a zero block and use it to update the whole state
    uint32_t L = 0, R = 0;
    for (int i = 0; i < 18; i += 2) {
        encrypt(&L, &R);
        P_array[i] = L;
        P_array[i+1] = R;
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 256; j += 2) {
            encrypt(&L, &R);
            S_boxes[i * 256 + j] = L;
            S_boxes[i * 256 + j + 1] = R;
        }
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <cost_factor> <seed>\n", argv[0]);
        exit(1);
    }

    cost_factor = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (cost_factor < 4 || cost_factor > 30) {
        fprintf(stderr, "Cost factor must be between 4 and 30.\n");
        exit(1);
    }

    // Allocate memory
    password = (char*)malloc(PW_LEN * sizeof(char));
    salt = (char*)malloc(SALT_LEN * sizeof(char));
    P_array = (uint32_t*)malloc(18 * sizeof(uint32_t));
    S_boxes = (uint32_t*)malloc(1024 * sizeof(uint32_t));
    if (!password || !salt || !P_array || !S_boxes) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate random password and salt
    for (int i = 0; i < PW_LEN; ++i) {
        password[i] = mt_rand() & 0xFF;
    }
    for (int i = 0; i < SALT_LEN; ++i) {
        salt[i] = mt_rand() & 0xFF;
    }
    
    // Initialize Blowfish state with Pi constants
    memcpy(P_array, PI_P, 18 * sizeof(uint32_t));
    memcpy(S_boxes, PI_S, 1024 * sizeof(uint32_t));
}

void run_computation() {
    long rounds = 1L << cost_factor;
    
    // First, mix in the salt and password to set up initial state
    expand_state(salt, SALT_LEN);
    expand_state(password, PW_LEN);

    // Main expensive loop in bcrypt
    for (long i = 0; i < rounds; i++) {
        // In real bcrypt, it alternates hashing the password and the salt.
        // We simulate this by calling expand_state for both.
        expand_state(password, PW_LEN);
        expand_state(salt, SALT_LEN);
    }

    // To prevent dead code elimination, compute a checksum of the final state.
    uint32_t checksum = 0;
    for (int i = 0; i < 18; ++i) {
        checksum ^= P_array[i];
    }
    for (int i = 0; i < 1024; ++i) {
        checksum ^= S_boxes[i];
    }
    final_checksum = checksum; 
}

void cleanup() {
    free(password);
    free(salt);
    free(P_array);
    free(S_boxes);
}

/*******************************************************************************
 * Main Function
 ******************************************************************************/

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%u\n", final_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}