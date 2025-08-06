#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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
// --- End of Mersenne Twister ---

// A small buffer to hold quaternions, designed to fit in cache
#define BUFFER_SIZE 8192

// Quaternion data structure
typedef struct {
    double w, x, y, z;
} Quaternion;

// Global variables for benchmark data
long num_multiplications;
Quaternion *q1_array;
Quaternion *q2_array;
double final_accumulator;


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_multiplications> <seed>\n", argv[0]);
        exit(1);
    }

    num_multiplications = atol(argv[1]);
    uint32_t seed = atoi(argv[2]);

    mt_seed(seed);

    q1_array = (Quaternion *)malloc(BUFFER_SIZE * sizeof(Quaternion));
    q2_array = (Quaternion *)malloc(BUFFER_SIZE * sizeof(Quaternion));

    if (!q1_array || !q2_array) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < BUFFER_SIZE; i++) {
        // Generate random doubles between -1.0 and 1.0
        q1_array[i].w = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        q1_array[i].x = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        q1_array[i].y = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        q1_array[i].z = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;

        q2_array[i].w = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        q2_array[i].x = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        q2_array[i].y = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
        q2_array[i].z = (mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
    }
}

void run_computation() {
    final_accumulator = 0.0;
    Quaternion qa, qb, result;

    for (long i = 0; i < num_multiplications; i++) {
        // Use cyclic access to the pre-filled buffers.
        // Bitwise AND is faster than modulo for powers of 2.
        int index = i & (BUFFER_SIZE - 1);
        qa = q1_array[index];
        qb = q2_array[index];

        // Quaternion multiplication: result = qa * qb
        result.w = qa.w * qb.w - qa.x * qb.x - qa.y * qb.y - qa.z * qb.z;
        result.x = qa.w * qb.x + qa.x * qb.w + qa.y * qb.z - qa.z * qb.y;
        result.y = qa.w * qb.y - qa.x * qb.z + qa.y * qb.w + qa.z * qb.x;
        result.z = qa.w * qb.z + qa.x * qb.y - qa.y * qb.x + qa.z * qb.w;
        
        // Accumulate one component of the result to prevent dead code elimination
        final_accumulator += result.w;
    }
}

void cleanup() {
    free(q1_array);
    free(q2_array);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final accumulated result to stdout
    printf("%f\n", final_accumulator);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
