#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Benchmark Globals ---
int NUM_DIGITS;
int *pi_digits; // To store the resulting digits of Pi
long long final_result_sum; // Sum of digits to prevent dead-code elimination

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_digits> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_DIGITS = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (NUM_DIGITS <= 0) {
        fprintf(stderr, "Number of digits must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate enough space for the digits, plus a small buffer for the spigot algorithm's nature
    pi_digits = (int *)malloc((NUM_DIGITS + 20) * sizeof(int));
    if (pi_digits == NULL) {
        fprintf(stderr, "Failed to allocate memory for pi digits.\n");
        exit(1);
    }

    final_result_sum = 0;
}

void run_computation() {
    // Implementation of the unbounded spigot algorithm for Pi by J. Gibbons.
    // This is computationally intensive and demonstrates arbitrary precision arithmetic.
    int len = (NUM_DIGITS * 10 / 3) + 2;
    int* a = (int*)malloc(len * sizeof(int));
    if (a == NULL) {
        fprintf(stderr, "Failed to allocate memory for state array.\n");
        exit(1);
    }

    for (int i = 0; i < len; i++) {
        a[i] = 2;
    }

    int digit_idx = 0;
    int nines = 0;
    int predigit = 0;

    for (int j = 0; j < NUM_DIGITS; j++) {
        int q = 0;
        for (int i = len - 1; i > 0; i--) {
            int x = 10 * a[i] + q * (i); // Note: Original formula is q * (i+1), but i is correct for 0-based array
            a[i] = x % (2 * i + 1);
            q = x / (2 * i + 1);
        }

        int x = 10 * a[0] + q; 
        a[0] = x % 10;
        int d = x / 10;

        if (digit_idx >= NUM_DIGITS + 20) break; // Safety break

        if (d < 9) {
            pi_digits[digit_idx++] = predigit;
            for (int k = 0; k < nines; k++) {
                 if (digit_idx < NUM_DIGITS + 20) pi_digits[digit_idx++] = 9;
            }
            predigit = d;
            nines = 0;
        } else if (d == 9) {
            nines++;
        } else { // d == 10
            pi_digits[digit_idx++] = predigit + 1;
            for (int k = 0; k < nines; k++) {
                 if (digit_idx < NUM_DIGITS + 20) pi_digits[digit_idx++] = 0;
            }
            predigit = 0;
            nines = 0;
        }
    }
    
    free(a);

    // Sum the generated digits to produce a final result
    // This prevents the compiler from optimizing away the main computation.
    for (int i = 0; i < NUM_DIGITS; ++i) {
        final_result_sum += pi_digits[i];
    }
}

void cleanup() {
    free(pi_digits);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%lld\n", final_result_sum);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
