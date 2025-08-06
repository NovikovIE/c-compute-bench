#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA ---
// Parameters
long long S; // The number to compute sqrt(S) for. Use long long for S-m*m
int MAX_TERMS;

// Data structures
int* terms; // To store the continued fraction terms
long long result_sum; // To accumulate a result and prevent dead-code elimination

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <number> <max_terms> <seed>\n", argv[0]);
        exit(1);
    }

    S = atoll(argv[1]);
    MAX_TERMS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed); // Seed the generator, though not used in the core logic

    // Basic validation
    if (S <= 0 || MAX_TERMS <= 0) {
        fprintf(stderr, "FATAL: Number and max_terms must be positive.\n");
        exit(1);
    }
    
    // Check if S is a perfect square. The algorithm terminates for perfect squares.
    long long root = round(sqrt((long double)S));
    if (root * root == S) {
         // The benchmark is designed for irrational roots, but we can allow this.
         // No warning is printed to keep stderr clean for the timer.
    }

    terms = (int*)malloc((size_t)MAX_TERMS * sizeof(int));
    if (terms == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for terms array.\n");
        exit(1);
    }

    result_sum = 0;
}

void run_computation() {
    // This function implements the PQa algorithm (sometimes attributed to Lagrange)
    // for finding the continued fraction of a quadratic irrational, specifically sqrt(S).
    // All calculations can be done with integer arithmetic, requiring long long
    // to prevent overflow for intermediate products.
    
    long long m = 0;
    long long d = 1;
    long long a0 = floor(sqrt((long double)S));
    long long a = a0;

    for (int i = 0; i < MAX_TERMS; ++i) {
        terms[i] = (int)a;

        // Calculate next state
        // m_next = d * a - m
        m = d * a - m;

        // d_next = (S - m_next^2) / d
        d = (S - m * m) / d;

        // If d becomes 0, it means S was a perfect square and the expansion terminates
        if (d == 0) {
            // Pad the rest of the terms with 0
            for (int j = i + 1; j < MAX_TERMS; ++j) {
                terms[j] = 0;
            }
            break; 
        }

        // a_next = floor((a0 + m_next) / d_next)
        a = floor((double)(a0 + m) / (double)d);
    }
    
    // Accumulate a result to ensure the computation is not optimized away
    for (int i = 0; i < MAX_TERMS; ++i) {
        result_sum += terms[i];
    }
}

void cleanup() {
    free(terms);
    terms = NULL;
}

// --- MAIN ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%lld\n", result_sum);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
