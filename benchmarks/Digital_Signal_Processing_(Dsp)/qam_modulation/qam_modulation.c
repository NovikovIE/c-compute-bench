#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator --- DO NOT MODIFY ---
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

// Benchmark state
typedef struct {
    long num_symbols;
    int constellation_size;
    int levels_per_dim;
    unsigned int *input_symbols;
    double *modulated_I;
    double *modulated_Q;
    double total_power;
} BenchmarkData;

BenchmarkData g_data;


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_symbols constellation_size seed\n", argv[0]);
        exit(1);
    }

    g_data.num_symbols = atol(argv[1]);
    g_data.constellation_size = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.num_symbols <= 0 || g_data.constellation_size <= 0) {
        fprintf(stderr, "Error: num_symbols and constellation_size must be positive.\n");
        exit(1);
    }

    g_data.levels_per_dim = (int)round(sqrt(g_data.constellation_size));
    if (g_data.levels_per_dim * g_data.levels_per_dim != g_data.constellation_size || g_data.constellation_size < 4) {
        fprintf(stderr, "Error: constellation_size must be a perfect square >= 4.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.input_symbols = (unsigned int *)malloc(g_data.num_symbols * sizeof(unsigned int));
    g_data.modulated_I = (double *)malloc(g_data.num_symbols * sizeof(double));
    g_data.modulated_Q = (double *)malloc(g_data.num_symbols * sizeof(double));

    if (!g_data.input_symbols || !g_data.modulated_I || !g_data.modulated_Q) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    for (long i = 0; i < g_data.num_symbols; ++i) {
        g_data.input_symbols[i] = mt_rand() % g_data.constellation_size;
    }

    g_data.total_power = 0.0;
}

void run_computation() {
    double power_sum = 0.0;
    int levels = g_data.levels_per_dim;
    long n_symbols = g_data.num_symbols;
    unsigned int* symbols = g_data.input_symbols;
    double* i_out = g_data.modulated_I;
    double* q_out = g_data.modulated_Q;

    for (long i = 0; i < n_symbols; ++i) {
        unsigned int symbol = symbols[i];
        
        // Decompose symbol index into I and Q indices
        int i_index = symbol % levels;
        int q_index = symbol / levels;
        
        // Map index to amplitude level, e.g., for 4 levels (0,1,2,3) -> (-3,-1,1,3)
        double i_val = (double)(2 * i_index - (levels - 1));
        double q_val = (double)(2 * q_index - (levels - 1));
        
        i_out[i] = i_val;
        q_out[i] = q_val;
        
        // Accumulate power to prevent dead code elimination
        power_sum += i_val * i_val + q_val * q_val;
    }

    g_data.total_power = power_sum;
}

void cleanup() {
    free(g_data.input_symbols);
    free(g_data.modulated_I);
    free(g_data.modulated_Q);
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
    printf("%.2f\n", g_data.total_power);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
