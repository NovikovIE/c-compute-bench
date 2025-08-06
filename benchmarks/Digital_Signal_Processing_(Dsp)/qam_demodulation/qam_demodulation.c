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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA AND CONFIGURATION ---
typedef struct {
    double real;
    double imag;
} Complex;

// Global struct to hold all benchmark data
struct {
    long num_samples;
    int constellation_size;

    Complex *received_signal;
    Complex *constellation_map;

    long final_result; // Accumulator to prevent dead code elimination
} g_data;

// --- BENCHMARK FUNCTIONS ---

// Helper to generate a random double in [0, 1)
double mt_rand_double() {
    return (double)mt_rand() / 4294967296.0; // 2^32
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_samples> <constellation_size> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atol(argv[1]);
    g_data.constellation_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    // Validate constellation_size: must be a perfect square
    int grid_dim = (int)round(sqrt(g_data.constellation_size));
    if (grid_dim * grid_dim != g_data.constellation_size) {
        fprintf(stderr, "FATAL: constellation_size must be a perfect square (e.g., 4, 16, 64).\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    g_data.constellation_map = (Complex *)malloc(g_data.constellation_size * sizeof(Complex));
    g_data.received_signal = (Complex *)malloc(g_data.num_samples * sizeof(Complex));

    if (!g_data.constellation_map || !g_data.received_signal) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate constellation map (e.g., for 16-QAM, points are on a 4x4 grid)
    int k = 0;
    for (int i = 0; i < grid_dim; i++) {
        for (int j = 0; j < grid_dim; j++) {
            g_data.constellation_map[k].real = (double)(2 * i - (grid_dim - 1));
            g_data.constellation_map[k].imag = (double)(2 * j - (grid_dim - 1));
            k++;
        }
    }

    // Generate noisy received signal
    // Each sample is a randomly chosen constellation point plus some noise
    for (long i = 0; i < g_data.num_samples; i++) {
        int symbol_index = mt_rand() % g_data.constellation_size;
        Complex ideal_point = g_data.constellation_map[symbol_index];

        // Add Gaussian-like noise (Box-Muller could be used, but this is simpler for a CPU benchmark)
        double noise_amplitude = 0.5;
        double u1 = mt_rand_double();
        double u2 = mt_rand_double();
        double noise_real = noise_amplitude * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        double noise_imag = noise_amplitude * sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);

        g_data.received_signal[i].real = ideal_point.real + noise_real;
        g_data.received_signal[i].imag = ideal_point.imag + noise_imag;
    }

    g_data.final_result = 0;
}

void run_computation() {
    long accumulated_symbols = 0;

    for (long i = 0; i < g_data.num_samples; i++) {
        Complex received_sample = g_data.received_signal[i];
        double min_dist_sq = -1.0;
        int best_symbol_index = -1;

        // Find the closest constellation point for the current sample
        for (int j = 0; j < g_data.constellation_size; j++) {
            Complex const_point = g_data.constellation_map[j];
            double dr = received_sample.real - const_point.real;
            double di = received_sample.imag - const_point.imag;
            double dist_sq = dr * dr + di * di;

            if (best_symbol_index == -1 || dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                best_symbol_index = j;
            }
        }
        accumulated_symbols += best_symbol_index;
    }

    g_data.final_result = accumulated_symbols;
}

void cleanup() {
    free(g_data.constellation_map);
    free(g_data.received_signal);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%ld\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
