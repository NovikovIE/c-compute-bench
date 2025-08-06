#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- BEGIN: Mersenne Twister (Do Not Modify) ---
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
// --- END: Mersenne Twister ---

// Benchmark data and parameters
typedef struct {
    int num_historical_bars;
    int num_assets;
    int strategy_complexity_score;
    double** price_history; // 2D array: [asset][bar]
} BenchmarkData;

BenchmarkData g_data;
double g_final_result = 0.0;

// Generates a random double between -0.5 and 0.5
double rand_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) - 0.5;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_historical_bars num_assets strategy_complexity_score seed\n", argv[0]);
        exit(1);
    }

    g_data.num_historical_bars = atoi(argv[1]);
    g_data.num_assets = atoi(argv[2]);
    g_data.strategy_complexity_score = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory for price history
    g_data.price_history = (double**)malloc(g_data.num_assets * sizeof(double*));
    if (g_data.price_history == NULL) {
        fprintf(stderr, "Memory allocation failed for price history pointers.\n");
        exit(1);
    }

    // Generate price history for each asset using a simple random walk
    for (int i = 0; i < g_data.num_assets; ++i) {
        g_data.price_history[i] = (double*)malloc(g_data.num_historical_bars * sizeof(double));
        if (g_data.price_history[i] == NULL) {
            fprintf(stderr, "Memory allocation failed for asset %d.\n", i);
            exit(1);
        }

        g_data.price_history[i][0] = 100.0 + (rand_double() * 20.0); // Start price around 100
        for (int j = 1; j < g_data.num_historical_bars; ++j) {
            double move = rand_double(); // Small random move
            g_data.price_history[i][j] = g_data.price_history[i][j-1] + move;
            // Ensure price doesn't go below a certain threshold
            if (g_data.price_history[i][j] < 1.0) {
                 g_data.price_history[i][j] = 1.0;
            }
        }
    }
}

void run_computation() {
    double total_pnl = 0.0;

    for (int i = 0; i < g_data.num_assets; ++i) {
        double asset_pnl = 0.0;
        int position = 0; // -1 for short, 0 for flat, 1 for long

        for (int j = 1; j < g_data.num_historical_bars; ++j) {
            double current_price = g_data.price_history[i][j];
            double prev_price = g_data.price_history[i][j - 1];
            double signal = 0.0;

            // Simulate a computationally intensive strategy calculation
            double temp_calc = 1.0;
            for (int k = 0; k < g_data.strategy_complexity_score; ++k) {
                // This calculation is arbitrary and designed to consume CPU cycles
                // while depending on the input data to prevent over-optimization.
                temp_calc += (prev_price * (k + 1.0)) / (current_price + k + 1.1);
                // Use fmod to keep the number within a manageable range
                temp_calc = fmod(temp_calc, 10000.0);
            }
            signal = temp_calc;

            // Update profit/loss from the previous bar's position
            if (position == 1) { // Was long
                asset_pnl += current_price - prev_price;
            } else if (position == -1) { // Was short
                asset_pnl -= current_price - prev_price;
            }

            // Simple trading logic based on the complex signal calculation
            if (signal > 5000.0 && position <= 0) {
                position = 1; // "Buy"
            } else if (signal < 5000.0 && position >= 0) {
                position = -1; // "Sell"
            }
        }
        total_pnl += asset_pnl;
    }

    g_final_result = total_pnl;
}

void cleanup() {
    for (int i = 0; i < g_data.num_assets; ++i) {
        free(g_data.price_history[i]);
    }
    free(g_data.price_history);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result to stdout to prevent dead code elimination
    printf("%f\n", g_final_result);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}