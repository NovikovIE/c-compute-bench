#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- BEGIN MERSENNE TWISTER (MT19937) ---
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
// --- END MERSENNE TWISTER (MT19937) ---

// --- BENCHMARK DATA AND PARAMETERS ---
static int SIGNAL_LENGTH;
static int NUM_FILTER_TAPS;
static double LEARNING_RATE;

static double *input_signal;
static double *desired_signal;
static double *filter_weights;
static double *error_signal; 

static double final_result;

// --- HELPER FUNCTIONS ---
double random_double_neg1_to_1() {
    return ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
}

// --- BENCHMARK FUNCTIONS ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s signal_length num_filter_taps learning_rate seed\n", argv[0]);
        exit(1);
    }

    SIGNAL_LENGTH = atoi(argv[1]);
    NUM_FILTER_TAPS = atoi(argv[2]);
    LEARNING_RATE = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    input_signal = (double *)malloc(SIGNAL_LENGTH * sizeof(double));
    desired_signal = (double *)malloc(SIGNAL_LENGTH * sizeof(double));
    filter_weights = (double *)malloc(NUM_FILTER_TAPS * sizeof(double));
    error_signal = (double *)malloc(SIGNAL_LENGTH * sizeof(double)); 

    if (!input_signal || !desired_signal || !filter_weights || !error_signal) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < SIGNAL_LENGTH; ++i) {
        input_signal[i] = random_double_neg1_to_1();
        desired_signal[i] = random_double_neg1_to_1(); // In a real scenario, this would be correlated to the input
        error_signal[i] = 0.0;
    }

    for (int i = 0; i < NUM_FILTER_TAPS; ++i) {
        filter_weights[i] = 0.0;
    }

    final_result = 0.0;
}

void run_computation() {
    // LMS aaptive filter algorithm
    for (int n = NUM_FILTER_TAPS - 1; n < SIGNAL_LENGTH; ++n) {
        // 1. Calculate filter output y[n]
        double y_n = 0.0;
        for (int k = 0; k < NUM_FILTER_TAPS; ++k) {
            y_n += filter_weights[k] * input_signal[n - k];
        }

        // 2. Calculate error e[n]
        double e_n = desired_signal[n] - y_n;
        error_signal[n] = e_n;

        // 3. Update filter weights w[k]
        for (int k = 0; k < NUM_FILTER_TAPS; ++k) {
            filter_weights[k] += LEARNING_RATE * e_n * input_signal[n - k];
        }
    }

    // Accumulate a final result to prevent dead code elimination
    double error_sum = 0.0;
    for (int i = 0; i < SIGNAL_LENGTH; ++i) {
        error_sum += fabs(error_signal[i]);
    }
    final_result = error_sum;
}

void cleanup() {
    free(input_signal);
    free(desired_signal);
    free(filter_weights);
    free(error_signal);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
