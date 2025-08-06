#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

/*
 * Mersenne Twister (MT19937)
 * Do Not Modify This Section
 */
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
/* End of Mersenne Twister */


// --- Benchmark Globals ---
long g_signal_length;
int g_window_size;
float *g_input_signal;
float *g_output_signal;
float *g_window_buffer;
double g_result_accumulator;

// Helper comparison function for qsort
int compare_floats(const void *a, const void *b) {
    float fa = *(const float *)a;
    float fb = *(const float *)b;
    if (fa < fb) return -1;
    if (fa > fb) return 1;
    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <signal_length> <window_size> <seed>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    g_signal_length = atol(argv[1]);
    g_window_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_signal_length <= 0 || g_window_size <= 0) {
        fprintf(stderr, "FATAL: signal_length and window_size must be positive.\n");
        exit(EXIT_FAILURE);
    }
    if (g_window_size > g_signal_length) {
        fprintf(stderr, "FATAL: window_size cannot be larger than signal_length.\n");
        exit(EXIT_FAILURE);
    }

    // A median filter requires an odd-sized window to have a unique middle element.
    // If an even size is given, we increment it to the next odd number.
    if (g_window_size % 2 == 0) {
        g_window_size++;
        // We must re-check this condition after possibly incrementing window_size
        if (g_window_size > g_signal_length) {
          fprintf(stderr, "FATAL: window_size (adjusted to be odd) is larger than signal_length.\n");
          exit(EXIT_FAILURE);
        }
    }
    
    mt_seed(seed);

    g_input_signal = (float *)malloc(g_signal_length * sizeof(float));
    g_output_signal = (float *)malloc(g_signal_length * sizeof(float)); // Allocate for max possible output
    g_window_buffer = (float *)malloc(g_window_size * sizeof(float));

    if (!g_input_signal || !g_output_signal || !g_window_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    // Generate a random signal
    for (long i = 0; i < g_signal_length; i++) {
        g_input_signal[i] = (float)mt_rand() / (float)UINT32_MAX;
    }
}

void run_computation() {
    long output_length = g_signal_length - g_window_size + 1;
    int median_index = g_window_size / 2;

    // Apply the median filter
    for (long i = 0; i < output_length; i++) {
        // Copy the current window of the signal into the temporary buffer
        for(int j = 0; j < g_window_size; j++) {
            g_window_buffer[j] = g_input_signal[i + j];
        }

        // Sort the buffer to find the median value
        qsort(g_window_buffer, g_window_size, sizeof(float), compare_floats);

        // The median is the middle element of the sorted window
        g_output_signal[i] = g_window_buffer[median_index];
    }
    
    // Accumulate the results into a single value to prevent dead code elimination
    // and provide a verifiable result.
    g_result_accumulator = 0.0;
    for (long i = 0; i < output_length; i++) {
        g_result_accumulator += g_output_signal[i];
    }
}

void cleanup() {
    free(g_input_signal);
    free(g_output_signal);
    free(g_window_buffer);
    g_input_signal = NULL;
    g_output_signal = NULL;
    g_window_buffer = NULL;
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
    printf("%f\n", g_result_accumulator);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
