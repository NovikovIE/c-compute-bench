#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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
typedef struct {
    int signal_length;
    int interpolation_factor;
    int decimation_factor;
    int num_taps_per_phase;
    long long output_signal_length;

    float *input_signal;
    float *output_signal;
    float *filter_coeffs; // Flattened 2D array: L x num_taps_per_phase

    float final_result; // For verification and preventing DCE
} BenchmarkData;

BenchmarkData *g_data = NULL;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

// --- Benchmark Implementation ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s signal_length interp_factor decim_factor taps_per_phase seed\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("Failed to allocate memory for BenchmarkData");
        exit(1);
    }

    // Parse arguments
    g_data->signal_length = atoi(argv[1]);
    g_data->interpolation_factor = atoi(argv[2]);
    g_data->decimation_factor = atoi(argv[3]);
    g_data->num_taps_per_phase = atoi(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 10);
    mt_seed(seed);

    // Validate parameters
    if (g_data->signal_length <= 0 || g_data->interpolation_factor <= 0 ||
        g_data->decimation_factor <= 0 || g_data->num_taps_per_phase <= 0) {
        fprintf(stderr, "FATAL: All parameters must be positive integers.\n");
        exit(1);
    }
    if (g_data->num_taps_per_phase > g_data->signal_length) {
        fprintf(stderr, "FATAL: num_taps_per_phase cannot be greater than signal_length.\n");
        exit(1);
    }

    // Calculate output length and allocate memory
    g_data->output_signal_length = ((long long)g_data->signal_length * g_data->interpolation_factor) / g_data->decimation_factor;
    
    g_data->input_signal = (float*)malloc((size_t)g_data->signal_length * sizeof(float));
    g_data->output_signal = (float*)malloc((size_t)g_data->output_signal_length * sizeof(float));
    size_t total_filter_taps = (size_t)g_data->interpolation_factor * g_data->num_taps_per_phase;
    g_data->filter_coeffs = (float*)malloc(total_filter_taps * sizeof(float));

    if (!g_data->input_signal || !g_data->output_signal || !g_data->filter_coeffs) {
        perror("Failed to allocate memory for data arrays");
        cleanup(); // Use cleanup to free any partially allocated memory
        exit(1);
    }

    // Initialize input signal and filter coefficients with random values
    for (int i = 0; i < g_data->signal_length; i++) {
        g_data->input_signal[i] = (float)mt_rand() / (float)UINT32_MAX;
    }
    for (size_t i = 0; i < total_filter_taps; i++) {
        g_data->filter_coeffs[i] = (float)mt_rand() / (float)UINT32_MAX;
    }

    g_data->final_result = 0.0f;
}

void run_computation() {
    const int L = g_data->interpolation_factor;
    const int M = g_data->decimation_factor;
    const int N_taps = g_data->num_taps_per_phase;
    
    float total_sum = 0.0f;

    // Iterate over each output sample
    for (long long k = 0; k < g_data->output_signal_length; k++) {
        // Calculate the corresponding time on the input signal's grid.
        // input_time = k * M / L
        long long input_base_index = (k * M) / L;
        
        // The phase of the filter to use is determined by the fractional part of the time,
        // which simplifies to (k * M) % L.
        int phase_index = (int)((k * M) % L);
        
        float accum = 0.0f;
        
        // Get a pointer to the start of the required filter phase coefficients.
        const float *phase_filter = g_data->filter_coeffs + (size_t)phase_index * N_taps;
        
        // Perform the convolution (FIR filtering) for this output sample.
        // The filter is applied "backwards" in time relative to the input signal.
        for (int i = 0; i < N_taps; i++) {
            long long current_input_index = input_base_index - i;
            
            // Handle boundary conditions: assume signal is zero before index 0.
            if (current_input_index >= 0 && current_input_index < g_data->signal_length) {
                accum += g_data->input_signal[current_input_index] * phase_filter[i];
            }
        }
        
        g_data->output_signal[k] = accum;
        total_sum += accum;
    }
    
    g_data->final_result = total_sum;
}

void cleanup() {
    if (g_data) {
        free(g_data->input_signal);
        g_data->input_signal = NULL;
        free(g_data->output_signal);
        g_data->output_signal = NULL;
        free(g_data->filter_coeffs);
        g_data->filter_coeffs = NULL;
        free(g_data);
        g_data = NULL;
    }
}

// --- Main and Timing ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", g_data->final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
