#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// --- Benchmark Globals and Functions ---

// A struct to hold all benchmark state, making it clear what is shared.
struct BenchmarkState {
    int lz77_window_size;
    unsigned char* input_data;
    size_t input_size_bytes;
    long long total_match_length;
};

static struct BenchmarkState state;

// Parses arguments, allocates memory, and generates compressible random data.
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <lz77_window_size> <seed>\n", argv[0]);
        exit(1);
    }

    // FIX: Use atof to correctly parse the floating-point MB value from arguments.
    // The original code used atol, which would truncate fractional parts (e.g., 1.04 -> 1).
    double input_size_mb = atof(argv[1]);
    state.lz77_window_size = atoi(argv[2]);
    uint32_t seed = (uint32_t)atol(argv[3]);

    mt_seed(seed);

    state.input_size_bytes = (size_t)(input_size_mb * 1024 * 1024);
    state.input_data = (unsigned char*)malloc(state.input_size_bytes);
    if (state.input_data == NULL) {
        fprintf(stderr, "Failed to allocate memory for input data.\n");
        exit(1);
    }
    
    state.total_match_length = 0;

    // Generate semi-compressible data. Pure random data is not compressible.
    // We introduce repetitions by copying blocks of already generated data.
    const size_t INITIAL_RANDOM_FILL = 65536;
    const int COPY_PROBABILITY = 50; // 50% chance to copy a block
    const int MIN_COPY_LEN = 3;
    const int MAX_COPY_LEN = 258; // A common max length in DEFLATE

    // Fill an initial portion with purely random bytes
    size_t i = 0;
    if (INITIAL_RANDOM_FILL < state.input_size_bytes) {
        for (; i < INITIAL_RANDOM_FILL; ++i) {
            state.input_data[i] = mt_rand() & 0xFF;
        }
    } else { // Handle small input sizes
        for (; i < state.input_size_bytes; ++i) {
            state.input_data[i] = mt_rand() & 0xFF;
        }
        return;
    }

    // Fill the rest of the data with a mix of random bytes and copied blocks
    while (i < state.input_size_bytes) {
        if ((mt_rand() % 100) < COPY_PROBABILITY && i > MIN_COPY_LEN) {
            size_t copy_len = MIN_COPY_LEN + (mt_rand() % (MAX_COPY_LEN - MIN_COPY_LEN + 1));
            if (i + copy_len > state.input_size_bytes) {
                copy_len = state.input_size_bytes - i;
            }

            size_t copy_dist = 1 + (mt_rand() % (i - MIN_COPY_LEN));
            size_t copy_from = i - copy_dist;
            
            // FIX: Use memmove instead of memcpy. The AddressSanitizer error
            // "memcpy-param-overlap" indicates that the source and destination
            // memory regions can overlap. memcpy has undefined behavior in this case,
            // while memmove is designed to handle it correctly.
            memmove(&state.input_data[i], &state.input_data[copy_from], copy_len);
            i += copy_len;
        } else {
            state.input_data[i] = mt_rand() & 0xFF;
            i++;
        }
    }
}

// Simulates the core computational work of an LZ77 compressor.
void run_computation() {
    const int MIN_MATCH_LEN = 3; // Minimum length to be considered a "match"

    for (size_t pos = 0; pos < state.input_size_bytes; ) {
        size_t search_start = (pos > (size_t)state.lz77_window_size) ? (pos - state.lz77_window_size) : 0;
        int best_len = 0;

        // Search for the longest match in the sliding window. This is the main workload.
        for (size_t i = search_start; i < pos; ++i) {
            int current_len = 0;
            while (pos + current_len < state.input_size_bytes &&
                   state.input_data[i + current_len] == state.input_data[pos + current_len]) {
                current_len++;
            }

            if (current_len > best_len) {
                best_len = current_len;
            }
        }

        // If a usable match is found, accumulate its length and skip the cursor ahead.
        // This simulates encoding a (distance, length) pair.
        if (best_len >= MIN_MATCH_LEN) {
            state.total_match_length += best_len;
            pos += best_len;
        } else {
            // No useful match. Simulate encoding a literal and move to the next byte.
            pos++;
        }
    }
}

// Frees all memory allocated in setup_benchmark.
void cleanup() {
    if (state.input_data) {
        free(state.input_data);
        state.input_data = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start_time, end_time;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start_time);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end_time);

    cleanup();

    double time_taken = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_nsec - start_time.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout.
    // This value is an aggregate of the work done, preventing dead code elimination.
    printf("%lld\n", state.total_match_length);

    // Print the timing information to stderr.
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
