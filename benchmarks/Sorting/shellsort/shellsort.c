#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Global variables for benchmark data ---
long num_elements;
int gap_sequence_type;
int *data_array;
int *gaps;
int num_gaps;
int final_result;

// --- Mersenne Twister (Provided) ---
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

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_elements> <gap_sequence_type> <seed>\n", argv[0]);
        fprintf(stderr, "  gap_sequence_type: 0=Shell, 1=Ciura, 2=Sedgewick(1,8,23,..)\n");
        exit(1);
    }

    num_elements = atol(argv[1]);
    gap_sequence_type = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_elements <= 0) {
        fprintf(stderr, "Number of elements must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    data_array = (int*)malloc(num_elements * sizeof(int));
    if (data_array == NULL) {
        fprintf(stderr, "Memory allocation failed for data_array\n");
        exit(1);
    }

    for (long i = 0; i < num_elements; i++) {
        data_array[i] = (int)(mt_rand()); // Use full random range
    }

    switch (gap_sequence_type) {
        case 0: { // Shell's original sequence: N/2, N/4, ..., 1
            int count = 0;
            for (long g = num_elements / 2; g > 0; g /= 2) {
                count++;
            }
            gaps = (int*)malloc(count * sizeof(int));
            if (gaps == NULL) { fprintf(stderr, "Gaps allocation failed\n"); exit(1); }
            num_gaps = count;
            int i = 0;
            for (long g = num_elements / 2; g > 0; g /= 2) {
                gaps[i++] = g;
            }
            break;
        }
        case 1: { // Ciura's sequence (empirically derived, good performance)
            const int ciura_seq[] = {1, 4, 10, 23, 57, 132, 301, 701, 1750, 3905, 8929, 19876, 44284, 98517};
            int max_gaps = sizeof(ciura_seq)/sizeof(int);
            int count = 0;
            for (int i = 0; i < max_gaps; i++) {
                if (ciura_seq[i] < num_elements) {
                    count++;
                } else {
                    break;
                }
            }
            if (count == 0 && num_elements > 1) { count = 1; } // At least gap 1
            if (count == 0 && num_elements <= 1) { num_gaps = 0; break; }

            gaps = (int*)malloc(count * sizeof(int));
            if (gaps == NULL) { fprintf(stderr, "Gaps allocation failed\n"); exit(1); }
            num_gaps = count;
            for (int i = 0; i < num_gaps; i++) {
                gaps[i] = ciura_seq[num_gaps - 1 - i];
            }
            break;
        }
        case 2: { // Sedgewick's sequence: 1, 8, 23, 77, ...
            int temp_gaps[32];
            int count = 0;
            for (int k = 1; k < 32; k++) {
                long gap = (1L << (2 * k)) + 3 * (1L << (k - 1)) + 1;
                if (gap >= num_elements || gap < 0) break; // Check for overflow
                if (count < 32) temp_gaps[count++] = (int)gap;
            }
            num_gaps = count + 1; // +1 for the gap '1'
            gaps = (int*)malloc(num_gaps * sizeof(int));
            if (gaps == NULL) { fprintf(stderr, "Gaps allocation failed\n"); exit(1); }
            for(int i = 0; i < count; i++) {
                gaps[i] = temp_gaps[count - 1 - i];
            }
            gaps[count] = 1;
            break;
        }
        default:
            fprintf(stderr, "Invalid gap sequence type\n");
            exit(1);
    }
}

void run_computation() {
    for (int i = 0; i < num_gaps; i++) {
        int gap = gaps[i];
        for (long j = gap; j < num_elements; j++) {
            int temp = data_array[j];
            long k;
            for (k = j; k >= gap && data_array[k - gap] > temp; k -= gap) {
                data_array[k] = data_array[k - gap];
            }
            data_array[k] = temp;
        }
    }
    // Result to prevent dead code elimination
    if (num_elements > 0) {
        final_result = data_array[num_elements / 2];
    } else {
        final_result = 0;
    }
}

void cleanup() {
    if (data_array != NULL) free(data_array);
    if (gaps != NULL) free(gaps);
    data_array = NULL;
    gaps = NULL;
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
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
