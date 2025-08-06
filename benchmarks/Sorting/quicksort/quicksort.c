#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// --- Mersenne Twister (MT19937) --- Do Not Modify ---
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
// --- End of Mersenne Twister ---

// --- Global Data ---
typedef int data_t;

static data_t *data_array;
static int num_elements;
static int final_result; // For preventing dead code elimination

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_elements> <initial_data_distribution> <seed>\n", argv[0]);
        exit(1);
    }

    num_elements = atoi(argv[1]);
    const char *distribution = argv[2];
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (num_elements <= 0) {
        fprintf(stderr, "FATAL: num_elements must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    data_array = (data_t *)malloc(num_elements * sizeof(data_t));
    if (data_array == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for data_array.\n");
        exit(1);
    }

    if (strcmp(distribution, "random") == 0) {
        for (int i = 0; i < num_elements; ++i) {
            data_array[i] = (data_t)(mt_rand() % 1000000);
        }
    } else if (strcmp(distribution, "sorted") == 0) {
        for (int i = 0; i < num_elements; ++i) {
            data_array[i] = (data_t)i;
        }
    } else if (strcmp(distribution, "reverse_sorted") == 0) {
        for (int i = 0; i < num_elements; ++i) {
            data_array[i] = (data_t)(num_elements - 1 - i);
        }
    } else if (strcmp(distribution, "nearly_sorted") == 0) {
        for (int i = 0; i < num_elements; ++i) {
            data_array[i] = (data_t)i;
        }
        // Perform a few swaps
        int swaps = num_elements / 100;
        for (int i = 0; i < swaps; ++i) {
            int idx1 = mt_rand() % num_elements;
            int idx2 = mt_rand() % num_elements;
            data_t temp = data_array[idx1];
            data_array[idx1] = data_array[idx2];
            data_array[idx2] = temp;
        }
    } else {
        fprintf(stderr, "FATAL: Unknown data distribution '%s'\n", distribution);
        exit(1);
    }
}

void swap(data_t* a, data_t* b) {
    data_t t = *a;
    *a = *b;
    *b = t;
}

// Lomuto partition scheme
int partition(data_t arr[], int low, int high) {
    data_t pivot = arr[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

// Iterative quicksort to prevent stack overflow on worst-case inputs
void quicksort_iterative(data_t arr[], int low, int high) {
    if (low >= high) return;

    // Create an auxiliary stack
    int stack_size = high - low + 1;
    int* stack = (int*)malloc(stack_size * sizeof(int));
    if (stack == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for quicksort stack.\n");
        exit(1);
    }

    // initialize top of stack
    int top = -1;

    // push initial values of low and high to stack
    stack[++top] = low;
    stack[++top] = high;

    // Keep popping from stack while is not empty
    while (top >= 0) {
        // Pop high and low
        high = stack[top--];
        low = stack[top--];

        // Set pivot element at its correct position in sorted array
        int p = partition(arr, low, high);

        // If there are elements on left side of pivot, then push left
        // side to stack
        if (p - 1 > low) {
            stack[++top] = low;
            stack[++top] = p - 1;
        }

        // If there are elements on right side of pivot, then push right
        // side to stack
        if (p + 1 < high) {
            stack[++top] = p + 1;
            stack[++top] = high;
        }
    }
    free(stack);
}

void run_computation() {
    quicksort_iterative(data_array, 0, num_elements - 1);

    // Calculate a checksum to prevent dead code elimination
    final_result = 0;
    int step = num_elements > 100 ? num_elements / 100 : 1;
    for (int i = 0; i < num_elements; i += step) {
      final_result += data_array[i];
    }
}

void cleanup() {
    free(data_array);
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
