#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA STRUCTURES ---
typedef struct {
    int id;
    int dimension;
    double birth_time;
} Simplex;

typedef struct {
    int num_simplices_in_filtration;
    int max_dimension;
    Simplex** filtration;
    long* low_indices;
    unsigned long final_result;
} BenchmarkData;

BenchmarkData* g_data;

// Comparison function for qsort
int compare_simplices(const void* a, const void* b) {
    Simplex* s1 = *(Simplex**)a;
    Simplex* s2 = *(Simplex**)b;
    if (s1->birth_time < s2->birth_time) return -1;
    if (s1->birth_time > s2->birth_time) return 1;
    // Tie-break by dimension for a valid filtration
    if (s1->dimension < s2->dimension) return -1;
    if (s1->dimension > s2->dimension) return 1;
    return 0;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_simplices_in_filtration> <max_dimension> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData*)malloc(sizeof(BenchmarkData));
    if (!g_data) {
        perror("Failed to allocate memory for benchmark data");
        exit(1);
    }

    g_data->num_simplices_in_filtration = atoi(argv[1]);
    g_data->max_dimension = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);

    mt_seed(seed);

    g_data->filtration = (Simplex**)malloc(g_data->num_simplices_in_filtration * sizeof(Simplex*));
    g_data->low_indices = (long*)malloc(g_data->num_simplices_in_filtration * sizeof(long));
    if (!g_data->filtration || !g_data->low_indices) {
        perror("Failed to allocate memory for filtration data");
        exit(1);
    }

    for (int i = 0; i < g_data->num_simplices_in_filtration; ++i) {
        g_data->filtration[i] = (Simplex*)malloc(sizeof(Simplex));
        if (!g_data->filtration[i]) {
            perror("Failed to allocate memory for a simplex");
            exit(1);
        }
        g_data->filtration[i]->id = i;
        g_data->filtration[i]->dimension = mt_rand() % (g_data->max_dimension + 1);
        g_data->filtration[i]->birth_time = (double)mt_rand() / (double)UINT32_MAX;
    }

    // Sort the simplices by birth time to create the filtration
    qsort(g_data->filtration, g_data->num_simplices_in_filtration, sizeof(Simplex*), compare_simplices);
}

void run_computation() {
    // This computation simulates the core loops of a persistent homology algorithm.
    // It's a two-pass process.

    // Pass 1: For each simplex, find its youngest facet in the filtration order.
    // This simulates computing the 'low(j)' entry for each column j in the boundary matrix.
    for (int j = 0; j < g_data->num_simplices_in_filtration; ++j) {
        Simplex* current_simplex = g_data->filtration[j];
        long pivot = -1;
        if (current_simplex->dimension > 0) {
            // Search backwards for the most recently added facet.
            for (int i = j - 1; i >= 0; --i) {
                if (g_data->filtration[i]->dimension == current_simplex->dimension - 1) {
                    pivot = i;
                    break; // Found the youngest facet
                }
            }
        }
        g_data->low_indices[j] = pivot;
    }

    // Pass 2: Simulate the matrix reduction. When two columns j and k have the
    // same low entry (low[j] == low[k]), a column addition is performed.
    // We simulate this with a computationally non-trivial accumulator.
    unsigned long work_accumulator = 0;
    for (int j = 1; j < g_data->num_simplices_in_filtration; ++j) {
        if (g_data->low_indices[j] != -1) {
            for (int k = 0; k < j; ++k) {
                if (g_data->low_indices[k] == g_data->low_indices[j]) {
                    // Collision found -> simulate work of column addition
                    work_accumulator += (unsigned long)(g_data->filtration[j]->id ^ g_data->filtration[k]->id);
                    // Add more dependent computation to prevent optimization
                    work_accumulator = (work_accumulator * 31) + g_data->filtration[j]->dimension;
                }
            }
        }
    }
    g_data->final_result = work_accumulator;
}

void cleanup() {
    for (int i = 0; i < g_data->num_simplices_in_filtration; ++i) {
        free(g_data->filtration[i]);
    }
    free(g_data->filtration);
    free(g_data->low_indices);
    free(g_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent dead code elimination
    printf("%lu\n", g_data->final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
