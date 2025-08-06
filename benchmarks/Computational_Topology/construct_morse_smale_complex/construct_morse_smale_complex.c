#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// --- BENCHMARK DATA AND GLOBALS ---

typedef struct {
    int dimension;
    double value; // Simulates a scalar field value on the simplex
    int num_cofaces; // Higher-dimensional neighbors
    int* cofaces; // Indices to coface simplices
} Simplex;

int g_num_simplices;
int g_max_dimension;
Simplex* g_simplicial_complex;
long long g_final_result; // Accumulated result to prevent dead code elimination

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_simplices max_dimension seed\n", argv[0]);
        exit(1);
    }

    g_num_simplices = atoi(argv[1]);
    g_max_dimension = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_num_simplices <= 0 || g_max_dimension <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters. Must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    g_simplicial_complex = (Simplex*)malloc(g_num_simplices * sizeof(Simplex));
    if (!g_simplicial_complex) {
        fprintf(stderr, "FATAL: Memory allocation failed for simplicial complex.\n");
        exit(1);
    }

    for (int i = 0; i < g_num_simplices; ++i) {
        Simplex* s = &g_simplicial_complex[i];
        s->dimension = mt_rand() % (g_max_dimension + 1);
        s->value = (double)mt_rand() / (double)UINT32_MAX;
        
        // Assign a random number of higher-dimensional neighbors (cofaces)
        // A simplex of dimension 'd' can have cofaces of dimension 'd+1'.
        // We simplify this by just having random connections.
        s->num_cofaces = (mt_rand() % (g_max_dimension + 2)) + 1; // at least one coface

        s->cofaces = (int*)malloc(s->num_cofaces * sizeof(int));
        if (!s->cofaces) {
            fprintf(stderr, "FATAL: Memory allocation failed for cofaces.\n");
            exit(1);
        }

        for (int j = 0; j < s->num_cofaces; ++j) {
            // In a real complex, cofaces must have dimension+1. Here, we form a random graph.
            s->cofaces[j] = mt_rand() % g_num_simplices;
        }
    }
}

void run_computation() {
    g_final_result = 0;
    
    // This algorithm simulates finding paths in a Morse-Smale complex.
    // For each simplex, we trace a path of "steepest ascent" by moving
    // to the coface with the highest value, until we reach a local maximum.
    // The total length of all such paths is the result.
    for (int i = 0; i < g_num_simplices; ++i) {
        int current_idx = i;
        long long path_len = 0;

        while (1) {
            Simplex* current_simplex = &g_simplicial_complex[current_idx];
            double max_value = current_simplex->value;
            int best_coface_idx = -1;

            for (int j = 0; j < current_simplex->num_cofaces; ++j) {
                int coface_idx = current_simplex->cofaces[j];
                if (g_simplicial_complex[coface_idx].value > max_value) {
                    max_value = g_simplicial_complex[coface_idx].value;
                    best_coface_idx = coface_idx;
                }
            }

            if (best_coface_idx == -1) {
                // No coface has a higher value, this is a local maximum for this path.
                break;
            }
            
            current_idx = best_coface_idx;
            path_len++;
        }
        g_final_result += path_len;
    }
}

void cleanup() {
    if (g_simplicial_complex) {
        for (int i = 0; i < g_num_simplices; ++i) {
            free(g_simplicial_complex[i].cofaces);
        }
        free(g_simplicial_complex);
    }
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
    printf("%lld\n", g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
