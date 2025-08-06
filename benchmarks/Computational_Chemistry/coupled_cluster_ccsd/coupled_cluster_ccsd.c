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

// Benchmark data structure
typedef struct {
    int num_occupied_orbitals;
    int num_virtual_orbitals;
    double**** v_vvvv; // Two-electron integrals over virtual orbitals <ab||cd>
    double**** t2_oovv; // T2 amplitudes t_ij^ab
    double**** t2_new_oovv; // Updated T2 amplitudes
    double final_result;
} BenchmarkData;

static BenchmarkData B;

// Helper to allocate a 4D array of doubles
double**** alloc_4d_double(int d1, int d2, int d3, int d4) {
    double**** array = (double****)malloc(d1 * sizeof(double***));
    if (!array) return NULL;

    for (int i = 0; i < d1; i++) {
        array[i] = (double***)malloc(d2 * sizeof(double**));
        if (!array[i]) return NULL; 
        for (int j = 0; j < d2; j++) {
            array[i][j] = (double**)malloc(d3 * sizeof(double*));
            if (!array[i][j]) return NULL;
            for (int k = 0; k < d3; k++) {
                array[i][j][k] = (double*)malloc(d4 * sizeof(double));
                if (!array[i][j][k]) return NULL;
            }
        }
    }
    return array;
}

// Helper to free a 4D array of doubles
void free_4d_double(double**** array, int d1, int d2, int d3) {
    if (!array) return;
    for (int i = 0; i < d1; i++) {
        if (array[i]) {
            for (int j = 0; j < d2; j++) {
                if (array[i][j]) {
                    for (int k = 0; k < d3; k++) {
                        free(array[i][j][k]);
                    }
                    free(array[i][j]);
                }
            }
            free(array[i]);
        }
    }
    free(array);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_occupied_orbitals> <num_virtual_orbitals> <seed>\n", argv[0]);
        exit(1);
    }

    B.num_occupied_orbitals = atoi(argv[1]);
    B.num_virtual_orbitals = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    int o = B.num_occupied_orbitals;
    int v = B.num_virtual_orbitals;

    // Allocate tensors
    B.v_vvvv = alloc_4d_double(v, v, v, v);
    B.t2_oovv = alloc_4d_double(o, o, v, v);
    B.t2_new_oovv = alloc_4d_double(o, o, v, v);
    
    if (!B.v_vvvv || !B.t2_oovv || !B.t2_new_oovv) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize tensors with random data
    for (int a = 0; a < v; a++) {
        for (int b = 0; b < v; b++) {
            for (int c = 0; c < v; c++) {
                for (int d = 0; d < v; d++) {
                    B.v_vvvv[a][b][c][d] = ((double)mt_rand() / (double)UINT32_MAX - 0.5) * 0.1;
                }
            }
        }
    }

    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    B.t2_oovv[i][j][a][b] = ((double)mt_rand() / (double)UINT32_MAX - 0.5) * 0.1;
                    B.t2_new_oovv[i][j][a][b] = 0.0;
                }
            }
        }
    }

    B.final_result = 0.0;
}

void run_computation() {
    int o = B.num_occupied_orbitals;
    int v = B.num_virtual_orbitals;
    double sum = 0.0;

    // This loop simulates the most expensive term in CCSD theory, scaling as O(o^2 * v^4).
    // new_t_ijab += Sum(c,d) <cd||ab> * t_ij^cd
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    double term_contribution = 0.0;
                    for (int c = 0; c < v; c++) {
                        for (int d = 0; d < v; d++) {
                            term_contribution += B.v_vvvv[c][d][a][b] * B.t2_oovv[i][j][c][d];
                        }
                    }
                    B.t2_new_oovv[i][j][a][b] = term_contribution; 
                }
            }
        }
    }
    
    // Accumulate the result to prevent dead code elimination
    for (int i = 0; i < o; i++) {
        for (int j = 0; j < o; j++) {
            for (int a = 0; a < v; a++) {
                for (int b = 0; b < v; b++) {
                    sum += B.t2_new_oovv[i][j][a][b];
                }
            }
        }
    }
    B.final_result = sum;
}

void cleanup() {
    int o = B.num_occupied_orbitals;
    int v = B.num_virtual_orbitals;
    free_4d_double(B.v_vvvv, v, v, v);
    free_4d_double(B.t2_oovv, o, o, v);
    free_4d_double(B.t2_new_oovv, o, o, v);
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
    printf("%f\n", B.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
