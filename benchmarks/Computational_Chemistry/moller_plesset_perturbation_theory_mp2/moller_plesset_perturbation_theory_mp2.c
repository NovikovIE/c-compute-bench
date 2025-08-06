#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- Utility Function for Random Doubles ---
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Data Structures ---
typedef struct {
    int num_occupied_orbitals;
    int num_virtual_orbitals;
    
    double* orbital_energies; // Combined occupied and virtual
    double**** two_electron_integrals; // 4D Tensor (ia|jb)

    double mp2_correlation_energy;
} BenchmarkData;

static BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_occupied_orbitals> <num_virtual_orbitals> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_occupied_orbitals = atoi(argv[1]);
    g_data.num_virtual_orbitals = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);
    mt_seed(seed);

    int o = g_data.num_occupied_orbitals;
    int v = g_data.num_virtual_orbitals;
    int total_orbitals = o + v;

    // Allocate and initialize orbital energies
    g_data.orbital_energies = (double*)malloc(total_orbitals * sizeof(double));
    if (!g_data.orbital_energies) { perror("malloc failed"); exit(1); }
    for (int i = 0; i < o; ++i) {
        g_data.orbital_energies[i] = -1.0 - rand_double() * 2.0; // Occupied energies are negative
    }
    for (int i = o; i < total_orbitals; ++i) {
        g_data.orbital_energies[i] = 1.0 + rand_double() * 2.0; // Virtual energies are positive
    }

    // Allocate 4D tensor for two-electron repulsion integrals (ERIs)
    g_data.two_electron_integrals = (double****)malloc(o * sizeof(double***));
    if (!g_data.two_electron_integrals) { perror("malloc failed"); exit(1); }
    
    for (int i = 0; i < o; ++i) {
        g_data.two_electron_integrals[i] = (double***)malloc(v * sizeof(double**));
        if (!g_data.two_electron_integrals[i]) { perror("malloc failed"); exit(1); }
        for (int a = 0; a < v; ++a) {
            g_data.two_electron_integrals[i][a] = (double**)malloc(o * sizeof(double*));
            if (!g_data.two_electron_integrals[i][a]) { perror("malloc failed"); exit(1); }
            for (int j = 0; j < o; ++j) {
                g_data.two_electron_integrals[i][a][j] = (double*)malloc(v * sizeof(double));
                if (!g_data.two_electron_integrals[i][a][j]) { perror("malloc failed"); exit(1); }
                for (int b = 0; b < v; ++b) {
                    // Initialize with small random values
                    g_data.two_electron_integrals[i][a][j][b] = (rand_double() - 0.5) * 0.1;
                }
            }
        }
    }
}

void run_computation() {
    int o = g_data.num_occupied_orbitals;
    int v = g_data.num_virtual_orbitals;
    double energy_sum = 0.0;

    // MP2 correlation energy formula:
    // E_MP2 = sum_{i,j,a,b} [ (ia|jb) * (2*(ia|jb) - (ib|ja)) ] / (e_i + e_j - e_a - e_b)
    // i,j run over occupied orbitals; a,b run over virtual orbitals.

    for (int i = 0; i < o; ++i) {
        for (int j = 0; j < o; ++j) {
            for (int a = 0; a < v; ++a) {
                for (int b = 0; b < v; ++b) {
                    double integral_iajb = g_data.two_electron_integrals[i][a][j][b];
                    double integral_ibja = g_data.two_electron_integrals[i][b][j][a];
                    
                    double energy_i = g_data.orbital_energies[i];
                    double energy_j = g_data.orbital_energies[j];
                    double energy_a = g_data.orbital_energies[o + a];
                    double energy_b = g_data.orbital_energies[o + b];

                    double denominator = energy_i + energy_j - energy_a - energy_b;
                    double numerator = integral_iajb * (2.0 * integral_iajb - integral_ibja);

                    energy_sum += numerator / denominator;
                }
            }
        }
    }

    g_data.mp2_correlation_energy = energy_sum;
}

void cleanup() {
    int o = g_data.num_occupied_orbitals;
    int v = g_data.num_virtual_orbitals;

    for (int i = 0; i < o; ++i) {
        for (int a = 0; a < v; ++a) {
            for (int j = 0; j < o; ++j) {
                free(g_data.two_electron_integrals[i][a][j]);
            }
            free(g_data.two_electron_integrals[i][a]);
        }
        free(g_data.two_electron_integrals[i]);
    }
    free(g_data.two_electron_integrals);
    free(g_data.orbital_energies);
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
    printf("%f\n", g_data.mp2_correlation_energy);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
