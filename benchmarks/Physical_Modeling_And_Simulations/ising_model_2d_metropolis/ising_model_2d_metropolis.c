#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// Mersenne Twister (MT19937) PRNG - DO NOT MODIFY
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
// End of Mersenne Twister

// --- Benchmark Globals ---
int lattice_size;
unsigned long long num_mc_steps;
double temperature;

// Lattice stored as a 1D array
int *lattice;
long long total_magnetization_sum;

// Pre-calculated Boltzmann factors for positive energy changes
double boltz_factor_4;
double boltz_factor_8;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <lattice_size> <num_mc_steps> <temperature> <seed>\n", argv[0]);
        exit(1);
    }

    lattice_size = atoi(argv[1]);
    num_mc_steps = strtoull(argv[2], NULL, 10);
    temperature = atof(argv[3]);
    uint32_t seed = atoi(argv[4]);

    if (lattice_size <= 0 || num_mc_steps == 0 || temperature <= 0) {
        fprintf(stderr, "Invalid arguments.\n");
        exit(1);
    }

    mt_seed(seed);

    lattice = (int *)malloc((size_t)lattice_size * lattice_size * sizeof(int));
    if (lattice == NULL) {
        fprintf(stderr, "Failed to allocate memory for the lattice.\n");
        exit(1);
    }

    for (int i = 0; i < lattice_size * lattice_size; ++i) {
        lattice[i] = (mt_rand() % 2) * 2 - 1; // Initialize with random spins: -1 or 1
    }

    // Pre-calculate Boltzmann factors for performance
    boltz_factor_4 = exp(-4.0 / temperature);
    boltz_factor_8 = exp(-8.0 / temperature);

    total_magnetization_sum = 0;
}

void run_computation() {
    const int L = lattice_size;
    const unsigned long long updates_per_step = (unsigned long long)L * L;
    const double max_rand_val = (double)UINT32_MAX;

    for (unsigned long long step = 0; step < num_mc_steps; ++step) {
        // Perform L*L Monte Carlo spin-flip attempts (one MC step)
        for (unsigned long long i = 0; i < updates_per_step; ++i) {
            // 1. Pick a random spin site
            int x = mt_rand() % L;
            int y = mt_rand() % L;
            int site_idx = x * L + y;

            // 2. Calculate sum of neighbor spins with periodic boundary conditions
            int xp1 = (x + 1) % L;
            int xm1 = (x - 1 + L) % L;
            int yp1 = (y + 1) % L;
            int ym1 = (y - 1 + L) % L;

            int neighbor_sum = lattice[xm1 * L + y] +
                               lattice[xp1 * L + y] +
                               lattice[x * L + ym1] +
                               lattice[x * L + yp1];

            // 3. Calculate energy change if spin is flipped
            int delta_E = 2 * lattice[site_idx] * neighbor_sum;

            // 4. Metropolis-Hastings acceptance rule
            if (delta_E <= 0) {
                lattice[site_idx] *= -1; // Flip the spin
            } else {
                double r = (double)mt_rand() / max_rand_val;
                if (delta_E == 4) {
                    if (r < boltz_factor_4) {
                        lattice[site_idx] *= -1;
                    }
                } else if (delta_E == 8) {
                    if (r < boltz_factor_8) {
                        lattice[site_idx] *= -1;
                    }
                }
            }
        }

        // Calculate and accumulate magnetization for this step
        long long current_magnetization = 0;
        for (int j = 0; j < L * L; ++j) {
            current_magnetization += lattice[j];
        }
        total_magnetization_sum += current_magnetization;
    }
}

void cleanup() {
    free(lattice);
}

// --- Main Function ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    // The final result to be printed.
    // We use the sum of magnetizations to prevent the compiler from optimizing away the computation.
    long long final_result = total_magnetization_sum;

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
