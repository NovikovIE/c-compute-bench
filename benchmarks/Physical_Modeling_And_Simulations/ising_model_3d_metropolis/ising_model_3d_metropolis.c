#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

/*
 * Mersenne Twister (MT19937)
 * Do Not Modify - Included Verbatim
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
/*
 * End of Mersenne Twister
 */


// Benchmark parameters
static int L;
static unsigned int NUM_MC_STEPS;
static double T;
static unsigned int SEED;

// Data structures
static int* lattice;
// Pre-calculated Boltzmann factors for positive energy changes
// delta_E can be 4, 8, 12
static double exp_vals[3];

// Final result
static long long final_result;

// Forward declarations of benchmark functions
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();


void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <lattice_size> <num_mc_steps> <temperature> <seed>\n", argv[0]);
        exit(1);
    }
    
    L = atoi(argv[1]);
    NUM_MC_STEPS = (unsigned int)atoi(argv[2]);
    T = atof(argv[3]);
    SEED = (unsigned int)atoi(argv[4]);
    
    if (L <= 1 || NUM_MC_STEPS <= 0 || T <= 0.0) {
        fprintf(stderr, "FATAL: Invalid parameters. L must be > 1.\n");
        exit(1);
    }
    
    mt_seed(SEED);
    
    long long N = (long long)L * L * L;
    lattice = (int*)malloc(N * sizeof(int));
    if (!lattice) {
        fprintf(stderr, "FATAL: Failed to allocate memory for lattice.\n");
        exit(1);
    }
    
    // Initialize lattice with random spins (-1 or +1)
    for (long long i = 0; i < N; ++i) {
        lattice[i] = (mt_rand() % 2) * 2 - 1;
    }
    
    // Pre-calculate Boltzmann factors for possible positive energy changes
    // delta_E = 4, 8, 12
    exp_vals[0] = exp(-4.0 / T);
    exp_vals[1] = exp(-8.0 / T);
    exp_vals[2] = exp(-12.0 / T);
}

void run_computation() {
    long long N = (long long)L * L * L;
    unsigned long long total_attempts = (unsigned long long)NUM_MC_STEPS * N;
    long long L_sq = (long long)L * L;

    for (unsigned long long i = 0; i < total_attempts; ++i) {
        // 1. Pick a random spin site by its 1D index
        long long idx = mt_rand() % N;
        
        // 2. Calculate sum of 6 neighbors with periodic boundary conditions
        int z = idx / L_sq;
        int y = (idx % L_sq) / L;
        int x = idx % L;
        
        int xp1 = (x == L - 1) ? 0 : x + 1;
        int xm1 = (x == 0) ? L - 1 : x - 1;
        int yp1 = (y == L - 1) ? 0 : y + 1;
        int ym1 = (y == 0) ? L - 1 : y - 1;
        int zp1 = (z == L - 1) ? 0 : z + 1;
        int zm1 = (z == 0) ? L - 1 : z - 1;
        
        long long idx_xp1 = (long long)z * L_sq + (long long)y * L + xp1;
        long long idx_xm1 = (long long)z * L_sq + (long long)y * L + xm1;
        long long idx_yp1 = (long long)z * L_sq + (long long)yp1 * L + x;
        long long idx_ym1 = (long long)z * L_sq + (long long)ym1 * L + x;
        long long idx_zp1 = (long long)zp1 * L_sq + (long long)y * L + x;
        long long idx_zm1 = (long long)zm1 * L_sq + (long long)y * L + x;

        int neighbor_sum = lattice[idx_xp1] + lattice[idx_xm1] +
                           lattice[idx_yp1] + lattice[idx_ym1] +
                           lattice[idx_zp1] + lattice[idx_zm1];
            
        // 3. Calculate energy change if spin is flipped
        int delta_E = 2 * lattice[idx] * neighbor_sum;
        
        // 4. Metropolis-Hastings acceptance rule
        if (delta_E <= 0) {
            lattice[idx] *= -1; // Accept flip
        } else {
            // Accept with probability exp(-delta_E / T)
            // Use pre-calculated values. delta_E can only be 4, 8, or 12 here.
            // Map delta_E to index: 4->0, 8->1, 12->2
            if ((mt_rand() / 4294967295.0) < exp_vals[delta_E/4 - 1]) {
                lattice[idx] *= -1; // Accept flip
            }
        }
    }
    
    // Calculate final total magnetization to prevent dead code elimination
    long long total_magnetization = 0;
    for (long long j = 0; j < N; ++j) {
        total_magnetization += lattice[j];
    }
    final_result = total_magnetization;
}

void cleanup() {
    free(lattice);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
