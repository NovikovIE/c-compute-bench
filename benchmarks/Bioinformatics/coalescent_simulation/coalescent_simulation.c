#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

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
// --- End of Mersenne Twister ---

// Global parameters
int g_sample_size;
int g_num_loci;
int g_num_generations;

// Global data structures
double* g_locus_mrca_time; // Time to most recent common ancestor for each locus
double g_final_result;     // Final accumulated result to prevent optimization

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <sample_size> <num_loci> <num_generations> <seed>\n", argv[0]);
        exit(1);
    }

    g_sample_size = atoi(argv[1]);
    g_num_loci = atoi(argv[2]);
    g_num_generations = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);

    if (g_sample_size <= 1 || g_num_loci <= 0 || g_num_generations <= 0) {
        fprintf(stderr, "FATAL: Invalid parameters. sample_size must be > 1, others must be > 0.\n");
        exit(1);
    }

    mt_seed(seed);

    g_locus_mrca_time = (double*)malloc(g_num_loci * sizeof(double));
    if (g_locus_mrca_time == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for locus times.\n");
        exit(1);
    }

    g_final_result = 0.0;
}

void run_computation() {
    double total_mrca_time = 0.0;

    for (int i = 0; i < g_num_loci; i++) {
        int current_lineages = g_sample_size;
        double current_time = 0.0;

        while (current_lineages > 1 && current_time < g_num_generations) {
            // Generate a uniform random double in the open interval (0, 1) to avoid log(0)
            double u = (mt_rand() + 1.0) / (UINT32_MAX + 2.0);

            // The rate of coalescence for k lineages is k*(k-1)/2 (in coalescent units where time is scaled by 2Ne)
            double lambda = (double)(current_lineages) * (current_lineages - 1) / 2.0;
            
            // Time to the next coalescence event is exponentially distributed with the calculated rate.
            double time_to_event = -log(u) / lambda;

            current_time += time_to_event;
            current_lineages--;
        }

        g_locus_mrca_time[i] = current_time;
        total_mrca_time += current_time; // Accumulate result
    }
    
    g_final_result = total_mrca_time;
}

void cleanup() {
    free(g_locus_mrca_time);
    g_locus_mrca_time = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    
    setup_benchmark(argc, argv);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print final result to stdout to prevent dead code elimination
    printf("%f\n", g_final_result);
    
    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);
    
    return 0;
}
