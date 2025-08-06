#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// Mersenne Twister (MT19937) - Do Not Modify
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

// Benchmark-specific global variables
int NUM_INDIVIDUALS;
int NUM_SNPS;
int NUM_ITERATIONS;
unsigned long long final_result = 0;

// Genotypes: 0 = hom_ref, 1 = het, 2 = hom_alt
uint8_t **genotypes;
// Haplotypes: 0 = ref allele, 1 = alt allele
uint8_t ***haplotypes;

// A simple metric to score how well a haplotype segment matches another
int score_match(const uint8_t* hap1, const uint8_t* hap2, int start, int end) {
    int score = 0;
    for (int k = start; k < end; k++) {
        if (hap1[k] == hap2[k]) {
            score++;
        }
    }
    return score;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_individuals> <num_snps> <num_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_INDIVIDUALS = atoi(argv[1]);
    NUM_SNPS = atoi(argv[2]);
    NUM_ITERATIONS = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory
    genotypes = (uint8_t **)malloc(NUM_INDIVIDUALS * sizeof(uint8_t *));
    haplotypes = (uint8_t ***)malloc(NUM_INDIVIDUALS * sizeof(uint8_t **));
    for (int i = 0; i < NUM_INDIVIDUALS; i++) {
        genotypes[i] = (uint8_t *)malloc(NUM_SNPS * sizeof(uint8_t));
        haplotypes[i] = (uint8_t **)malloc(2 * sizeof(uint8_t *));
        haplotypes[i][0] = (uint8_t *)malloc(NUM_SNPS * sizeof(uint8_t));
        haplotypes[i][1] = (uint8_t *)malloc(NUM_SNPS * sizeof(uint8_t));
    }

    // Initialize data: Create random genotypes and an initial random phasing for heterozygous sites
    for (int i = 0; i < NUM_INDIVIDUALS; i++) {
        for (int j = 0; j < NUM_SNPS; j++) {
            genotypes[i][j] = mt_rand() % 3; // 0, 1, or 2
            
            if (genotypes[i][j] == 0) { // Homozygous reference
                haplotypes[i][0][j] = 0;
                haplotypes[i][1][j] = 0;
            } else if (genotypes[i][j] == 2) { // Homozygous alternate
                haplotypes[i][0][j] = 1;
                haplotypes[i][1][j] = 1;
            } else { // Heterozygous, assign a random initial phase
                if (mt_rand() % 2 == 0) {
                    haplotypes[i][0][j] = 0;
                    haplotypes[i][1][j] = 1;
                } else {
                    haplotypes[i][0][j] = 1;
                    haplotypes[i][1][j] = 0;
                }
            }
        }
    }
}

void run_computation() {
    unsigned long long checksum = 0;
    const int WINDOW_SIZE = 20; // Haplotype segment to compare
    const int NUM_REF_INDS = 10; // Number of other individuals to compare against

    for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
        for (int i = 0; i < NUM_INDIVIDUALS; i++) {
            for (int j = 0; j < NUM_SNPS; j++) {
                if (genotypes[i][j] != 1) continue; // Only phase heterozygous sites

                // Test phase 1: (0, 1)
                haplotypes[i][0][j] = 0;
                haplotypes[i][1][j] = 1;
                int score1 = 0;
                for (int k_idx = 0; k_idx < NUM_REF_INDS; k_idx++) {
                    int k = mt_rand() % NUM_INDIVIDUALS;
                    if (i == k) continue;

                    int start = (j > WINDOW_SIZE / 2) ? (j - WINDOW_SIZE / 2) : 0;
                    int end = (j < NUM_SNPS - WINDOW_SIZE / 2) ? (j + WINDOW_SIZE / 2) : NUM_SNPS;

                    int match_a0 = score_match(haplotypes[i][0], haplotypes[k][0], start, end);
                    int match_a1 = score_match(haplotypes[i][0], haplotypes[k][1], start, end);
                    int match_b0 = score_match(haplotypes[i][1], haplotypes[k][0], start, end);
                    int match_b1 = score_match(haplotypes[i][1], haplotypes[k][1], start, end);
                    score1 += (match_a0 > match_a1 ? match_a0 : match_a1) + (match_b0 > match_b1 ? match_b0 : match_b1);
                }

                // Test phase 2: (1, 0)
                haplotypes[i][0][j] = 1;
                haplotypes[i][1][j] = 0;
                int score2 = 0;
                for (int k_idx = 0; k_idx < NUM_REF_INDS; k_idx++) {
                    int k = mt_rand() % NUM_INDIVIDUALS;
                    if (i == k) continue;

                    int start = (j > WINDOW_SIZE / 2) ? (j - WINDOW_SIZE / 2) : 0;
                    int end = (j < NUM_SNPS - WINDOW_SIZE / 2) ? (j + WINDOW_SIZE / 2) : NUM_SNPS;

                    int match_a0 = score_match(haplotypes[i][0], haplotypes[k][0], start, end);
                    int match_a1 = score_match(haplotypes[i][0], haplotypes[k][1], start, end);
                    int match_b0 = score_match(haplotypes[i][1], haplotypes[k][0], start, end);
                    int match_b1 = score_match(haplotypes[i][1], haplotypes[k][1], start, end);
                    score2 += (match_a0 > match_a1 ? match_a0 : match_a1) + (match_b0 > match_b1 ? match_b0 : match_b1);
                }

                // Choose the better phasing (deterministically for simplicity)
                if (score2 > score1) {
                    checksum += score2; // Phase (1, 0) is already set
                } else {
                    haplotypes[i][0][j] = 0; // Revert to phase (0, 1)
                    haplotypes[i][1][j] = 1;
                    checksum += score1;
                }
            }
        }
    }
    final_result = checksum;
}

void cleanup() {
    for (int i = 0; i < NUM_INDIVIDUALS; i++) {
        free(genotypes[i]);
        free(haplotypes[i][0]);
        free(haplotypes[i][1]);
        free(haplotypes[i]);
    }
    free(genotypes);
    free(haplotypes);
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
    printf("%llu\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
