#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- Mersenne Twister (Do Not Modify) ---
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
// --- END Mersenne Twister ---

// --- Benchmark Configuration ---
#define READ_LENGTH 100
#define SNP_RATE 0.01

// --- Global Data Structures ---
long param_genome_region_length;
int param_read_depth;
long total_reads;
long total_variants;

char* reference_genome;
char** reads;
int* read_starts;
int** pileup_counts; // [genome_position][base_index] where A=0, C=1, G=2, T=3

const char BASES[] = {'A', 'C', 'G', 'T'};

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
static inline char random_base() {
    return BASES[mt_rand() % 4];
}

// --- Benchmark Implementation ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s genome_region_length read_depth seed\n", argv[0]);
        exit(1);
    }

    param_genome_region_length = atol(argv[1]);
    param_read_depth = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    if (param_genome_region_length < READ_LENGTH) {
        fprintf(stderr, "FATAL: genome_region_length must be >= READ_LENGTH (%d)\n", READ_LENGTH);
        exit(1);
    }

    // 1. Generate reference genome
    reference_genome = (char*)malloc(param_genome_region_length * sizeof(char));
    if (!reference_genome) {
        fprintf(stderr, "FATAL: Memory allocation failed for reference genome.\n");
        exit(1);
    }
    for (long i = 0; i < param_genome_region_length; ++i) {
        reference_genome[i] = random_base();
    }

    // 2. Calculate number of reads and allocate containers
    total_reads = (param_genome_region_length * param_read_depth) / READ_LENGTH;
    reads = (char**)malloc(total_reads * sizeof(char*));
    read_starts = (int*)malloc(total_reads * sizeof(int));
    if (!reads || !read_starts) {
        fprintf(stderr, "FATAL: Memory allocation failed for reads array.\n");
        exit(1);
    }
    
    // 3. Generate reads based on reference
    long max_start_pos = param_genome_region_length - READ_LENGTH;
    for (long i = 0; i < total_reads; ++i) {
        reads[i] = (char*)malloc((READ_LENGTH + 1) * sizeof(char));
        if (!reads[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for a read string.\n");
            exit(1);
        }
        
        read_starts[i] = mt_rand() % (max_start_pos + 1);

        memcpy(reads[i], &reference_genome[read_starts[i]], READ_LENGTH);
        reads[i][READ_LENGTH] = '\0';

        // Introduce SNPs (mutations)
        for (int j = 0; j < READ_LENGTH; ++j) {
            double random_val = (double)mt_rand() / UINT32_MAX;
            if (random_val < SNP_RATE) {
                char original_base = reads[i][j];
                char new_base;
                do {
                    new_base = random_base();
                } while (new_base == original_base);
                reads[i][j] = new_base;
            }
        }
    }

    // 4. Allocate and initialize pileup counts table
    pileup_counts = (int**)malloc(param_genome_region_length * sizeof(int*));
    if (!pileup_counts) {
        fprintf(stderr, "FATAL: Memory allocation failed for pileup_counts pointers.\n");
        exit(1);
    }
    for (long i = 0; i < param_genome_region_length; ++i) {
        pileup_counts[i] = (int*)malloc(4 * sizeof(int));
        if (!pileup_counts[i]) {
            fprintf(stderr, "FATAL: Memory allocation failed for pileup_counts row.\n");
            exit(1);
        }
        memset(pileup_counts[i], 0, 4 * sizeof(int));
    }
    
    total_variants = 0;
}

void run_computation() {
    // Phase 1: Build the pileup table from all reads
    for (long i = 0; i < total_reads; ++i) {
        int start_pos = read_starts[i];
        for (int j = 0; j < READ_LENGTH; ++j) {
            long genome_pos = start_pos + j;
            char base = reads[i][j];
            switch (base) {
                case 'A': pileup_counts[genome_pos][0]++; break;
                case 'C': pileup_counts[genome_pos][1]++; break;
                case 'G': pileup_counts[genome_pos][2]++; break;
                case 'T': pileup_counts[genome_pos][3]++; break;
            }
        }
    }

    // Phase 2: Call variants from the pileup table
    long variant_count = 0;
    for (long i = 0; i < param_genome_region_length; ++i) {
        char ref_base = reference_genome[i];
        int max_count = -1;
        int consensus_base_idx = -1;

        for (int j = 0; j < 4; ++j) {
            if (pileup_counts[i][j] > max_count) {
                max_count = pileup_counts[i][j];
                consensus_base_idx = j;
            }
        }
        
        if (consensus_base_idx != -1) {
            char consensus_base = BASES[consensus_base_idx];
            if (consensus_base != ref_base) {
                variant_count++;
            }
        }
    }
    total_variants = variant_count;
}

void cleanup() {
    if (pileup_counts) {
        for (long i = 0; i < param_genome_region_length; ++i) {
            free(pileup_counts[i]);
        }
        free(pileup_counts);
    }
    
    if (reads) {
        for (long i = 0; i < total_reads; ++i) {
            free(reads[i]);
        }
        free(reads);
    }
    
    free(read_starts);
    free(reference_genome);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%ld\n", total_variants);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
