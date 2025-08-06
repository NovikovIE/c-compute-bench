#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---

// --- Benchmark Globals ---
typedef struct {
    uint64_t key;
    int count;
} HashTableEntry;

// Parameters
long long GENOME_SIZE;
int KMER_SIZE;

// Data
char* genome;
HashTableEntry* hash_table;

// Hash table specifics
const long long HASH_TABLE_SIZE = 1LL << 29; // 536,870,912 entries
uint64_t KMER_MASK;

// Result
long long final_result;

// --- Helper Functions ---
// Convert base 'A', 'C', 'G', 'T' to 2-bit representation
static inline uint64_t encode_base(char base) {
    switch(base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    }
    return 0; // Should not happen with generated data
}

// Simple integer mixer hash function
static inline uint64_t hash_func(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <genome_size_mb> <kmer_size> <seed>\n", argv[0]);
        exit(1);
    }

    long long genome_size_mb = atoll(argv[1]);
    KMER_SIZE = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if(genome_size_mb <= 0 || KMER_SIZE <= 0 || KMER_SIZE > 31) {
        fprintf(stderr, "Invalid arguments. genome_size_mb > 0, 0 < kmer_size <= 31\n");
        exit(1);
    }

    GENOME_SIZE = genome_size_mb * 1024 * 1024;
    KMER_MASK = (1ULL << (2 * KMER_SIZE)) - 1;

    mt_seed(seed);

    // Allocate genome sequence
    genome = (char*)malloc(GENOME_SIZE * sizeof(char));
    if (!genome) {
        fprintf(stderr, "Failed to allocate memory for genome.\n");
        exit(1);
    }

    // Generate random genome
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (long long i = 0; i < GENOME_SIZE; ++i) {
        genome[i] = bases[mt_rand() % 4];
    }
    
    // Allocate and initialize hash table
    hash_table = (HashTableEntry*)calloc(HASH_TABLE_SIZE, sizeof(HashTableEntry));
    if (!hash_table) {
        fprintf(stderr, "Failed to allocate memory for hash table.\n");
        free(genome);
        exit(1);
    }
}

void run_computation() {
    if (GENOME_SIZE < KMER_SIZE) {
        final_result = 0;
        return;
    }

    // --- Initial k-mer (of length k-1) ---
    uint64_t current_kmer = 0;
    for (int i = 0; i < KMER_SIZE - 1; ++i) {
        current_kmer = (current_kmer << 2) | encode_base(genome[i]);
    }
    
    // --- K-mer counting with rolling window ---
    long long num_kmers = GENOME_SIZE - KMER_SIZE + 1;
    long long hash_mask = HASH_TABLE_SIZE - 1;

    for (long long i = 0; i < num_kmers; ++i) {
        // Roll forward to build the full k-mer
        current_kmer = ((current_kmer << 2) & KMER_MASK) | encode_base(genome[i + KMER_SIZE - 1]);

        // Find slot in hash table using open addressing with linear probing
        uint64_t idx = hash_func(current_kmer) & hash_mask;

        while (1) {
            if (hash_table[idx].count == 0) {
                // Found empty slot, insert new k-mer
                hash_table[idx].key = current_kmer;
                hash_table[idx].count = 1;
                break;
            } else if (hash_table[idx].key == current_kmer) {
                // Found existing k-mer, increment count
                hash_table[idx].count++;
                break;
            }
            // Collision, move to next slot
            idx = (idx + 1) & hash_mask;
        }
    }
    
    // Calculate a final result to prevent dead code elimination
    long long checksum = 0;
    for (long long i = 0; i < HASH_TABLE_SIZE; ++i) {
        if (hash_table[i].count > 0) {
            checksum += (hash_table[i].key % 997) + hash_table[i].count;
        }
    }
    final_result = checksum;
}

void cleanup() {
    free(genome);
    free(hash_table);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    
    final_result = 0; // Initialize result

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
