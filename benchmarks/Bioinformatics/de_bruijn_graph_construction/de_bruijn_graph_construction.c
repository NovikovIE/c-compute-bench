#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// Mersenne Twister (MT19937) PRNG
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

// Parameters
long num_reads;
int read_length;
int kmer_size;

// Input Data
char **reads;

// Hash Table for De Bruijn graph nodes (k-mers)
// Using a fixed-size hash table with chaining for collision resolution.
#define HASH_TABLE_SIZE (1 << 22) // Approx. 4.2 million slots
typedef struct KmerNode {
    const char *kmer;
    struct KmerNode *next;
} KmerNode;

KmerNode **hash_table;
long unique_kmers_count;

// Memory arenas for performant allocation of nodes and k-mer strings.
// This avoids millions of small malloc calls within the timed section.
char *kmer_arena;
size_t kmer_arena_size;
size_t kmer_arena_offset;

KmerNode *node_arena;
size_t node_arena_size;
size_t node_arena_offset;

// --- Utility Function ---

// DJB2 hash function for a string of a given length
static inline unsigned long hash_kmer(const char *str, int len) {
    unsigned long hash = 5381;
    const char *end = str + len;
    while (str < end) {
        hash = ((hash << 5) + hash) + *str++; // hash * 33 + c
    }
    return hash;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_reads read_length kmer_size seed\n", argv[0]);
        exit(1);
    }

    num_reads = atol(argv[1]);
    read_length = atoi(argv[2]);
    kmer_size = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (num_reads <= 0 || read_length <= 0 || kmer_size <= 0) {
        fprintf(stderr, "FATAL: Parameters must be positive.\n");
        exit(1);
    }
    if (kmer_size >= read_length) {
        fprintf(stderr, "FATAL: kmer_size must be smaller than read_length.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate and generate random DNA reads
    const char bases[] = "ACGT";
    reads = (char **)malloc(num_reads * sizeof(char *));
    if (!reads) { fprintf(stderr, "FATAL: Memory allocation for reads failed.\n"); exit(1); }
    for (long i = 0; i < num_reads; ++i) {
        reads[i] = (char *)malloc((read_length + 1) * sizeof(char));
        if (!reads[i]) { fprintf(stderr, "FATAL: Memory allocation for a read failed.\n"); exit(1); }
        for (int j = 0; j < read_length; ++j) {
            reads[i][j] = bases[mt_rand() % 4];
        }
        reads[i][read_length] = '\0';
    }

    // Allocate the hash table (array of pointers), initialized to NULL
    hash_table = (KmerNode **)calloc(HASH_TABLE_SIZE, sizeof(KmerNode *));
    if (!hash_table) { fprintf(stderr, "FATAL: Memory allocation for hash_table failed.\n"); exit(1); }

    // Allocate memory arenas. Sized to handle all possible k-mers if they were unique.
    long total_kmers = num_reads * (long)(read_length - kmer_size + 1);
    kmer_arena_size = total_kmers * (size_t)(kmer_size + 1);
    node_arena_size = total_kmers * sizeof(KmerNode);
    
    kmer_arena = (char*)malloc(kmer_arena_size);
    node_arena = (KmerNode*)malloc(node_arena_size);

    if (!kmer_arena || !node_arena) { fprintf(stderr, "FATAL: Memory allocation for arenas failed.\n"); exit(1); }
    
    kmer_arena_offset = 0;
    node_arena_offset = 0;
}

void run_computation() {
    unique_kmers_count = 0;

    for (long i = 0; i < num_reads; ++i) {
        for (int j = 0; j <= read_length - kmer_size; ++j) {
            const char *kmer_ptr = reads[i] + j;

            // Compute hash for the current k-mer
            unsigned long h = hash_kmer(kmer_ptr, kmer_size) % HASH_TABLE_SIZE;
            
            // Check if k-mer already exists in the hash table
            KmerNode *node = hash_table[h];
            int found = 0;
            while (node != NULL) {
                if (strncmp(node->kmer, kmer_ptr, kmer_size) == 0) {
                    found = 1;
                    break;
                }
                node = node->next;
            }

            // If not found, add it to the table
            if (!found) {
                // Check if arenas have enough space
                if (node_arena_offset >= node_arena_size / sizeof(KmerNode) ||
                    kmer_arena_offset + kmer_size + 1 >= kmer_arena_size) {
                    // For a benchmark, we pre-allocate enough and can stop here.
                    continue; 
                }

                // Get space for the k-mer string from the k-mer arena
                char *new_kmer_str = kmer_arena + kmer_arena_offset;
                strncpy(new_kmer_str, kmer_ptr, kmer_size);
                new_kmer_str[kmer_size] = '\0';
                kmer_arena_offset += (kmer_size + 1);
                
                // Get a new node from the node arena
                KmerNode *new_node = node_arena + node_arena_offset++;
                new_node->kmer = new_kmer_str;
                
                // Insert the new node at the head of the chain
                new_node->next = hash_table[h];
                hash_table[h] = new_node;

                unique_kmers_count++;
            }
        }
    }
}

void cleanup() {
    for (long i = 0; i < num_reads; ++i) {
        free(reads[i]);
    }
    free(reads);

    free(hash_table);
    free(kmer_arena);
    free(node_arena);
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
    printf("%ld\n", unique_kmers_count);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
