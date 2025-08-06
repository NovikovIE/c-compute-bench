#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) a_i32 Generator ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---

// A struct to hold 3D coordinates of an amino acid residue.
struct AminoAcid {
    double x, y, z;
};

int PROTEIN_LEN;
struct AminoAcid *protein_chain;
int *contact_map;
long total_contacts; // Final result accumulator


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s protein_sequence_length seed\n", argv[0]);
        exit(1);
    }

    PROTEIN_LEN = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);

    if (PROTEIN_LEN <= 0) {
        fprintf(stderr, "ERROR: protein_sequence_length must be a positive integer.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory for the protein chain (array of amino acid coordinates)
    protein_chain = (struct AminoAcid*)malloc(PROTEIN_LEN * sizeof(struct AminoAcid));
    if (!protein_chain) {
        fprintf(stderr, "FATAL: Memory allocation failed for protein_chain.\n");
        exit(1);
    }

    // Allocate memory for the contact map (a 2D matrix)
    // Use size_t for multiplication to prevent overflow with large PROTEIN_LEN
    size_t map_size = (size_t)PROTEIN_LEN * PROTEIN_LEN;
    contact_map = (int*)malloc(map_size * sizeof(int));
    if (!contact_map) {
        fprintf(stderr, "FATAL: Memory allocation failed for contact_map.\n");
        free(protein_chain);
        exit(1);
    }

    // Populate the protein chain with random 3D coordinates
    for (int i = 0; i < PROTEIN_LEN; i++) {
        // Scale coordinates to a typical simulation box size (e.g., 0-100 Angstroms)
        protein_chain[i].x = ((double)mt_rand() / (double)UINT32_MAX) * 100.0;
        protein_chain[i].y = ((double)mt_rand() / (double)UINT32_MAX) * 100.0;
        protein_chain[i].z = ((double)mt_rand() / (double)UINT32_MAX) * 100.0;
    }

    total_contacts = 0; // Initialize result
}

void run_computation() {
    long contacts = 0;
    const double CONTACT_THRESHOLD = 8.0; // Contact distance in Angstroms
    const double THRESHOLD_SQ = CONTACT_THRESHOLD * CONTACT_THRESHOLD;

    // Iterate over all unique pairs of amino acids
    for (int i = 0; i < PROTEIN_LEN; ++i) {
        // Start j from i + 1 to calculate for each pair only once
        for (int j = i + 1; j < PROTEIN_LEN; ++j) {
            double dx = protein_chain[i].x - protein_chain[j].x;
            double dy = protein_chain[i].y - protein_chain[j].y;
            double dz = protein_chain[i].z - protein_chain[j].z;
            double dist_sq = dx * dx + dy * dy + dz * dz;

            // Compare squared distance to avoid expensive sqrt()
            if (dist_sq < THRESHOLD_SQ) {
                // The pair is in contact
                contact_map[(size_t)i * PROTEIN_LEN + j] = 1;
                contact_map[(size_t)j * PROTEIN_LEN + i] = 1; // Map is symmetric
                contacts++;
            } else {
                contact_map[(size_t)i * PROTEIN_LEN + j] = 0;
                contact_map[(size_t)j * PROTEIN_LEN + i] = 0;
            }
        }
    }
    // The final accumulated result is the total number of contacting pairs
    total_contacts = contacts;
}

void cleanup() {
    free(protein_chain);
    free(contact_map);
}


// --- Main Driver ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (total number of contacts) to stdout
    printf("%ld\n", total_contacts);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
