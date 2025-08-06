#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) psychic-snail --- DO NOT MODIFY ---
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

// Benchmark Parameters
int SEARCH_DEPTH;
int BRANCHING_FACTOR;
int BEAM_WIDTH;

// Data structures
typedef struct {
    long score; // Represents the evaluation of a state in the search
} Node;

Node* current_beam; // Holds the best nodes from the previous level
Node* candidate_nodes; // Holds all generated nodes for the current level

// Final result accumulator
long final_result;

// Comparison function for qsort to sort nodes by score in descending order
int compare_nodes(const void* a, const void* b) {
    long score_a = ((const Node*)a)->score;
    long score_b = ((const Node*)b)->score;
    if (score_b > score_a) return 1;
    if (score_b < score_a) return -1;
    return 0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <search_depth> <branching_factor> <beam_width> <seed>\n", argv[0]);
        exit(1);
    }

    SEARCH_DEPTH = atoi(argv[1]);
    BRANCHING_FACTOR = atoi(argv[2]);
    BEAM_WIDTH = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    // Allocate memory for the beams
    current_beam = (Node*)malloc(BEAM_WIDTH * sizeof(Node));
    if (!current_beam) {
        fprintf(stderr, "Failed to allocate memory for current_beam\n");
        exit(1);
    }

    size_t candidate_size = (size_t)BEAM_WIDTH * BRANCHING_FACTOR;
    candidate_nodes = (Node*)malloc(candidate_size * sizeof(Node));
    if (!candidate_nodes) {
        fprintf(stderr, "Failed to allocate memory for candidate_nodes\n");
        free(current_beam);
        exit(1);
    }

    // Initialize the search with a single root node
    for (int i = 0; i < BEAM_WIDTH; ++i) {
        current_beam[i].score = 0; // Inactive slot
    }
    current_beam[0].score = 1; // Active root node with a starting score
}

void run_computation() {
    for (int depth = 0; depth < SEARCH_DEPTH; ++depth) {
        size_t candidate_count = 0;
        size_t max_candidates = (size_t)BEAM_WIDTH * BRANCHING_FACTOR;

        // 1. Generate candidate nodes from the current beam
        for (int i = 0; i < BEAM_WIDTH; ++i) {
            if (current_beam[i].score <= 0) continue; // Skip inactive nodes

            for (int j = 0; j < BRANCHING_FACTOR; ++j) {
                if (candidate_count >= max_candidates) break; 

                long score_perturbation = (mt_rand() % 201) - 100; //-100 to +100
                candidate_nodes[candidate_count].score = current_beam[i].score + score_perturbation;
                candidate_count++;
            }
        }

        // If no more active nodes, the search ends
        if (candidate_count == 0) break;

        // 2. Sort all generated candidates to find the best ones
        qsort(candidate_nodes, candidate_count, sizeof(Node), compare_nodes);

        // 3. Select the top BEAM_WIDTH nodes to form the next beam
        size_t num_to_copy = (candidate_count < (size_t)BEAM_WIDTH) ? candidate_count : (size_t)BEAM_WIDTH;
        for (size_t i = 0; i < num_to_copy; ++i) {
            current_beam[i] = candidate_nodes[i];
        }
        
        // Mark remaining slots in the beam as inactive
        for (size_t i = num_to_copy; i < (size_t)BEAM_WIDTH; ++i) {
            current_beam[i].score = 0;
        }
    }

    // Accumulate the final scores to prevent dead code elimination
    final_result = 0;
    for (int i = 0; i < BEAM_WIDTH; ++i) {
        final_result += current_beam[i].score;
    }
}

void cleanup() {
    free(current_beam);
    free(candidate_nodes);
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
    printf("%ld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
