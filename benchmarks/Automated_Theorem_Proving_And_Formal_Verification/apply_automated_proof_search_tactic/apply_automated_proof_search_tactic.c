#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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
        int i;
        for (i = 0; i < MT_N - MT_M; i++) {
            y = (mt[i] & MT_UPPER_MASK) | (mt[i + 1] & MT_LOWER_MASK);
            mt[i] = mt[i + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (; i < MT_N - 1; i++) {
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

// --- BENCHMARK IMPLEMENTATION ---

// Represents a logical formula, lemma, or goal state.
// The 'terms' are abstract components, and 'complexity' is their count.
typedef struct {
    uint32_t* terms;
    int complexity;
} Formula;

// Global parameters and data structures
static int g_num_initial_lemmas;
static int g_goal_complexity;
static int g_search_depth;

static Formula* g_lemmas;       // Knowledge base of existing proofs/lemmas
static Formula g_goal;           // The main theorem to be proven
static long long g_proofs_found; // Accumulator for final result

// Recursive helper function to simulate the proof search.
// It navigates the state space by applying lemmas as tactics.
void apply_tactic_recursive(const Formula* current_goal, int depth) {
    // Base case 1: Terminate if the goal is proven (simplified enough)
    if (current_goal->complexity <= 2) {
        g_proofs_found++;
        return;
    }

    // Base case 2: Terminate if maximum search depth is reached
    if (depth >= g_search_depth) {
        return;
    }

    // Use a Variable Length Array (VLA) on the stack for the next goal's terms.
    // This avoids heap allocation overhead (malloc/free) in the hot computation loop.
    // `g_goal_complexity` sets its maximum size.
    uint32_t next_goal_terms[g_goal_complexity];

    // Recursive step: Try to apply each lemma as a transformation tactic
    for (int i = 0; i < g_num_initial_lemmas; ++i) {
        Formula next_goal;
        next_goal.terms = next_goal_terms;

        // Simulate applying a tactic: Create a new formula by combining the current
        // goal and a lemma. The new complexity is a deterministic combination of the two.
        next_goal.complexity = (current_goal->complexity + g_lemmas[i].complexity) % g_goal_complexity;
        if (next_goal.complexity == 0) {
            next_goal.complexity = 1;
        }

        // The new terms are derived by XORing the current terms and the lemma's terms.
        // This simulates a logical transformation.
        for (int j = 0; j < next_goal.complexity; ++j) {
            next_goal.terms[j] = current_goal->terms[j % current_goal->complexity] ^ g_lemmas[i].terms[j % g_lemmas[i].complexity];
        }

        // Recurse with the new goal and increased depth
        apply_tactic_recursive(&next_goal, depth + 1);
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_initial_lemmas> <goal_complexity> <search_depth> <seed>\n", argv[0]);
        exit(1);
    }

    g_num_initial_lemmas = atoi(argv[1]);
    g_goal_complexity = atoi(argv[2]);
    g_search_depth = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (g_num_initial_lemmas <= 0 || g_goal_complexity <= 0 || g_search_depth <= 0) {
        fprintf(stderr, "FATAL: Parameters must be positive integers.\n");
        exit(1);
    }
    // Safety check to prevent stack overflow from VLA in the recursive function
    if (g_goal_complexity > 8192) {
        fprintf(stderr, "FATAL: goal_complexity over 8192 is too large and may cause stack overflow.\n");
        exit(1);
    }

    mt_seed(seed);

    g_lemmas = (Formula*)malloc(g_num_initial_lemmas * sizeof(Formula));
    if (!g_lemmas) { perror("malloc failed"); exit(1); }

    for (int i = 0; i < g_num_initial_lemmas; ++i) {
        g_lemmas[i].complexity = (mt_rand() % g_goal_complexity) + 1;
        g_lemmas[i].terms = (uint32_t*)malloc(g_lemmas[i].complexity * sizeof(uint32_t));
        if (!g_lemmas[i].terms) { perror("malloc failed"); exit(1); }
        for (int j = 0; j < g_lemmas[i].complexity; ++j) {
            g_lemmas[i].terms[j] = mt_rand();
        }
    }

    g_goal.complexity = g_goal_complexity;
    g_goal.terms = (uint32_t*)malloc(g_goal.complexity * sizeof(uint32_t));
    if (!g_goal.terms) { perror("malloc failed"); exit(1); }
    for (int i = 0; i < g_goal.complexity; ++i) {
        g_goal.terms[i] = mt_rand();
    }
}

void run_computation() {
    g_proofs_found = 0;
    apply_tactic_recursive(&g_goal, 0);
}

void cleanup() {
    for (int i = 0; i < g_num_initial_lemmas; ++i) {
        free(g_lemmas[i].terms);
    }
    free(g_lemmas);
    free(g_goal.terms);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    long long final_result = g_proofs_found;

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
