#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- MERSENNE TWISTER (from prompt) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA STRUCTURES ---

// A Term is a sequence of integers, terminated by 0.
// Positive integers represent function symbols.
typedef int* Term;

// A Rule is a pair of terms, lhs -> rhs.
// For simplicity, we enforce that length(lhs) == length(rhs).
typedef struct {
    Term lhs;
    Term rhs;
    int len; // Cache rule length for performance
} Rule;

// Global state
int num_rules;
int max_term_length;
int num_pairs_to_check;

Rule* rules = NULL;
Term* initial_pairs1 = NULL;
Term* initial_pairs2 = NULL;

int final_result = 0;

// --- HELPER FUNCTIONS ---

// Helper to get term length, stops at the 0 terminator.
int get_term_len(Term t) {
    int len = 0;
    while (t[len] != 0) {
        len++;
    }
    return len;
}

// --- COMPUTATION KERNEL ---

// Rewrite a term until it reaches a normal form (no more rules apply),
// then return a simple hash of the resulting term.
unsigned long normalize_and_hash(Term t) {
    Term work_term = (Term)malloc((max_term_length + 1) * sizeof(int));
    if (!work_term) {
        fprintf(stderr, "FATAL: Failed to allocate work_term.\n");
        exit(1);
    }
    memcpy(work_term, t, (max_term_length + 1) * sizeof(int));

    int work_term_len = get_term_len(work_term);
    int changed = 1;

    while (changed) {
        changed = 0;
        for (int i = 0; i < num_rules; ++i) {
            int rule_len = rules[i].len;
            for (int j = 0; j <= work_term_len - rule_len; ++j) {
                if (memcmp(work_term + j, rules[i].lhs, rule_len * sizeof(int)) == 0) {
                    memcpy(work_term + j, rules[i].rhs, rule_len * sizeof(int));
                    changed = 1;
                    goto next_pass; // Apply first possible rule and restart scan
                }
            }
        }
    next_pass:;
    }

    unsigned long hash = 5381;
    for (int i = 0; i < work_term_len; ++i) {
        hash = ((hash << 5) + hash) + work_term[i]; // djb2 hash
    }

    free(work_term);
    return hash;
}


// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_rules> <max_term_length> <seed>\n", argv[0]);
        exit(1);
    }

    num_rules = atoi(argv[1]);
    max_term_length = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    // Allocate rewrite rules
    rules = (Rule*)malloc(num_rules * sizeof(Rule));
    if (!rules) { fprintf(stderr, "FATAL: malloc failed for rules\n"); exit(1); }

    for (int i = 0; i < num_rules; ++i) {
        // Ensure rules are not too short or too long
        int len = (mt_rand() % (max_term_length / 2)) + 4;
        if (len >= max_term_length) len = max_term_length - 1;

        rules[i].len = len;
        rules[i].lhs = (Term)malloc((len + 1) * sizeof(int));
        rules[i].rhs = (Term)malloc((len + 1) * sizeof(int));
        if (!rules[i].lhs || !rules[i].rhs) { fprintf(stderr, "FATAL: malloc failed for rule terms\n"); exit(1); }

        for (int j = 0; j < len; ++j) {
            rules[i].lhs[j] = mt_rand() % 1000 + 1; // 1-1000 are function symbols
            rules[i].rhs[j] = mt_rand() % 1000 + 1;
        }
        rules[i].lhs[len] = 0; // Null terminator
        rules[i].rhs[len] = 0; // Null terminator
    }

    // Generate initial term pairs to check for confluence
    num_pairs_to_check = num_rules * 2;
    initial_pairs1 = (Term*)malloc(num_pairs_to_check * sizeof(Term));
    initial_pairs2 = (Term*)malloc(num_pairs_to_check * sizeof(Term));
    if (!initial_pairs1 || !initial_pairs2) { fprintf(stderr, "FATAL: malloc failed for initial pairs\n"); exit(1); }

    for (int i = 0; i < num_pairs_to_check; ++i) {
        initial_pairs1[i] = (Term)malloc((max_term_length + 1) * sizeof(int));
        initial_pairs2[i] = (Term)malloc((max_term_length + 1) * sizeof(int));
        if (!initial_pairs1[i] || !initial_pairs2[i]) { fprintf(stderr, "FATAL: malloc failed for pair terms\n"); exit(1); }

        int len1 = (mt_rand() % (max_term_length - 1)) + 1;
        int len2 = (mt_rand() % (max_term_length - 1)) + 1;

        for (int j = 0; j < len1; ++j) initial_pairs1[i][j] = mt_rand() % 1000 + 1;
        initial_pairs1[i][len1] = 0;

        for (int j = 0; j < len2; ++j) initial_pairs2[i][j] = mt_rand() % 1000 + 1;
        initial_pairs2[i][len2] = 0;
    }
}

void run_computation() {
    int joinable_pairs_count = 0;
    for (int i = 0; i < num_pairs_to_check; ++i) {
        unsigned long hash1 = normalize_and_hash(initial_pairs1[i]);
        unsigned long hash2 = normalize_and_hash(initial_pairs2[i]);
        if (hash1 == hash2) {
            joinable_pairs_count++;
        }
    }
    final_result = joinable_pairs_count;
}

void cleanup() {
    if (rules) {
        for (int i = 0; i < num_rules; ++i) {
            free(rules[i].lhs);
            free(rules[i].rhs);
        }
        free(rules);
    }
    if (initial_pairs1) {
        for (int i = 0; i < num_pairs_to_check; ++i) {
            free(initial_pairs1[i]);
        }
        free(initial_pairs1);
    }
    if (initial_pairs2) {
        for (int i = 0; i < num_pairs_to_check; ++i) {
            free(initial_pairs2[i]);
        }
        free(initial_pairs2);
    }
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
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
