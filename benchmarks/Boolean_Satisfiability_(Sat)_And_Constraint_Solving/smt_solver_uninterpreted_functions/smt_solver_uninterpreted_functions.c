#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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

// A term is either a variable or a unary function application f(arg).
// If arg is -1, it's a variable.
typedef struct {
    int arg;
} Term;

typedef struct {
    int t1;
    int t2;
} Constraint;

// Global structure to hold all benchmark data
typedef struct {
    int num_variables;
    int num_function_applications;
    int num_equalities;
    int num_disequalities;
    int total_terms;

    Term* terms;
    Constraint* equalities;
    Constraint* disequalities;

    // For computation
    int* parent;   // Union-find parent array
    int* rank;     // Union-find rank array
    int* func_app_indices; // Indices of function application terms
} BenchmarkData;

static BenchmarkData g_data;
static int g_final_result = 0;

// --- Union-Find Helper Functions ---
int find_set(int i) {
    if (g_data.parent[i] == i)
        return i;
    return g_data.parent[i] = find_set(g_data.parent[i]); // Path compression
}

int unite_sets(int i, int j) {
    int root_i = find_set(i);
    int root_j = find_set(j);
    if (root_i != root_j) {
        // Union by rank
        if (g_data.rank[root_i] < g_data.rank[root_j]) {
            int temp = root_i;
            root_i = root_j;
            root_j = temp;
        }
        g_data.parent[root_j] = root_i;
        if (g_data.rank[root_i] == g_data.rank[root_j]) {
            g_data.rank[root_i]++;
        }
        return 1; // A merge happened
    }
    return 0; // No merge happened
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_variables num_function_applications num_equalities num_disequalities seed\n", argv[0]);
        exit(1);
    }

    g_data.num_variables = atoi(argv[1]);
    g_data.num_function_applications = atoi(argv[2]);
    g_data.num_equalities = atoi(argv[3]);
    g_data.num_disequalities = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);

    mt_seed(seed);

    g_data.total_terms = g_data.num_variables + g_data.num_function_applications;

    // Allocate memory
    g_data.terms = (Term*)malloc(g_data.total_terms * sizeof(Term));
    g_data.equalities = (Constraint*)malloc(g_data.num_equalities * sizeof(Constraint));
    g_data.disequalities = (Constraint*)malloc(g_data.num_disequalities * sizeof(Constraint));
    g_data.parent = (int*)malloc(g_data.total_terms * sizeof(int));
    g_data.rank = (int*)malloc(g_data.total_terms * sizeof(int));
    g_data.func_app_indices = (int*)malloc(g_data.num_function_applications * sizeof(int));

    if (!g_data.terms || !g_data.equalities || !g_data.disequalities || !g_data.parent || !g_data.rank || !g_data.func_app_indices) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Generate terms
    int func_app_idx = 0;
    for (int i = 0; i < g_data.total_terms; ++i) {
        if (i < g_data.num_variables) {
            g_data.terms[i].arg = -1; // It's a variable
        } else {
            // It's a function application f(x), where x is a variable
            g_data.terms[i].arg = mt_rand() % g_data.num_variables;
            g_data.func_app_indices[func_app_idx++] = i;
        }
    }

    // Generate equalities
    for (int i = 0; i < g_data.num_equalities; ++i) {
        g_data.equalities[i].t1 = mt_rand() % g_data.total_terms;
        g_data.equalities[i].t2 = mt_rand() % g_data.total_terms;
    }

    // Generate disequalities
    for (int i = 0; i < g_data.num_disequalities; ++i) {
        g_data.disequalities[i].t1 = mt_rand() % g_data.total_terms;
        g_data.disequalities[i].t2 = mt_rand() % g_data.total_terms;
    }
}

void run_computation() {
    // Initialize Union-Find structure
    for (int i = 0; i < g_data.total_terms; ++i) {
        g_data.parent[i] = i;
        g_data.rank[i] = 0;
    }

    // 1. Process initial equalities
    for (int i = 0; i < g_data.num_equalities; ++i) {
        unite_sets(g_data.equalities[i].t1, g_data.equalities[i].t2);
    }

    // 2. Iteratively apply congruence rule: a=b => f(a)=f(b)
    int changed;
    do {
        changed = 0;
        for (int i = 0; i < g_data.num_function_applications; ++i) {
            for (int j = i + 1; j < g_data.num_function_applications; ++j) {
                int term_idx1 = g_data.func_app_indices[i];
                int term_idx2 = g_data.func_app_indices[j];

                int arg1 = g_data.terms[term_idx1].arg;
                int arg2 = g_data.terms[term_idx2].arg;

                if (find_set(arg1) == find_set(arg2)) {
                    if (unite_sets(term_idx1, term_idx2)) {
                        changed = 1;
                    }
                }
            }
        }
    } while (changed);

    // 3. Check disequalities for contradictions
    int contradictions = 0;
    for (int i = 0; i < g_data.num_disequalities; ++i) {
        if (find_set(g_data.disequalities[i].t1) == find_set(g_data.disequalities[i].t2)) {
            contradictions++;
        }
    }

    g_final_result = contradictions;
}

void cleanup() {
    free(g_data.terms);
    free(g_data.equalities);
    free(g_data.disequalities);
    free(g_data.parent);
    free(g_data.rank);
    free(g_data.func_app_indices);
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", g_final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
