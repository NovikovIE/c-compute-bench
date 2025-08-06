#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Benchmark Parameters
static int NUM_VARIABLES;
static int NUM_CLAUSES;
static int VAR_ORDERING_HEURISTIC;

// Data Structures
// A 3-SAT formula: array of clauses, each with 3 literals.
// A literal `v` is `v+1`, `~v` is `-(v+1)`. `v` is 0-indexed.
static int** clauses;
static const int CLAUSE_SIZE = 3;

// Variable ordering: maps a BDD level to a variable ID and vice versa.
static int* var_order; // var_order[level] = var_id
static int* var_level; // var_level[var_id] = level

// BDD Node: level is the variable index in the ordering.
typedef struct {
    int level;
    int low;   // Node ID for the 'false' branch
    int high;  // Node ID for the 'true' branch
} BDDNode;

#define BDD_FALSE 0
#define BDD_TRUE  1

// Unique table for BDD nodes to enforce reduction (ROBDD).
static BDDNode* unique_table;
static int unique_table_size;
static int unique_table_capacity;

// Hash table for fast lookups in the unique table.
#define U_TABLE_HASH_CAPACITY (1 << 22) // 4,194,304 slots
static int* u_table_hash;

// Computed table (cache) for the 'apply' operation.
typedef struct {
    uint64_t key; // (op << 62) | (u1_idx << 31) | u2_idx
    int res;
} ComputedEntry;

#define C_TABLE_HASH_CAPACITY (1 << 23) // 8,388,608 slots
static ComputedEntry* c_table_hash;

#define OP_AND 0
#define OP_OR  1

static long long final_result;

// Forward declarations for core BDD functions
int mk(int level, int low, int high);
int apply(int op, int u1_idx, int u2_idx);

// --- Hashing Functions ---
static inline uint32_t hash_u_table(int level, int low, int high) {
    uint32_t hash = 5381;
    hash = ((hash << 5) + hash) ^ level;
    hash = ((hash << 5) + hash) ^ low;
    hash = ((hash << 5) + hash) ^ high;
    return hash & (U_TABLE_HASH_CAPACITY - 1);
}

static inline uint64_t make_c_table_key(int op, int u1, int u2) {
    if ((op == OP_AND || op == OP_OR) && u1 > u2) {
        int temp = u1; u1 = u2; u2 = temp;
    }
    return ((uint64_t)op << 62) | ((uint64_t)u1 << 31) | (uint64_t)u2;
}

static inline uint32_t hash_c_table(uint64_t key) {
    key = (~key) + (key << 21);
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8);
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4);
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key & (C_TABLE_HASH_CAPACITY - 1);
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char* argv[]) {
    NUM_VARIABLES = atoi(argv[1]);
    NUM_CLAUSES = atoi(argv[2]);
    VAR_ORDERING_HEURISTIC = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // Setup variable ordering
    var_order = (int*)malloc(NUM_VARIABLES * sizeof(int));
    var_level = (int*)malloc(NUM_VARIABLES * sizeof(int));
    for (int i = 0; i < NUM_VARIABLES; i++) var_order[i] = i;

    if (VAR_ORDERING_HEURISTIC == 1) { // Random order shuffle
        for (int i = NUM_VARIABLES - 1; i > 0; i--) {
            int j = mt_rand() % (i + 1);
            int temp = var_order[i]; var_order[i] = var_order[j]; var_order[j] = temp;
        }
    }
    for (int i = 0; i < NUM_VARIABLES; i++) var_level[var_order[i]] = i;

    // Generate random 3-SAT formula
    clauses = (int**)malloc(NUM_CLAUSES * sizeof(int*));
    char* used_vars = (char*)malloc(NUM_VARIABLES);
    for (int i = 0; i < NUM_CLAUSES; i++) {
        clauses[i] = (int*)malloc(CLAUSE_SIZE * sizeof(int));
        memset(used_vars, 0, NUM_VARIABLES);
        for (int j = 0; j < CLAUSE_SIZE; j++) {
            int var;
            do {
                var = mt_rand() % NUM_VARIABLES;
            } while (used_vars[var]);
            used_vars[var] = 1;
            clauses[i][j] = (var + 1) * ((mt_rand() % 2) ? 1 : -1);
        }
    }
    free(used_vars);

    // Initialize BDD system
    unique_table_capacity = 1 << 21; // 2,097,152 nodes
    unique_table = (BDDNode*)malloc(unique_table_capacity * sizeof(BDDNode));
    u_table_hash = (int*)malloc(U_TABLE_HASH_CAPACITY * sizeof(int));
    memset(u_table_hash, -1, U_TABLE_HASH_CAPACITY * sizeof(int));

    c_table_hash = (ComputedEntry*)calloc(C_TABLE_HASH_CAPACITY, sizeof(ComputedEntry));
    
    // Create terminal nodes (FALSE and TRUE)
    unique_table[BDD_FALSE].level = NUM_VARIABLES; // Max level to ensure they are terminals
    unique_table[BDD_FALSE].low = BDD_FALSE;
    unique_table[BDD_FALSE].high = BDD_FALSE;
    unique_table[BDD_TRUE].level = NUM_VARIABLES;
    unique_table[BDD_TRUE].low = BDD_TRUE;
    unique_table[BDD_TRUE].high = BDD_TRUE;
    unique_table_size = 2;
}

void cleanup() {
    for (int i = 0; i < NUM_CLAUSES; i++) free(clauses[i]);
    free(clauses);
    free(var_order);
    free(var_level);
    free(unique_table);
    free(u_table_hash);
    free(c_table_hash);
}

// Create a BDD node, ensuring uniqueness.
int mk(int level, int low, int high) {
    if (low == high) return low;

    uint32_t index = hash_u_table(level, low, high);
    int node_idx;
    while ((node_idx = u_table_hash[index]) != -1) {
        if (unique_table[node_idx].level == level &&
            unique_table[node_idx].low == low &&
            unique_table[node_idx].high == high) {
            return node_idx;
        }
        index = (index + 1) & (U_TABLE_HASH_CAPACITY - 1); // Linear probe
    }
    
    if (unique_table_size >= unique_table_capacity) {
        fprintf(stderr, "FATAL: Unique table capacity exceeded.\n");
        exit(1);
    }
    int new_id = unique_table_size++;
    unique_table[new_id].level = level;
    unique_table[new_id].low = low;
    unique_table[new_id].high = high;
    
    u_table_hash[index] = new_id;
    return new_id;
}

// Apply a logical operation (AND/OR) to two BDDs.
int apply(int op, int u1_idx, int u2_idx) {
    // Terminal cases for short-circuiting
    if (op == OP_AND) {
        if (u1_idx == BDD_FALSE || u2_idx == BDD_FALSE) return BDD_FALSE;
        if (u1_idx == BDD_TRUE) return u2_idx;
        if (u2_idx == BDD_TRUE) return u1_idx;
    } else { // OP_OR
        if (u1_idx == BDD_TRUE || u2_idx == BDD_TRUE) return BDD_TRUE;
        if (u1_idx == BDD_FALSE) return u2_idx;
        if (u2_idx == BDD_FALSE) return u1_idx;
    }
    if (u1_idx == u2_idx) return u1_idx;

    uint64_t key = make_c_table_key(op, u1_idx, u2_idx);
    uint32_t index = hash_c_table(key);
    uint32_t start_index = index;
    do {
        if (c_table_hash[index].key == key) return c_table_hash[index].res;
        if (c_table_hash[index].key == 0) break; // Empty slot
        index = (index + 1) & (C_TABLE_HASH_CAPACITY - 1);
    } while(index != start_index);

    BDDNode u1 = unique_table[u1_idx];
    BDDNode u2 = unique_table[u2_idx];
    
    int top_level;
    int r_low, r_high;

    if (u1.level < u2.level) {
        top_level = u1.level;
        r_low = apply(op, u1.low, u2_idx);
        r_high = apply(op, u1.high, u2_idx);
    } else if (u2.level < u1.level) {
        top_level = u2.level;
        r_low = apply(op, u1_idx, u2.low);
        r_high = apply(op, u1_idx, u2.high);
    } else { // u1.level == u2.level
        top_level = u1.level;
        r_low = apply(op, u1.low, u2.low);
        r_high = apply(op, u1.high, u2.high);
    }

    int result = mk(top_level, r_low, r_high);
    
    c_table_hash[index].key = key;
    c_table_hash[index].res = result;

    return result;
}

void run_computation() {
    int formula_bdd = BDD_TRUE;

    for (int i = 0; i < NUM_CLAUSES; i++) {
        int clause_bdd = BDD_FALSE;
        for (int j = 0; j < CLAUSE_SIZE; j++) {
            int literal = clauses[i][j];
            int var = abs(literal) - 1;
            int level = var_level[var];
            int is_positive = literal > 0;
            
            int literal_bdd = is_positive ? mk(level, BDD_FALSE, BDD_TRUE) : mk(level, BDD_TRUE, BDD_FALSE);
            clause_bdd = apply(OP_OR, clause_bdd, literal_bdd);
        }
        formula_bdd = apply(OP_AND, formula_bdd, clause_bdd);
    }

    // The number of unique nodes created is a measure of the total complexity.
    // If the final BDD is the FALSE terminal, the formula is unsatisfiable.
    final_result = unique_table_size;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_variables num_clauses variable_ordering_heuristic seed\n", argv[0]);
        return 1;
    }

    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
