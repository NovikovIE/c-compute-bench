#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator --- Do Not Modify ---
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
// --- End of MT19937 code ---

// --- Benchmark Data Structures ---

// Forward declaration
struct FP_Node_s;

typedef struct {
    int item_id;
    int support;
    struct FP_Node_s* node_link_head;
} HeaderTableEntry;

typedef struct FP_Node_s {
    int item_id;
    int count;
    struct FP_Node_s* parent;
    struct FP_Node_s* children_head;
    struct FP_Node_s* next_sibling;
    struct FP_Node_s* node_link_next;
} FP_Node;

typedef struct { int id; int count; } TempItemCount;

// --- Global Variables ---

// Parameters
static int NUM_TRANSACTIONS;
static int NUM_UNIQUE_ITEMS;
static float MIN_SUPPORT_THRESHOLD;
static int MIN_SUPPORT_COUNT;

// Input Data
static int** g_transactions = NULL;
static int* g_transaction_lengths = NULL;

// FP-Growth Structures
static FP_Node* g_fp_tree_root = NULL;
static HeaderTableEntry* g_header_table = NULL;
static int g_frequent_item_count = 0;
static int* g_item_to_rank_map = NULL; // Maps original item_id to its rank

// Final Result
static long long g_final_result = 0;

// --- Helper Functions ---

int compare_item_counts(const void* a, const void* b) {
    return ((TempItemCount*)b)->count - ((TempItemCount*)a)->count;
}

int compare_transaction_items(const void* a, const void* b) {
    int rank_a = g_item_to_rank_map[*(const int*)a];
    int rank_b = g_item_to_rank_map[*(const int*)b];
    if (rank_a < rank_b) return -1;
    if (rank_a > rank_b) return 1;
    return 0;
}

FP_Node* create_fp_node(int item_id, FP_Node* parent) {
    FP_Node* node = (FP_Node*)malloc(sizeof(FP_Node));
    if (!node) { perror("malloc failed"); exit(1); }
    node->item_id = item_id;
    node->count = 1;
    node->parent = parent;
    node->children_head = NULL;
    node->next_sibling = NULL;
    node->node_link_next = NULL;
    return node;
}

void insert_tree(int* items, int num_items) {
    FP_Node* current_node = g_fp_tree_root;
    for (int i = 0; i < num_items; ++i) {
        int item_id = items[i];
        FP_Node* child = current_node->children_head;
        FP_Node* last_child = NULL;
        while (child != NULL) {
            if (child->item_id == item_id) {
                child->count++;
                current_node = child;
                goto next_item;
            }
            last_child = child;
            child = child->next_sibling;
        }

        // Item not found, create a new node
        FP_Node* new_node = create_fp_node(item_id, current_node);
        if (last_child == NULL) {
            current_node->children_head = new_node;
        } else {
            last_child->next_sibling = new_node;
        }

        // Link it in the header table
        int rank = g_item_to_rank_map[item_id];
        new_node->node_link_next = g_header_table[rank].node_link_head;
        g_header_table[rank].node_link_head = new_node;
        
        current_node = new_node;
    next_item:
        continue;
    }
}

void free_fp_tree(FP_Node* node) {
    if (node == NULL) return;
    FP_Node* child = node->children_head;
    while (child != NULL) {
        FP_Node* next = child->next_sibling;
        free_fp_tree(child);
        child = next;
    }
    free(node);
}


// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_transactions num_unique_items min_support_threshold seed\n", argv[0]);
        exit(1);
    }

    NUM_TRANSACTIONS = atoi(argv[1]);
    NUM_UNIQUE_ITEMS = atoi(argv[2]);
    MIN_SUPPORT_THRESHOLD = atof(argv[3]);
    uint32_t seed = atoi(argv[4]);

    mt_seed(seed);

    MIN_SUPPORT_COUNT = (int)(MIN_SUPPORT_THRESHOLD * NUM_TRANSACTIONS);
    if (MIN_SUPPORT_COUNT < 1) MIN_SUPPORT_COUNT = 1;

    g_transactions = (int**)malloc(NUM_TRANSACTIONS * sizeof(int*));
    g_transaction_lengths = (int*)malloc(NUM_TRANSACTIONS * sizeof(int));
    int* item_present = (int*)malloc(NUM_UNIQUE_ITEMS * sizeof(int));
    
    if (!g_transactions || !g_transaction_lengths || !item_present) {
        perror("malloc failed");
        exit(1);
    }

    for (int i = 0; i < NUM_TRANSACTIONS; ++i) {
        int len = 10 + (mt_rand() % 41); // Length from 10 to 50
        g_transaction_lengths[i] = len;
        g_transactions[i] = (int*)malloc(len * sizeof(int));
        if (!g_transactions[i]) { perror("malloc failed"); exit(1); }

        for (int k = 0; k < NUM_UNIQUE_ITEMS; ++k) item_present[k] = 0;

        int count = 0;
        while(count < len) {
            double r = (double)mt_rand() / (double)UINT32_MAX;
            int item_id = (int)(pow(r, 3.0) * NUM_UNIQUE_ITEMS);
            if (item_id >= NUM_UNIQUE_ITEMS) item_id = NUM_UNIQUE_ITEMS - 1;
            
            if (!item_present[item_id]) {
                g_transactions[i][count++] = item_id;
                item_present[item_id] = 1;
            }
        }
    }
    free(item_present);
}

void run_computation() {
    // --- Pass 1: Count item support and find frequent items ---
    int* item_counts = (int*)calloc(NUM_UNIQUE_ITEMS, sizeof(int));
    if(!item_counts) { perror("calloc failed"); exit(1); }
    
    for (int i = 0; i < NUM_TRANSACTIONS; ++i) {
        for (int j = 0; j < g_transaction_lengths[i]; ++j) {
            item_counts[g_transactions[i][j]]++;
        }
    }

    TempItemCount* temp_items = (TempItemCount*)malloc(NUM_UNIQUE_ITEMS * sizeof(TempItemCount));
    if (!temp_items) { perror("malloc failed"); exit(1); }
    
    g_frequent_item_count = 0;
    for (int i = 0; i < NUM_UNIQUE_ITEMS; ++i) {
        if (item_counts[i] >= MIN_SUPPORT_COUNT) {
            temp_items[g_frequent_item_count].id = i;
            temp_items[g_frequent_item_count].count = item_counts[i];
            g_frequent_item_count++;
        }
    }

    qsort(temp_items, g_frequent_item_count, sizeof(TempItemCount), compare_item_counts);

    g_item_to_rank_map = (int*)malloc(NUM_UNIQUE_ITEMS * sizeof(int));
    g_header_table = (HeaderTableEntry*)malloc(g_frequent_item_count * sizeof(HeaderTableEntry));
    if (!g_item_to_rank_map || !g_header_table) { perror("malloc failed"); exit(1); }

    for (int i = 0; i < NUM_UNIQUE_ITEMS; ++i) g_item_to_rank_map[i] = -1;

    for (int i = 0; i < g_frequent_item_count; ++i) {
        g_header_table[i].item_id = temp_items[i].id;
        g_header_table[i].support = temp_items[i].count;
        g_header_table[i].node_link_head = NULL;
        g_item_to_rank_map[temp_items[i].id] = i;
    }
    free(temp_items);
    free(item_counts);

    // --- Pass 2: Build FP-Tree ---
    g_fp_tree_root = create_fp_node(-1, NULL);
    int* temp_transaction_buffer = (int*)malloc(NUM_UNIQUE_ITEMS * sizeof(int));
    if(!temp_transaction_buffer) { perror("malloc failed"); exit(1); }

    for (int i = 0; i < NUM_TRANSACTIONS; ++i) {
        int frequent_count = 0;
        for (int j = 0; j < g_transaction_lengths[i]; ++j) {
            if (g_item_to_rank_map[g_transactions[i][j]] != -1) {
                temp_transaction_buffer[frequent_count++] = g_transactions[i][j];
            }
        }
        qsort(temp_transaction_buffer, frequent_count, sizeof(int), compare_transaction_items);
        insert_tree(temp_transaction_buffer, frequent_count);
    }
    free(temp_transaction_buffer);

    // --- Simplified Mining: Generate conditional pattern bases and accumulate a result ---
    g_final_result = 0;
    for (int i = 0; i < g_frequent_item_count; ++i) {
        HeaderTableEntry* header = &g_header_table[i];
        FP_Node* current_node = header->node_link_head;
        while (current_node != NULL) {
            int path_support = current_node->count;
            FP_Node* path_node = current_node->parent;
            while (path_node != NULL && path_node->item_id != -1) {
                g_final_result += path_node->item_id + path_support;
                path_node = path_node->parent;
            }
            current_node = current_node->node_link_next;
        }
    }
}

void cleanup() {
    for (int i = 0; i < NUM_TRANSACTIONS; ++i) {
        free(g_transactions[i]);
    }
    free(g_transactions);
    free(g_transaction_lengths);
    
    free_fp_tree(g_fp_tree_root);
    free(g_header_table);
    free(g_item_to_rank_map);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
