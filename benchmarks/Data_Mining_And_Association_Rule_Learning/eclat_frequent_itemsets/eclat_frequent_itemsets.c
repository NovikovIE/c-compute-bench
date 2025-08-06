#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>


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

// --- BENCHMARK DATA STRUCTURES AND GLOBALS ---

// Benchmark parameters
int NUM_TRANSACTIONS;
int NUM_UNIQUE_ITEMS;
int MIN_SUPPORT_THRESHOLD;

// Data structure for a Transaction ID list (TID-list)
typedef struct {
    int* tids;      // Dynamically allocated array of transaction IDs
    int support;    // Current number of TIDs in the list
    int capacity;   // Allocated capacity of the list
} TidList;

// Data structure for an item in the Eclat algorithm context
typedef struct {
    int item_id;
    TidList tids;
} EclatItem;

// Global array of initial frequent 1-itemsets, input to the computation
EclatItem* g_initial_frequent_items = NULL;
int g_num_initial_frequent_items = 0;

// Global accumulator for the final result
long long g_frequent_itemset_count = 0;

// --- HELPER FUNCTIONS ---

// Adds a transaction ID to a TID list, resizing if necessary
void add_tid_to_list(TidList* list, int tid) {
    if (list->support >= list->capacity) {
        list->capacity = (list->capacity == 0) ? 16 : list->capacity * 2;
        list->tids = (int*)realloc(list->tids, list->capacity * sizeof(int));
        if (!list->tids) {
            fprintf(stderr, "FATAL: realloc failed in add_tid_to_list\n");
            exit(1);
        }
    }
    list->tids[list->support++] = tid;
}

// Intersects two sorted TID lists, creating a new one.
// The new TidList's tids array must be freed by the caller.
TidList intersect_tidlists(const TidList* t1, const TidList* t2) {
    TidList result;
    // Pre-allocate memory with a reasonable guess for the capacity.
    result.capacity = (t1->support < t2->support) ? t1->support : t2->support;
    if (result.capacity > 0) {
        result.tids = (int*)malloc(result.capacity * sizeof(int));
        if (!result.tids) {
            fprintf(stderr, "FATAL: malloc failed in intersect_tidlists\n");
            exit(1);
        }
    } else {
        result.tids = NULL;
    }
    result.support = 0;

    int i = 0, j = 0;
    while (i < t1->support && j < t2->support) {
        if (t1->tids[i] < t2->tids[j]) {
            i++;
        } else if (t2->tids[j] < t1->tids[i]) {
            j++;
        } else {
            // Add the common TID to the result
            result.tids[result.support++] = t1->tids[i];
            i++; j++;
        }
    }
    return result;
}

// Forward declaration of the recursive computation function
void eclat_recursive(const EclatItem* item_set, int num_items);

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_transactions num_unique_items min_support_threshold seed\n", argv[0]);
        exit(1);
    }

    NUM_TRANSACTIONS = atoi(argv[1]);
    NUM_UNIQUE_ITEMS = atoi(argv[2]);
    MIN_SUPPORT_THRESHOLD = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // 1. Create a temporary vertical database format directly
    TidList* vertical_db = (TidList*)calloc(NUM_UNIQUE_ITEMS, sizeof(TidList));
    if (!vertical_db) {
        fprintf(stderr, "FATAL: Failed to allocate vertical_db\n");
        exit(1);
    }

    // 2. Generate transactions and populate the vertical database
    // Each transaction has a random length between 5 and 24 items
    for (int t = 0; t < NUM_TRANSACTIONS; ++t) {
        int transaction_length = 5 + (mt_rand() % 20);
        // To avoid duplicates, use a bitmap for items in the current transaction
        char* items_in_transaction = (char*)calloc(NUM_UNIQUE_ITEMS, sizeof(char));
         if (!items_in_transaction) {
            fprintf(stderr, "FATAL: Failed to allocate items_in_transaction bitmap\n");
            exit(1);
        }
        for (int i = 0; i < transaction_length; ++i) {
            int item_id = mt_rand() % NUM_UNIQUE_ITEMS;
            if (!items_in_transaction[item_id]) {
                add_tid_to_list(&vertical_db[item_id], t);
                items_in_transaction[item_id] = 1;
            }
        }
        free(items_in_transaction);
    }

    // 3. Filter for frequent 1-itemsets and create the initial set for Eclat
    int frequent_count = 0;
    for (int i = 0; i < NUM_UNIQUE_ITEMS; ++i) {
        if (vertical_db[i].support >= MIN_SUPPORT_THRESHOLD) {
            frequent_count++;
        }
    }
    
    g_num_initial_frequent_items = frequent_count;
    g_initial_frequent_items = (EclatItem*)malloc(frequent_count * sizeof(EclatItem));
    if (!g_initial_frequent_items) {
        fprintf(stderr, "FATAL: Failed to allocate g_initial_frequent_items\n");
        exit(1);
    }

    int current_frequent_item = 0;
    for (int i = 0; i < NUM_UNIQUE_ITEMS; ++i) {
        if (vertical_db[i].support >= MIN_SUPPORT_THRESHOLD) {
            g_initial_frequent_items[current_frequent_item].item_id = i;
            // Transfer ownership of TID list from temp db to global structure
            g_initial_frequent_items[current_frequent_item].tids = vertical_db[i];
            current_frequent_item++;
        } else {
            // Free TID lists of non-frequent items
            free(vertical_db[i].tids);
        }
    }

    // Free the temporary database structure (not the TID lists which were moved)
    free(vertical_db);
}

void run_computation() {
    // Initialize count with all frequent 1-itemsets
    g_frequent_itemset_count = g_num_initial_frequent_items;
    // Start the recursive mining process
    eclat_recursive(g_initial_frequent_items, g_num_initial_frequent_items);
}

// The recursive part of the Eclat algorithm
// Finds frequent itemsets by extending prefixes with items from the current item_set
void eclat_recursive(const EclatItem* item_set, int num_items) {
    if (num_items < 2) {
        return;
    }

    for (int i = 0; i < num_items; ++i) {
        const EclatItem* p1 = &item_set[i];

        // Build a conditional database for p1
        // We only need to consider items that appear after p1 to avoid duplicates
        int conditional_db_size = 0;
        int conditional_db_capacity = num_items - 1 - i;
        if (conditional_db_capacity <= 0) continue;

        EclatItem* conditional_db = (EclatItem*)malloc(conditional_db_capacity * sizeof(EclatItem));
        if (!conditional_db) {
            fprintf(stderr, "FATAL: Failed to allocate conditional_db\n");
            exit(1);
        }

        for (int j = i + 1; j < num_items; ++j) {
            const EclatItem* p2 = &item_set[j];
            TidList intersected_tids = intersect_tidlists(&p1->tids, &p2->tids);

            if (intersected_tids.support >= MIN_SUPPORT_THRESHOLD) {
                // A new frequent itemset is found. The set is {prefix..., p1->item_id, p2->item_id}
                g_frequent_itemset_count++;
                
                // Add this to the conditional database for the next recursive level
                conditional_db[conditional_db_size].item_id = p2->item_id;
                conditional_db[conditional_db_size].tids = intersected_tids; // transfer ownership
                conditional_db_size++;
            } else {
                // This intersection is not frequent, so free its resources.
                free(intersected_tids.tids);
            }
        }

        if (conditional_db_size > 0) {
            eclat_recursive(conditional_db, conditional_db_size);
        }

        // Cleanup: free the tids for all items in the conditional_db produced at this level.
        for (int k = 0; k < conditional_db_size; ++k) {
            free(conditional_db[k].tids.tids);
        }
        free(conditional_db);
    }
}

void cleanup() {
    if (g_initial_frequent_items) {
        for (int i = 0; i < g_num_initial_frequent_items; ++i) {
            free(g_initial_frequent_items[i].tids.tids);
        }
        free(g_initial_frequent_items);
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

    // Print final result to stdout
    printf("%lld\n", g_frequent_itemset_count);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
