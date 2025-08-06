#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// START of Mersenne Twister (DO NOT MODIFY)
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
// END of Mersenne Twister

// Node structure for linked lists in buckets
struct Node {
    int data;
    struct Node* next;
};

// Global data structure
typedef struct {
    int num_elements;
    int num_buckets;
    int max_val;
    int* data;
    struct Node** buckets;
    int final_result;
} BenchmarkData;

static BenchmarkData g_data;

// Helper function to insert a node into a sorted linked list
void sortedInsert(struct Node** head_ref, struct Node* new_node) {
    struct Node* current;
    if (*head_ref == NULL || (*head_ref)->data >= new_node->data) {
        new_node->next = *head_ref;
        *head_ref = new_node;
    } else {
        current = *head_ref;
        while (current->next != NULL && current->next->data < new_node->data) {
            current = current->next;
        }
        new_node->next = current->next;
        current->next = new_node;
    }
}

// Sorts a linked list using insertion sort
void insertionSort(struct Node** head_ref) {
    struct Node* sorted = NULL;
    struct Node* current = *head_ref;
    while (current != NULL) {
        struct Node* next = current->next;
        sortedInsert(&sorted, current);
        current = next;
    }
    *head_ref = sorted;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_elements> <num_buckets> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_elements = atoi(argv[1]);
    g_data.num_buckets = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (g_data.num_elements <= 0 || g_data.num_buckets <= 0) {
        fprintf(stderr, "FATAL: num_elements and num_buckets must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    g_data.max_val = 1000000;
    g_data.data = (int*)malloc(g_data.num_elements * sizeof(int));
    if (!g_data.data) {
        fprintf(stderr, "FATAL: Memory allocation failed for data array.\n");
        exit(1);
    }
    for (int i = 0; i < g_data.num_elements; i++) {
        g_data.data[i] = mt_rand() % (g_data.max_val + 1);
    }

    g_data.buckets = (struct Node**)calloc(g_data.num_buckets, sizeof(struct Node*));
    if (!g_data.buckets) {
        fprintf(stderr, "FATAL: Memory allocation failed for buckets array.\n");
        free(g_data.data);
        exit(1);
    }
}

void run_computation() {
    // 1. Scatter: Distribute elements into buckets.
    for (int i = 0; i < g_data.num_elements; i++) {
        int val = g_data.data[i];
        int bucket_idx = (int)(((long long)val * g_data.num_buckets) / (g_data.max_val + 1));

        struct Node* newNode = (struct Node*)malloc(sizeof(struct Node));
        if (!newNode) {
            fprintf(stderr, "FATAL: Memory a-llocation failed for list node.\n");
            exit(1);
        }
        newNode->data = val;
        newNode->next = g_data.buckets[bucket_idx];
        g_data.buckets[bucket_idx] = newNode;
    }

    // 2. Sort each bucket using insertion sort
    for (int i = 0; i < g_data.num_buckets; i++) {
        insertionSort(&g_data.buckets[i]);
    }

    // 3. Gather: Concatenate sorted buckets back into the original array
    int index = 0;
    for (int i = 0; i < g_data.num_buckets; i++) {
        struct Node* current = g_data.buckets[i];
        while (current != NULL) {
            g_data.data[index++] = current->data;
            current = current->next;
        }
    }

    // 4. Calculate final result to prevent dead code elimination
    g_data.final_result = 0;
    int step = g_data.num_elements > 10 ? g_data.num_elements / 10 : 1;
    for (int i = 0; i < g_data.num_elements; i += step) {
        g_data.final_result ^= g_data.data[i];
    }
}

void cleanup() {
    for (int i = 0; i < g_data.num_buckets; i++) {
        struct Node* current = g_data.buckets[i];
        while (current != NULL) {
            struct Node* temp = current;
            current = current->next;
            free(temp);
        }
    }
    free(g_data.buckets);
    free(g_data.data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    printf("%d\n", g_data.final_result);

    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
