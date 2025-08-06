#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
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

// Benchmark parameters and data
int num_elements;
char initial_data_distribution[32];
int *data_array;
int final_result;

// Binary tree node structure
typedef struct Node {
    int key;
    struct Node *left, *right;
} Node;

// Helper to create a new tree node
Node* newNode(int item) {
    Node* temp = (Node*)malloc(sizeof(Node));
    if (!temp) {
        fprintf(stderr, "Failed to allocate memory for tree node.\n");
        exit(1);
    }
    temp->key = item;
    temp->left = temp->right = NULL;
    return temp;
}

// Helper to insert a new key in BST
Node* insert(Node* node, int key) {
    if (node == NULL) return newNode(key);
    if (key < node->key)
        node->left = insert(node->left, key);
    else
        node->right = insert(node->right, key);
    return node;
}

// Helper to do inorder traversal of BST and store in array
void storeSorted(Node *root, int arr[], int *i) {
    if (root != NULL) {
        storeSorted(root->left, arr, i);
        arr[(*i)++] = root->key;
        storeSorted(root->right, arr, i);
    }
}

// Helper to free tree nodes (post-order traversal)
void freeTree(Node* node) {
    if (node == NULL) return;
    freeTree(node->left);
    freeTree(node->right);
    free(node); 
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_elements initial_data_distribution seed\n", argv[0]);
        exit(1);
    }

    num_elements = atoi(argv[1]);
    strncpy(initial_data_distribution, argv[2], sizeof(initial_data_distribution) - 1);
    initial_data_distribution[sizeof(initial_data_distribution) - 1] = '\0';
    uint32_t seed = atoi(argv[3]);

    if (num_elements <= 0) {
        fprintf(stderr, "Number of elements must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    data_array = (int*)malloc(num_elements * sizeof(int));
    if (!data_array) {
        fprintf(stderr, "Failed to allocate memory for data array.\n");
        exit(1);
    }

    if (strcmp(initial_data_distribution, "random") == 0) {
        for (int i = 0; i < num_elements; i++) {
            data_array[i] = mt_rand();
        }
    } else if (strcmp(initial_data_distribution, "sorted") == 0) {
        for (int i = 0; i < num_elements; i++) {
            data_array[i] = i;
        }
    } else if (strcmp(initial_data_distribution, "reverse_sorted") == 0) {
        for (int i = 0; i < num_elements; i++) {
            data_array[i] = num_elements - 1 - i;
        }
    } else if (strcmp(initial_data_distribution, "few_unique") == 0) {
        int max_val = num_elements > 100 ? (num_elements / 10) : 10;
        for (int i = 0; i < num_elements; i++) {
            data_array[i] = mt_rand() % max_val;
        }
    } else {
        fprintf(stderr, "Invalid data distribution: %s\n", initial_data_distribution);
        free(data_array);
        exit(1);
    }
}

void run_computation() {
    if (num_elements == 0) {
        final_result = 0;
        return;
    }

    // 1. Build the BST
    Node *root = NULL;
    root = insert(root, data_array[0]);
    for (int i = 1; i < num_elements; i++) {
        insert(root, data_array[i]);
    }

    // 2. Store the in-order traversal back into the array
    int i = 0;
    storeSorted(root, data_array, &i);

    // 3. Free the tree memory
    freeTree(root);

    // 4. Calculate a final result to prevent dead code elimination
    final_result = 0;
    for (int j = 0; j < num_elements; j++) {
        final_result ^= data_array[j];
    }
}

void cleanup() {
    free(data_array);
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
