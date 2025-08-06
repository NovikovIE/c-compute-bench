#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (Do Not Modify - Included Verbatim) ---
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

// --- Benchmark Specific Code ---

// Node for the Huffman Tree
struct Node {
    int data; // Using int to accommodate alphabet_size > 255
    unsigned long freq;
    struct Node *left, *right;
};

// Min-Heap for storing tree nodes
struct MinHeap {
    unsigned size;
    unsigned capacity;
    struct Node **array;
};

// Global state structure
struct {
    size_t data_size;
    int alphabet_size;
    unsigned int *input_text;
    unsigned long *freq_map;
    struct Node *huffman_root;
    struct MinHeap* min_heap; // The heap itself needs to be cleaned up
    int final_result;
} state;

// --- Utility Functions for Min-Heap and Node ---

struct Node* newNode(int data, unsigned long freq) {
    struct Node* temp = (struct Node*)malloc(sizeof(struct Node));
    if (!temp) {
        fprintf(stderr, "Memory allocation failed for new node\n");
        exit(1);
    }
    temp->left = temp->right = NULL;
    temp->data = data;
    temp->freq = freq;
    return temp;
}

struct MinHeap* createMinHeap(unsigned capacity) {
    struct MinHeap* minHeap = (struct MinHeap*)malloc(sizeof(struct MinHeap));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (struct Node**)malloc(minHeap->capacity * sizeof(struct Node*));
    if (!minHeap || !minHeap->array) {
        fprintf(stderr, "Memory allocation failed for min heap\n");
        exit(1);
    }
    return minHeap;
}

void swapNodes(struct Node** a, struct Node** b) {
    struct Node* t = *a;
    *a = *b;
    *b = t;
}

void minHeapify(struct MinHeap* minHeap, int idx) {
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;

    if (left < minHeap->size && minHeap->array[left]->freq < minHeap->array[smallest]->freq)
        smallest = left;

    if (right < minHeap->size && minHeap->array[right]->freq < minHeap->array[smallest]->freq)
        smallest = right;

    if (smallest != idx) {
        swapNodes(&minHeap->array[smallest], &minHeap->array[idx]);
        minHeapify(minHeap, smallest);
    }
}

struct Node* extractMin(struct MinHeap* minHeap) {
    struct Node* temp = minHeap->array[0];
    minHeap->array[0] = minHeap->array[minHeap->size - 1];
    --minHeap->size;
    minHeapify(minHeap, 0);
    return temp;
}

void insertMinHeap(struct MinHeap* minHeap, struct Node* node) {
    ++minHeap->size;
    int i = minHeap->size - 1;
    while (i && node->freq < minHeap->array[(i - 1) / 2]->freq) {
        minHeap->array[i] = minHeap->array[(i - 1) / 2];
        i = (i - 1) / 2;
    }
    minHeap->array[i] = node;
}

int countNodes(struct Node* root) {
    if (root == NULL) {
        return 0;
    }
    return 1 + countNodes(root->left) + countNodes(root->right);
}

void freeTree(struct Node* root) {
    if (root == NULL) {
        return;
    }
    freeTree(root->left);
    freeTree(root->right);
    free(root);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_data_size_kb> <alphabet_size> <seed>\n", argv[0]);
        exit(1);
    }

    long input_data_size_kb = atol(argv[1]);
    state.alphabet_size = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    state.data_size = input_data_size_kb * 1024;
    
    if (input_data_size_kb <= 0 || state.alphabet_size <= 1) {
        fprintf(stderr, "Invalid input sizes. Must be > 0 and alphabet > 1.\n");
        exit(1);
    }

    mt_seed(seed);

    state.input_text = (unsigned int*)malloc(state.data_size * sizeof(unsigned int));
    state.freq_map = (unsigned long*)calloc(state.alphabet_size, sizeof(unsigned long));
    if (!state.input_text || !state.freq_map) {
        fprintf(stderr, "Memory allocation failed in setup\n");
        exit(1);
    }

    // Generate random text data and build frequency map.
    // The frequency map is the input to the computation.
    for (size_t i = 0; i < state.data_size; ++i) {
        unsigned int val = mt_rand() % state.alphabet_size;
        state.input_text[i] = val;
        state.freq_map[val]++;
    }

    state.huffman_root = NULL;
    state.min_heap = NULL;
    state.final_result = 0;
}

void run_computation() {
    int non_zero_freq_count = 0;
    for (int i = 0; i < state.alphabet_size; ++i) {
        if (state.freq_map[i] > 0) {
            non_zero_freq_count++;
        }
    }
    
    if (non_zero_freq_count == 0) {
        state.final_result = 0;
        return;
    }

    state.min_heap = createMinHeap(non_zero_freq_count);

    for (int i = 0; i < state.alphabet_size; ++i) {
        if (state.freq_map[i] > 0) {
            insertMinHeap(state.min_heap, newNode(i, state.freq_map[i]));
        }
    }

    // Build the Huffman tree
    while (state.min_heap->size > 1) {
        struct Node *left = extractMin(state.min_heap);
        struct Node *right = extractMin(state.min_heap);
        
        // '$' is a special value for internal nodes
        struct Node *top = newNode('$', left->freq + right->freq);
        top->left = left;
        top->right = right;
        insertMinHeap(state.min_heap, top);
    }

    state.huffman_root = extractMin(state.min_heap);

    // To prevent D.C.E., compute a result from the tree.
    state.final_result = countNodes(state.huffman_root);
}

void cleanup() {
    free(state.input_text);
    free(state.freq_map);

    freeTree(state.huffman_root);

    if(state.min_heap) {
        free(state.min_heap->array);
        free(state.min_heap);
    }
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%d\n", state.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
