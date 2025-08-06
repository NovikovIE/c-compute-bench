#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

// --- START of Mersenne Twister ---
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
// --- END of Mersenne Twister ---

// --- START Huffman-specific structures and functions ---

// Huffman tree node
typedef struct MinHeapNode {
    unsigned char data;
    unsigned freq;
    struct MinHeapNode *left, *right;
} MinHeapNode;

// Min heap for the priority queue
typedef struct MinHeap {
    unsigned size;
    unsigned capacity;
    MinHeapNode **array;
} MinHeap;

// Global data structure for the benchmark
struct BenchmarkData {
    size_t input_size_bytes;
    int alphabet_size;
    unsigned char* input_data;
    unsigned* freq;
    char** huffman_codes;
    MinHeapNode* huffman_tree_root;
    unsigned long final_result; // Accumulated value to prevent dead code elimination
};

static struct BenchmarkData g_data;

MinHeapNode* newNode(unsigned char data, unsigned freq) {
    MinHeapNode* temp = (MinHeapNode*)malloc(sizeof(MinHeapNode));
    temp->left = temp->right = NULL;
    temp->data = data;
    temp->freq = freq;
    return temp;
}

MinHeap* createMinHeap(unsigned capacity) {
    MinHeap* minHeap = (MinHeap*)malloc(sizeof(MinHeap));
    minHeap->size = 0;
    minHeap->capacity = capacity;
    minHeap->array = (MinHeapNode**)malloc(minHeap->capacity * sizeof(MinHeapNode*));
    return minHeap;
}

void swapMinHeapNode(MinHeapNode** a, MinHeapNode** b) {
    MinHeapNode* t = *a;
    *a = *b;
    *b = t;
}

void minHeapify(MinHeap* minHeap, int idx) {
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;

    if (left < (int)minHeap->size && minHeap->array[left]->freq < minHeap->array[smallest]->freq)
        smallest = left;

    if (right < (int)minHeap->size && minHeap->array[right]->freq < minHeap->array[smallest]->freq)
        smallest = right;

    if (smallest != idx) {
        swapMinHeapNode(&minHeap->array[smallest], &minHeap->array[idx]);
        minHeapify(minHeap, smallest);
    }
}

int isSizeOne(MinHeap* minHeap) {
    return (minHeap->size == 1);
}

MinHeapNode* extractMin(MinHeap* minHeap) {
    MinHeapNode* temp = minHeap->array[0];
    minHeap->array[0] = minHeap->array[minHeap->size - 1];
    --minHeap->size;
    minHeapify(minHeap, 0);
    return temp;
}

void insertMinHeap(MinHeap* minHeap, MinHeapNode* minHeapNode) {
    ++minHeap->size;
    int i = minHeap->size - 1;
    while (i && minHeapNode->freq < minHeap->array[(i - 1) / 2]->freq) {
        minHeap->array[i] = minHeap->array[(i - 1) / 2];
        i = (i - 1) / 2;
    }
    minHeap->array[i] = minHeapNode;
}

void buildMinHeap(MinHeap* minHeap) {
    int n = minHeap->size - 1;
    for (int i = (n - 1) / 2; i >= 0; --i)
        minHeapify(minHeap, i);
}

MinHeapNode* buildHuffmanTree() {
    MinHeapNode *left, *right, *top;
    MinHeap* minHeap = createMinHeap(g_data.alphabet_size);

    for (int i = 0; i < g_data.alphabet_size; ++i) {
        if (g_data.freq[i])
            minHeap->array[minHeap->size++] = newNode(i, g_data.freq[i]);
    }
    
    if (minHeap->size == 0) {
        free(minHeap->array);
        free(minHeap);
        return NULL;
    }

    buildMinHeap(minHeap);

    while (!isSizeOne(minHeap)) {
        left = extractMin(minHeap);
        right = extractMin(minHeap);

        top = newNode('$', left->freq + right->freq);
        top->left = left;
        top->right = right;
        insertMinHeap(minHeap, top);
    }

    MinHeapNode* root = extractMin(minHeap);
    free(minHeap->array);
    free(minHeap);
    return root;
}

void generateCodes(MinHeapNode* root, int arr[], int top) {
    if (root->left) {
        arr[top] = 0;
        generateCodes(root->left, arr, top + 1);
    }
    if (root->right) {
        arr[top] = 1;
        generateCodes(root->right, arr, top + 1);
    }
    if (!(root->left) && !(root->right)) { // Leaf node
        g_data.huffman_codes[root->data] = (char*)malloc((top + 1) * sizeof(char));
        int i;
        for (i = 0; i < top; ++i)
            g_data.huffman_codes[root->data][i] = arr[i] + '0';
        g_data.huffman_codes[root->data][i] = '\0';
    }
}

void freeTree(MinHeapNode* node) {
    if (node == NULL) return;
    freeTree(node->left);
    freeTree(node->right);
    free(node);
}

// --- END Huffman-specific ---


// --- START Benchmark functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <input_size_mb> <alphabet_size> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.input_size_bytes = (size_t)atoi(argv[1]) * 1024 * 1024;
    g_data.alphabet_size = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.alphabet_size <= 0 || g_data.alphabet_size > 256) {
        fprintf(stderr, "FATAL: alphabet_size must be between 1 and 256.\n");
        exit(1);
    }
    
    mt_seed(seed);

    g_data.input_data = (unsigned char*)malloc(g_data.input_size_bytes);
    if (!g_data.input_data) { fprintf(stderr, "FATAL: Malloc failed for input_data\n"); exit(1); }
    
    g_data.freq = (unsigned*)calloc(g_data.alphabet_size, sizeof(unsigned));
    if (!g_data.freq) { fprintf(stderr, "FATAL: Calloc failed for freq table\n"); exit(1); }

    g_data.huffman_codes = (char**)calloc(g_data.alphabet_size, sizeof(char*));
    if (!g_data.huffman_codes) { fprintf(stderr, "FATAL: Calloc failed for huffman_codes\n"); exit(1); }

    for (size_t i = 0; i < g_data.input_size_bytes; ++i) {
        g_data.input_data[i] = mt_rand() % g_data.alphabet_size;
    }
    
    g_data.huffman_tree_root = NULL;
    g_data.final_result = 0;
}

void run_computation() {
    // 1. Calculate character frequencies
    for (size_t i = 0; i < g_data.input_size_bytes; ++i) {
        g_data.freq[g_data.input_data[i]]++;
    }

    // 2. Build Huffman Tree
    g_data.huffman_tree_root = buildHuffmanTree();

    // 3. Generate Huffman codes from the tree
    if (g_data.huffman_tree_root) {
        if (!g_data.huffman_tree_root->left && !g_data.huffman_tree_root->right) {
            // Edge case: only one unique character in input
            g_data.huffman_codes[g_data.huffman_tree_root->data] = (char*)malloc(2);
            strcpy(g_data.huffman_codes[g_data.huffman_tree_root->data], "0");
        } else {
            int arr[g_data.alphabet_size], top = 0;
            generateCodes(g_data.huffman_tree_root, arr, top);
        }
    }

    // 4. Calculate total size of compressed data in bits
    // This loop uses the generated codes and ensures previous work is not optimized out.
    unsigned long total_bits = 0;
    for (size_t i = 0; i < g_data.input_size_bytes; ++i) {
        char* code = g_data.huffman_codes[g_data.input_data[i]];
        if(code) {
           total_bits += strlen(code);
        }
    }
    g_data.final_result = total_bits;
}

void cleanup() {
    free(g_data.input_data);
    free(g_data.freq);
    
    if (g_data.huffman_codes) {
        for (int i = 0; i < g_data.alphabet_size; ++i) {
            if (g_data.huffman_codes[i])
                free(g_data.huffman_codes[i]);
        }
        free(g_data.huffman_codes);
    }
    
    freeTree(g_data.huffman_tree_root);
}
// --- END Benchmark functions ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%lu\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
