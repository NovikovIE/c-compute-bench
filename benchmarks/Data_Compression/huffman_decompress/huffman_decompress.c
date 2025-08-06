#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

// --- Mersenne Twister (verbatim) ---
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
// --- end Mersenne Twister ---

// --- Benchmark Specific Code ---

// Huffman tree node structure
typedef struct HuffNode {
    unsigned char data;
    unsigned int freq;
    struct HuffNode *left, *right;
} HuffNode;

// --- Global variables for benchmark data ---
static unsigned char *compressed_data = NULL;
static size_t compressed_size_bits = 0;
static HuffNode *huffman_tree_root = NULL;
static size_t original_size = 0;

// Data produced by the computation
static unsigned char *decompressed_data = NULL;
static unsigned long long final_result = 0; // Checksum for stdout

// --- Helper functions for setup ---

// Creates a new Huffman tree node
HuffNode* newNode(unsigned char data, unsigned int freq) {
    HuffNode* temp = (HuffNode*)malloc(sizeof(HuffNode));
    if (!temp) {
        fprintf(stderr, "FATAL: Memory allocation failed for HuffNode\n");
        exit(1);
    }
    temp->left = temp->right = NULL;
    temp->data = data;
    temp->freq = freq;
    return temp;
}

// Recursively generate Huffman codes from the tree
void generateCodes(HuffNode* root, unsigned int code, int depth, unsigned int* codes, int* depths) {
    if (!root) return;

    if (!root->left && !root->right) { // Leaf node
        codes[root->data] = code;
        depths[root->data] = (depth == 0) ? 1 : depth; // Handle single-node tree
        return;
    }
    
    generateCodes(root->left, code << 1, depth + 1, codes, depths);
    generateCodes(root->right, (code << 1) | 1, depth + 1, codes, depths);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <input_size_mb> <seed>\n", argv[0]);
        exit(1);
    }

    int input_size_mb = atoi(argv[1]);
    unsigned int seed = atoi(argv[2]);

    original_size = (size_t)input_size_mb * 1024 * 1024;
    mt_seed(seed);

    if (original_size == 0) {
        huffman_tree_root = newNode(0, 0);
        return;
    }

    unsigned char* original_data = (unsigned char*)malloc(original_size);
    if (!original_data) {
        fprintf(stderr, "FATAL: Failed to allocate memory for original data.\n");
        exit(1);
    }

    // 1. Generate random original data
    for (size_t i = 0; i < original_size; ++i) {
        original_data[i] = mt_rand() % 256;
    }

    // 2. Calculate character frequencies
    unsigned int freqs[256] = {0};
    for (size_t i = 0; i < original_size; ++i) {
        freqs[original_data[i]]++;
    }

    // 3. Build the Huffman Tree
    HuffNode* nodes[256];
    int node_count = 0;
    for (int i = 0; i < 256; ++i) {
        if (freqs[i] > 0) {
            nodes[node_count++] = newNode((unsigned char)i, freqs[i]);
        }
    }
    
    if (node_count <= 1) {
        huffman_tree_root = (node_count == 1) ? nodes[0] : newNode(0,0);
    } else {
        for (int i = 0; i < node_count - 1; ++i) {
            long min1_idx = -1, min2_idx = -1;
            for (int j = 0; j < node_count; ++j) {
                if (nodes[j]) {
                    if (min1_idx == -1 || nodes[j]->freq < nodes[min1_idx]->freq) {
                        min2_idx = min1_idx;
                        min1_idx = j;
                    } else if (min2_idx == -1 || nodes[j]->freq < nodes[min2_idx]->freq) {
                        min2_idx = j;
                    }
                }
            }
            
            HuffNode* left = nodes[min1_idx];
            HuffNode* right = nodes[min2_idx];
            HuffNode* parent = newNode(0, left->freq + right->freq);
            parent->left = left;
            parent->right = right;

            nodes[min1_idx] = parent;
            nodes[min2_idx] = NULL;
        }
        for (int i = 0; i < node_count; ++i) {
            if (nodes[i]) {
                huffman_tree_root = nodes[i];
                break;
            }
        }
    }

    // 4. Generate Huffman codes from the tree
    unsigned int codes[256];
    int depths[256];
    memset(codes, 0, sizeof(codes));
    memset(depths, 0, sizeof(depths));
    generateCodes(huffman_tree_root, 0, 0, codes, depths);

    // 5. Compress the data
    compressed_size_bits = 0;
    for (size_t i = 0; i < original_size; ++i) {
        compressed_size_bits += depths[original_data[i]];
    }

    size_t compressed_size_bytes = (compressed_size_bits + 7) / 8;
    compressed_data = (unsigned char*)malloc(compressed_size_bytes);
    if (!compressed_data) {
        fprintf(stderr, "FATAL: Failed to allocate memory for compressed data.\n");
        exit(1);
    }
    memset(compressed_data, 0, compressed_size_bytes);

    size_t bit_pos = 0;
    for (size_t i = 0; i < original_size; ++i) {
        unsigned char c = original_data[i];
        int depth = depths[c];
        unsigned int code = codes[c];
        for (int j = depth - 1; j >= 0; --j) {
            if ((code >> j) & 1) {
                compressed_data[bit_pos / 8] |= (1 << (7 - (bit_pos % 8)));
            }
            bit_pos++;
        }
    }
    
    free(original_data);
}

void run_computation() {
    if (original_size == 0) {
        final_result = 0;
        return;
    }

    decompressed_data = (unsigned char*)malloc(original_size);
    if (!decompressed_data) {
        fprintf(stderr, "FATAL: Memory allocation failed for decompressed data\n");
        exit(1);
    }

    size_t bit_pos = 0;
    for (size_t i = 0; i < original_size; ++i) {
        HuffNode* current = huffman_tree_root;
        while (current->left != NULL) { // Internal nodes always have both children
            int bit = (compressed_data[bit_pos / 8] >> (7 - (bit_pos % 8))) & 1;
            bit_pos++;
            if (bit) {
                current = current->right;
            } else {
                current = current->left;
            }
        }
        decompressed_data[i] = current->data;
    }

    // Calculate a checksum to prevent dead code elimination
    unsigned long long checksum = 0;
    for (size_t i = 0; i < original_size; ++i) {
        checksum += decompressed_data[i];
    }
    final_result = checksum;
}

void free_tree(HuffNode* node) {
    if (node == NULL) {
        return;
    }
    free_tree(node->left);
    free_tree(node->right);
    free(node);
}

void cleanup() {
    if (compressed_data) {
        free(compressed_data);
    }
    if (decompressed_data) {
        free(decompressed_data);
    }
    if (huffman_tree_root) {
        free_tree(huffman_tree_root);
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
    printf("%llu\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
