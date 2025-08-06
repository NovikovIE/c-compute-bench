#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

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

// --- BENCHMARK IMPLEMENTATION ---

#define ALPHABET_SIZE 5 // 'a','b','c','d' and '$' as terminator

// Forward declaration for Node struct
struct Node;

// Node in the Suffix Tree
typedef struct Node {
    struct Node *children[ALPHABET_SIZE];
    struct Node *suffix_link;
    int start; 
    int end; // For leaves, end = -1 to signify global end
} Node;

// --- GLOBAL STATE ---
char *text = NULL;
int text_length_param;
long long final_result = 0;

Node *root = NULL;
Node *last_new_node = NULL;
Node *active_node = NULL;
int active_edge_idx = -1;
int active_length = 0;
int remaining_suffix_count = 0;
int N; // Effective text length, including terminator
long long node_count = 0;

// --- HELPER FUNCTIONS ---

static inline int char_to_index(char c) {
    if (c == '$') return ALPHABET_SIZE - 1;
    return c - 'a';
}

// Create a new node for the suffix tree
Node *newNode(int start, int end) {
    Node *node = (Node *)malloc(sizeof(Node));
    if (!node) {
        fprintf(stderr, "FATAL: Memory allocation failed for new node.\n");
        exit(1);
    }
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        node->children[i] = NULL;
    }
    node->suffix_link = root; // Default suffix link to root
    node->start = start;
    node->end = end;
    node_count++;
    return node;
}

// Get edge length
int edge_length(Node* n, int current_pos) {
    int real_end = (n->end == -1) ? current_pos : n->end;
    return real_end - n->start + 1;
}

// The core of Ukkonen's algorithm: extend the tree with a new character
void extend_suffix_tree(int pos) {
    remaining_suffix_count++;
    last_new_node = NULL;

    while (remaining_suffix_count > 0) {
        if (active_length == 0) {
            active_edge_idx = char_to_index(text[pos]);
        }

        Node *next_node = active_node->children[active_edge_idx];

        if (next_node == NULL) { // Rule 2: New leaf creation
            active_node->children[active_edge_idx] = newNode(pos, -1);
            if (last_new_node != NULL) {
                last_new_node->suffix_link = active_node;
                last_new_node = NULL;
            }
        } else {
            int len = edge_length(next_node, pos);
            if (active_length >= len) { // Traverse to next node
                active_node = next_node;
                active_length -= len;
                active_edge_idx = char_to_index(text[pos - active_length]);
                continue; // Retry from the new active node
            }

            if (text[next_node->start + active_length] == text[pos]) { // Rule 3: Character already exists
                if (last_new_node != NULL && active_node != root) {
                    last_new_node->suffix_link = active_node;
                    last_new_node = NULL;
                }
                active_length++;
                break; // Stop this phase
            }

            // Rule 2: Split edge
            int split_end = next_node->start + active_length - 1;
            Node* split_node = newNode(next_node->start, split_end);
            active_node->children[active_edge_idx] = split_node;
            
            split_node->children[char_to_index(text[pos])] = newNode(pos, -1);
            next_node->start += active_length;
            split_node->children[char_to_index(text[next_node->start])] = next_node;
            
            if (last_new_node != NULL) {
                last_new_node->suffix_link = split_node;
            }
            last_new_node = split_node;
        }

        remaining_suffix_count--;
        if (active_node == root && active_length > 0) {
            active_length--;
            active_edge_idx = char_to_index(text[pos - remaining_suffix_count + 1]);
        } else if (active_node != root) {
            active_node = active_node->suffix_link;
        }
    }
}

// Recursive helper to free tree nodes
void free_tree_recursive(Node* node) {
    if (node == NULL) return;
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        if (node->children[i] != NULL) {
            free_tree_recursive(node->children[i]);
        }
    }
    free(node);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <text_length> <seed>\n", argv[0]);
        exit(1);
    }

    text_length_param = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (text_length_param <= 0) {
        fprintf(stderr, "FATAL: text_length must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    N = text_length_param + 1; // +1 for the terminator character `$`
    text = (char *)malloc((N + 1) * sizeof(char));
    if (text == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for text.\n");
        exit(1);
    }

    for (int i = 0; i < text_length_param; i++) {
        text[i] = 'a' + (mt_rand() % (ALPHABET_SIZE - 1)); // 'a' to 'd'
    }
    text[text_length_param] = '$'; // Terminator character
    text[N] = '\0';
}

void run_computation() {
    // Initialize the Suffix Tree structure
    root = newNode(-1, -1); 
    active_node = root;

    // Build the suffix tree character by character
    for(int i = 0; i < N; i++) {
        extend_suffix_tree(i);
    }

    final_result = node_count;
}

void cleanup() {
    free_tree_recursive(root);
    free(text);
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
    printf("%lld\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
