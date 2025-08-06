#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

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

// --- BENCHMARK SPECIFIC CODE ---

// Configuration
#define LEAF_LEN 32
#define MAX_INSERT_LEN 16

// Data Structures
typedef struct RopeNode {
    struct RopeNode *left, *right;
    char* str;      // NULL for internal nodes
    int weight;     // For leaves: string length. For internal: total length of left subtree.
} RopeNode;

typedef enum {
    OP_INSERT,
    OP_DELETE
} OpType;

typedef struct {
    OpType type;
    int index;
    union {
        const char* str_to_insert;
        int len_to_delete;
    };
    int insert_len;
} Operation;

// Global variables
static int INITIAL_STRING_LENGTH;
static int NUM_OPERATIONS;
static RopeNode* root = NULL;
static Operation* operations = NULL;
static char* string_arena = NULL;
static size_t arena_pos = 0;
static size_t arena_size = 0;
static char* insert_strings_buffer = NULL;
static unsigned long final_checksum = 0;

// Forward declarations
void free_rope(RopeNode* node);
RopeNode* rope_concat(RopeNode* r1, RopeNode* r2);

// Arena allocator for string fragments
char* arena_copy_str(const char* src, size_t len) {
    if (arena_pos + len + 1 > arena_size) {
        fprintf(stderr, "FATAL: String arena overflow.\n");
        exit(1);
    }
    char* dest = string_arena + arena_pos;
    memcpy(dest, src, len);
    dest[len] = '\0';
    arena_pos += len + 1;
    return dest;
}

// Destructive split: consumes `node` and produces `l_out` and `r_out`
void rope_split(RopeNode* node, int i, RopeNode** l_out, RopeNode** r_out) {
    if (!node) {
        *l_out = *r_out = NULL;
        return;
    }

    if (node->str) { // Leaf node
        if (i >= node->weight) {
            *l_out = node;
            *r_out = NULL;
        } else {
            char* s1 = arena_copy_str(node->str, i);
            char* s2 = arena_copy_str(node->str + i, node->weight - i);

            RopeNode* n1 = (RopeNode*)malloc(sizeof(RopeNode));
            n1->left = n1->right = NULL; n1->str = s1; n1->weight = i;
            
            RopeNode* n2 = (RopeNode*)malloc(sizeof(RopeNode));
            n2->left = n2->right = NULL; n2->str = s2; n2->weight = node->weight - i;

            *l_out = n1;
            *r_out = n2;
            free(node);
        }
    } else { // Internal node
        if (i < node->weight) {
            RopeNode* temp_r;
            rope_split(node->left, i, l_out, &temp_r);
            *r_out = rope_concat(temp_r, node->right);
            free(node);
        } else {
            RopeNode* temp_l;
            rope_split(node->right, i - node->weight, &temp_l, r_out);
            *l_out = rope_concat(node->left, temp_l);
            free(node);
        }
    }
}

RopeNode* rope_concat(RopeNode* r1, RopeNode* r2) {
    if (!r1) return r2;
    if (!r2) return r1;

    RopeNode* new_root = (RopeNode*)malloc(sizeof(RopeNode));
    if (!new_root) { fprintf(stderr, "FATAL: malloc failed"); exit(1); }
    
    new_root->left = r1;
    new_root->right = r2;
    new_root->str = NULL;

    int left_len = 0;
    if (r1) {
        RopeNode* temp = r1;
        while(temp && !temp->str) { // Traverse right spine to find total length
            left_len += temp->weight;
            temp = temp->right;
        }
        if (temp) left_len += temp->weight;
    }
    new_root->weight = left_len;

    return new_root;
}

void rope_insert(int i, const char* s, int s_len) {
    char* s_in_arena = arena_copy_str(s, s_len);
    RopeNode* new_leaf = (RopeNode*)malloc(sizeof(RopeNode));
    new_leaf->left = new_leaf->right = NULL;
    new_leaf->str = s_in_arena;
    new_leaf->weight = s_len;

    RopeNode *L, *R;
    rope_split(root, i, &L, &R);
    root = rope_concat(rope_concat(L, new_leaf), R);
}

void rope_delete(int i, int n) {
    RopeNode* L, *mid, *R;
    rope_split(root, i, &L, &mid);
    rope_split(mid, n, &mid, &R); 
    root = rope_concat(L, R);
    free_rope(mid);
}

void inorder_checksum(RopeNode* node, unsigned long* checksum, int* pos) {
    if (!node) return;
    if (node->str) {
        for (int i = 0; i < node->weight; ++i) {
            *checksum += (unsigned long)((*pos)++) * node->str[i];
        }
    } else {
        inorder_checksum(node->left, checksum, pos);
        inorder_checksum(node->right, checksum, pos);
    }
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <initial_string_length> <num_operations> <seed>\n", argv[0]);
        exit(1);
    }

    INITIAL_STRING_LENGTH = atoi(argv[1]);
    NUM_OPERATIONS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    // Allocate memory
    size_t arena_estimated_size = INITIAL_STRING_LENGTH + (size_t)NUM_OPERATIONS * MAX_INSERT_LEN * 2;
    string_arena = (char*)malloc(arena_estimated_size);
    operations = (Operation*)malloc(NUM_OPERATIONS * sizeof(Operation));
    size_t insert_buffer_size = (NUM_OPERATIONS / 2 + 2) * (MAX_INSERT_LEN + 1);
    insert_strings_buffer = (char*)malloc(insert_buffer_size);

    if (!string_arena || !operations || !insert_strings_buffer) {
        fprintf(stderr, "FATAL: Memory allocation failed in setup.\n");
        exit(1);
    }
    arena_size = arena_estimated_size;

    // Create initial temporary string
    char* temp_initial_string = (char*)malloc(INITIAL_STRING_LENGTH);
    if (!temp_initial_string) { exit(1); }
    for (int i = 0; i < INITIAL_STRING_LENGTH; ++i) {
        temp_initial_string[i] = (mt_rand() % 94) + 33;
    }

    // Build initial rope with small leaves
    for (int i = 0; i < INITIAL_STRING_LENGTH; i += LEAF_LEN) {
        int chunk_len = (i + LEAF_LEN > INITIAL_STRING_LENGTH) ? (INITIAL_STRING_LENGTH - i) : LEAF_LEN;
        char* chunk_in_arena = arena_copy_str(temp_initial_string + i, chunk_len);

        RopeNode* leaf = (RopeNode*)malloc(sizeof(RopeNode));
        leaf->left = leaf->right = NULL;
        leaf->str = chunk_in_arena;
        leaf->weight = chunk_len;

        root = rope_concat(root, leaf);
    }
    free(temp_initial_string);

    // Generate operations
    char* current_insert_ptr = insert_strings_buffer;
    long current_rope_len = INITIAL_STRING_LENGTH;
    for (int i = 0; i < NUM_OPERATIONS; ++i) {
        int op_choice = mt_rand() % 100;
        if (current_rope_len > MAX_INSERT_LEN && (op_choice < 50 || current_rope_len > INITIAL_STRING_LENGTH * 2)) {
            operations[i].type = OP_DELETE;
            operations[i].index = mt_rand() % current_rope_len;
            operations[i].len_to_delete = (mt_rand() % MAX_INSERT_LEN) + 1;
            if (operations[i].index + operations[i].len_to_delete > current_rope_len) {
                operations[i].len_to_delete = current_rope_len - operations[i].index;
            }
             if (operations[i].len_to_delete <= 0) { --i; continue; } // retry if invalid op
            current_rope_len -= operations[i].len_to_delete;
        } else {
            operations[i].type = OP_INSERT;
            operations[i].index = current_rope_len > 0 ? (mt_rand() % (current_rope_len + 1)) : 0;
            operations[i].insert_len = (mt_rand() % (MAX_INSERT_LEN - 1)) + 1;

            operations[i].str_to_insert = current_insert_ptr;
            for (int j = 0; j < operations[i].insert_len; ++j) {
                *current_insert_ptr++ = (mt_rand() % 94) + 33;
            }
            *current_insert_ptr++ = '\0';
            current_rope_len += operations[i].insert_len;
        }
    }
}

void run_computation() {
    for (int i = 0; i < NUM_OPERATIONS; ++i) {
        Operation* op = &operations[i];
        if (op->type == OP_INSERT) {
            rope_insert(op->index, op->str_to_insert, op->insert_len);
        } else { // OP_DELETE
            rope_delete(op->index, op->len_to_delete);
        }
    }
    
    final_checksum = 0;
    int pos = 1;
    inorder_checksum(root, &final_checksum, &pos);
}

void cleanup() {
    free(operations);
    free(insert_strings_buffer);
    free_rope(root);
    free(string_arena);
}

void free_rope(RopeNode* node) {
    if (!node) return;
    // Leaf strings are in the arena, do not free `node->str`
    free_rope(node->left);
    free_rope(node->right);
    free(node);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lu\n", final_checksum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
