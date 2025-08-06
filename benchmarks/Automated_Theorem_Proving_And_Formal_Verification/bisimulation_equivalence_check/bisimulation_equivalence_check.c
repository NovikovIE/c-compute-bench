#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>

// --- Mersenne Twister (MT19937) --- DO NOT MODIFY ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Specific Definitions ---
#define NUM_LABELS 8

typedef struct {
    int to;
    int label;
} Transition;

typedef struct {
    Transition* transitions;
    int count;
} State;

// --- Globals for Benchmark Data ---
int num_states_a, num_states_b;
long num_transitions_a, num_transitions_b;
State* lts_a;
State* lts_b;
Transition* all_transitions_a;
Transition* all_transitions_b;

int num_total_states;
int* partition;
int num_blocks;

// Using volatile to prevent the compiler from optimizing away the result calculation
volatile int final_result;

// --- Utility Functions ---
int compare_ints(const void* a, const void* b) {
    int val_a = *(const int*)a;
    int val_b = *(const int*)b;
    if (val_a < val_b) return -1;
    if (val_a > val_b) return 1;
    return 0;
}

State* get_state_from_combined_idx(int combined_idx) {
    if (combined_idx < num_states_a) {
        return &lts_a[combined_idx];
    }
    return &lts_b[combined_idx - num_states_a];
}

bool is_lts_a_state(int combined_idx) {
    return combined_idx < num_states_a;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_states_a num_transitions_a num_states_b num_transitions_b seed\n", argv[0]);
        exit(1);
    }
    num_states_a = atoi(argv[1]);
    num_transitions_a = atol(argv[2]);
    num_states_b = atoi(argv[3]);
    num_transitions_b = atol(argv[4]);
    mt_seed(atoi(argv[5]));

    // --- LTS A generation ---
    lts_a = (State*)calloc(num_states_a, sizeof(State));
    all_transitions_a = (Transition*)malloc(num_transitions_a * sizeof(Transition));
    int* from_states_a = (int*)malloc(num_transitions_a * sizeof(int));
    if (!lts_a || !all_transitions_a || !from_states_a) { fprintf(stderr, "alloc failed"); exit(1); }

    for (long i = 0; i < num_transitions_a; ++i) {
        int from = mt_rand() % num_states_a;
        lts_a[from].count++;
        from_states_a[i] = from;
    }

    Transition* current_trans_ptr_a = all_transitions_a;
    for (int i = 0; i < num_states_a; ++i) {
        lts_a[i].transitions = current_trans_ptr_a;
        current_trans_ptr_a += lts_a[i].count;
        lts_a[i].count = 0; // Reset to use as write index
    }

    for (long i = 0; i < num_transitions_a; ++i) {
        int from = from_states_a[i];
        int to = mt_rand() % num_states_a;
        int label = mt_rand() % NUM_LABELS;
        lts_a[from].transitions[lts_a[from].count++] = (Transition){to, label};
    }
    free(from_states_a);

    // --- LTS B generation ---
    lts_b = (State*)calloc(num_states_b, sizeof(State));
    all_transitions_b = (Transition*)malloc(num_transitions_b * sizeof(Transition));
    int* from_states_b = (int*)malloc(num_transitions_b * sizeof(int));
    if (!lts_b || !all_transitions_b || !from_states_b) { fprintf(stderr, "alloc failed"); exit(1); }
    
    for (long i = 0; i < num_transitions_b; ++i) {
        int from = mt_rand() % num_states_b;
        lts_b[from].count++;
        from_states_b[i] = from;
    }

    Transition* current_trans_ptr_b = all_transitions_b;
    for (int i = 0; i < num_states_b; ++i) {
        lts_b[i].transitions = current_trans_ptr_b;
        current_trans_ptr_b += lts_b[i].count;
        lts_b[i].count = 0; // Reset
    }

    for (long i = 0; i < num_transitions_b; ++i) {
        int from = from_states_b[i];
        int to = mt_rand() % num_states_b;
        int label = mt_rand() % NUM_LABELS;
        lts_b[from].transitions[lts_b[from].count++] = (Transition){to, label};
    }
    free(from_states_b);

    // --- Partition setup for computation ---
    num_total_states = num_states_a + num_states_b;
    partition = (int*)malloc(num_total_states * sizeof(int));
    if (!partition) { fprintf(stderr, "alloc failed"); exit(1); }
    for (int i = 0; i < num_total_states; ++i) {
        partition[i] = 0;
    }
    num_blocks = 1;
}

// --- Data structures for partition refinement ---
typedef struct {
    int* block_list;
    int count;
} SubSignature;

typedef struct {
    SubSignature parts[NUM_LABELS];
} Signature;

typedef struct {
    Signature* sig;
    int original_index;
} StateSignature;

void free_signature(Signature* sig) {
    for (int i = 0; i < NUM_LABELS; i++) {
        free(sig->parts[i].block_list);
    }
}

int compare_signatures(const Signature* a, const Signature* b) {
    for (int i = 0; i < NUM_LABELS; i++) {
        if (a->parts[i].count != b->parts[i].count) {
            return a->parts[i].count - b->parts[i].count;
        }
        for (int j = 0; j < a->parts[i].count; j++) {
            if (a->parts[i].block_list[j] != b->parts[i].block_list[j]) {
                return a->parts[i].block_list[j] - b->parts[i].block_list[j];
            }
        }
    }
    return 0;
}

int compare_state_signatures(const void* a, const void* b) {
    const StateSignature* ss_a = (const StateSignature*)a;
    const StateSignature* ss_b = (const StateSignature*)b;
    return compare_signatures(ss_a->sig, ss_b->sig);
}

void run_computation() {
    bool changed = true;
    int* temp_target_blocks = (int*)malloc(num_total_states * sizeof(int)); // Over-provisioned buffer

    while (changed) {
        changed = false;
        int next_block_id = num_blocks;
        int temp_blocks[num_blocks];
        for (int i=0; i < num_blocks; ++i) temp_blocks[i] = i; 

        for (int i_block = 0; i_block < num_blocks; ++i_block) {
            int b_id = temp_blocks[i_block];

            int states_in_block_count = 0;
            int* block_state_indices = (int*)malloc(num_total_states * sizeof(int));
            for (int i = 0; i < num_total_states; i++) {
                if (partition[i] == b_id) {
                    block_state_indices[states_in_block_count++] = i;
                }
            }

            if (states_in_block_count <= 1) {
                free(block_state_indices);
                continue;
            }

            StateSignature* state_sigs = (StateSignature*)malloc(states_in_block_count * sizeof(StateSignature));
            Signature* all_sigs = (Signature*)malloc(states_in_block_count * sizeof(Signature));

            for (int i = 0; i < states_in_block_count; i++) {
                int s_idx = block_state_indices[i];
                state_sigs[i].original_index = s_idx;
                state_sigs[i].sig = &all_sigs[i];
                State* s = get_state_from_combined_idx(s_idx);

                for (int l = 0; l < NUM_LABELS; l++) {
                    int target_count = 0;
                    for (int t = 0; t < s->count; t++) {
                        if (s->transitions[t].label == l) {
                            int target_idx = is_lts_a_state(s_idx) ? s->transitions[t].to : num_states_a + s->transitions[t].to;
                            temp_target_blocks[target_count++] = partition[target_idx];
                        }
                    }
                    if (target_count > 0) {
                        qsort(temp_target_blocks, target_count, sizeof(int), compare_ints);
                        int unique_count = (target_count > 0) ? 1 : 0;
                        for (int k = 1; k < target_count; k++) {
                            if (temp_target_blocks[k] != temp_target_blocks[k - 1]) {
                                temp_target_blocks[unique_count++] = temp_target_blocks[k];
                            }
                        }
                        state_sigs[i].sig->parts[l].count = unique_count;
                        state_sigs[i].sig->parts[l].block_list = (int*)malloc(unique_count * sizeof(int));
                        memcpy(state_sigs[i].sig->parts[l].block_list, temp_target_blocks, unique_count * sizeof(int));
                    } else {
                        state_sigs[i].sig->parts[l] = (SubSignature){NULL, 0};
                    }
                }
            }

            qsort(state_sigs, states_in_block_count, sizeof(StateSignature), compare_state_signatures);

            if (states_in_block_count > 0) {
                for (int i = 1; i < states_in_block_count; i++) {
                    if (compare_signatures(state_sigs[i - 1].sig, state_sigs[i].sig) != 0) {
                        changed = true;
                        int new_id = next_block_id++;
                        partition[state_sigs[i].original_index] = new_id;
                        for(int j = i + 1; j < states_in_block_count; j++) {
                            if (compare_signatures(state_sigs[j].sig, state_sigs[i].sig) == 0) {
                                partition[state_sigs[j].original_index] = new_id;
                            }
                        }
                    }
                }
            }

            for (int i = 0; i < states_in_block_count; i++) free_signature(state_sigs[i].sig);
            free(state_sigs);
            free(all_sigs);
            free(block_state_indices);
        }
        if (changed) num_blocks = next_block_id;
    }

    free(temp_target_blocks);

    long long sum = 0;
    for (int i = 0; i < num_total_states; i++) {
        sum += (long long)partition[i] * (i + 1);
    }
    final_result = (int)(sum % 2147483647);
}

void cleanup() {
    free(all_transitions_a);
    free(all_transitions_b);
    free(lts_a);
    free(lts_b);
    free(partition);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}