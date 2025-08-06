#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <float.h>

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

// --- BENCHMARK GLOBALS ---
int num_taxa;
float** dist_matrix;
int* active_nodes; // 1 if active, 0 if not
float final_result;

// --- BENCHMARK FUNCTIONS ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_taxa> <seed>\n", argv[0]);
        exit(1);
    }
    num_taxa = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);
    mt_seed(seed);

    if (num_taxa < 3) {
        fprintf(stderr, "FATAL: num_taxa must be at least 3.\n");
        exit(1);
    }

    dist_matrix = (float**)malloc(num_taxa * sizeof(float*));
    if (!dist_matrix) { fprintf(stderr, "Memory allocation failed.\n"); exit(1); }
    for (int i = 0; i < num_taxa; i++) {
        dist_matrix[i] = (float*)malloc(num_taxa * sizeof(float));
        if (!dist_matrix[i]) { fprintf(stderr, "Memory allocation failed.\n"); exit(1); }
    }

    active_nodes = (int*)malloc(num_taxa * sizeof(int));
    if (!active_nodes) { fprintf(stderr, "Memory allocation failed.\n"); exit(1); }

    for (int i = 0; i < num_taxa; i++) {
        active_nodes[i] = 1;
        dist_matrix[i][i] = 0.0f;
        for (int j = i + 1; j < num_taxa; j++) {
            float dist = ((float)mt_rand() / (float)UINT32_MAX) * 99.0f + 1.0f;
            dist_matrix[i][j] = dist;
            dist_matrix[j][i] = dist;
        }
    }
}

void run_computation() {
    int current_nodes = num_taxa;
    float total_branch_length = 0.0f;
    float* r = (float*)malloc(num_taxa * sizeof(float));
    if (!r) { exit(1); }

    while (current_nodes > 2) {
        // Step 1: Calculate sum of distances for each node
        for (int i = 0; i < num_taxa; i++) {
            if (active_nodes[i]) {
                r[i] = 0.0f;
                for (int k = 0; k < num_taxa; k++) {
                    if (active_nodes[k] && i != k) {
                        r[i] += dist_matrix[i][k];
                    }
                }
            }
        }

        // Step 2: Find pair (i, j) that minimizes Q_ij
        float min_q = FLT_MAX;
        int min_i = -1, min_j = -1;
        float n_minus_2 = (float)(current_nodes - 2);

        for (int i = 0; i < num_taxa; i++) {
            if (!active_nodes[i]) continue;
            for (int j = i + 1; j < num_taxa; j++) {
                if (!active_nodes[j]) continue;

                float q_val = n_minus_2 * dist_matrix[i][j] - r[i] - r[j];
                if (q_val < min_q) {
                    min_q = q_val;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        // Step 3: Calculate limb lengths for new node and add to total tree length
        float limb_i = 0.5f * dist_matrix[min_i][min_j] + (r[min_i] - r[min_j]) / (2.0f * n_minus_2);
        float limb_j = dist_matrix[min_i][min_j] - limb_i;
        total_branch_length += limb_i + limb_j;

        // Step 4: Update distance matrix. min_i becomes the new combined node.
        for (int k = 0; k < num_taxa; k++) {
            if (active_nodes[k] && k != min_i && k != min_j) {
                float new_dist = 0.5f * (dist_matrix[min_i][k] + dist_matrix[min_j][k] - dist_matrix[min_i][min_j]);
                dist_matrix[min_i][k] = new_dist;
                dist_matrix[k][min_i] = new_dist;
            }
        }

        // Step 5: Deactivate node min_j from consideration
        active_nodes[min_j] = 0;
        current_nodes--;
    }
    
    free(r);

    // Step 6: Add the last branch connecting the final two nodes
    int last_i = -1, last_j = -1;
    for (int i = 0; i < num_taxa; i++) {
        if (active_nodes[i]) {
            if (last_i == -1) {
                last_i = i;
            } else {
                last_j = i;
                break;
            }
        }
    }
    if (last_i != -1 && last_j != -1) {
        total_branch_length += dist_matrix[last_i][last_j];
    }
    
    final_result = total_branch_length;
}

void cleanup() {
    for (int i = 0; i < num_taxa; i++) {
        free(dist_matrix[i]);
    }
    free(dist_matrix);
    free(active_nodes);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
