#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>
#include <limits.h>

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

// Benchmark configuration
#define LOCAL_SEARCH_ITERATIONS 150000
#define MAX_PROCESSING_TIME 100

#define MAX(a, b) ((a) > (b) ? (a) : (b))

// Represents a single operation (job_id, op_idx)
typedef struct {
    int job_id;
    int op_idx; // Index of operation within the job
} Operation;

typedef struct {
    int num_jobs;
    int num_machines;
    int **processing_times;   // [job][op_idx] -> time
    int **machine_sequence;   // [job][op_idx] -> machine_id
    int **op_map_to_machine;  // [job][machine_id] -> op_idx
    int **schedule;           // [machine_id][order_idx] -> job_id
    long final_makespan;      // Result
} BenchmarkData;

static BenchmarkData *g_data;

// Function prototypes for helpers
static long calculate_makespan();
static void shuffle_array(int *array, int n);

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_jobs> <num_machines> <seed>\n", argv[0]);
        exit(1);
    }

    g_data = (BenchmarkData *)malloc(sizeof(BenchmarkData));
    if (!g_data) { perror("Failed to allocate benchmark data"); exit(1); }

    g_data->num_jobs = atoi(argv[1]);
    g_data->num_machines = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    // Allocate memory
    g_data->processing_times = (int **)malloc(g_data->num_jobs * sizeof(int *));
    g_data->machine_sequence = (int **)malloc(g_data->num_jobs * sizeof(int *));
    g_data->op_map_to_machine = (int **)malloc(g_data->num_jobs * sizeof(int *));
    g_data->schedule = (int **)malloc(g_data->num_machines * sizeof(int *));

    for (int i = 0; i < g_data->num_jobs; i++) {
        g_data->processing_times[i] = (int *)malloc(g_data->num_machines * sizeof(int));
        g_data->machine_sequence[i] = (int *)malloc(g_data->num_machines * sizeof(int));
        g_data->op_map_to_machine[i] = (int *)malloc(g_data->num_machines * sizeof(int));
    }
    for (int i = 0; i < g_data->num_machines; i++) {
        g_data->schedule[i] = (int *)malloc(g_data->num_jobs * sizeof(int));
    }

    // Generate data
    int *machine_perm_buffer = (int *)malloc(g_data->num_machines * sizeof(int));
    for (int i = 0; i < g_data->num_jobs; i++) {
        for (int j = 0; j < g_data->num_machines; j++) {
            g_data->processing_times[i][j] = (mt_rand() % MAX_PROCESSING_TIME) + 1;
        }

        for (int j = 0; j < g_data->num_machines; j++) machine_perm_buffer[j] = j;
        shuffle_array(machine_perm_buffer, g_data->num_machines);
        for (int j = 0; j < g_data->num_machines; j++) {
            g_data->machine_sequence[i][j] = machine_perm_buffer[j];
            g_data->op_map_to_machine[i][machine_perm_buffer[j]] = j;
        }
    }
    free(machine_perm_buffer);

    for (int i = 0; i < g_data->num_machines; i++) {
        for (int j = 0; j < g_data->num_jobs; j++) {
            g_data->schedule[i][j] = j;
        }
    }
}

void run_computation() {
    long best_makespan = calculate_makespan();

    for (int i = 0; i < LOCAL_SEARCH_ITERATIONS; ++i) {
        int m_swap = mt_rand() % g_data->num_machines;
        int j_idx_swap = mt_rand() % (g_data->num_jobs - 1);

        int job1 = g_data->schedule[m_swap][j_idx_swap];
        int job2 = g_data->schedule[m_swap][j_idx_swap + 1];
        g_data->schedule[m_swap][j_idx_swap] = job2;
        g_data->schedule[m_swap][j_idx_swap + 1] = job1;

        long current_makespan = calculate_makespan();

        if (current_makespan < best_makespan) {
            best_makespan = current_makespan;
        } else {
            g_data->schedule[m_swap][j_idx_swap] = job1;
            g_data->schedule[m_swap][j_idx_swap + 1] = job2;
        }
    }
    g_data->final_makespan = best_makespan;
}

void cleanup() {
    for (int i = 0; i < g_data->num_jobs; i++) {
        free(g_data->processing_times[i]);
        free(g_data->machine_sequence[i]);
        free(g_data->op_map_to_machine[i]);
    }
    for (int i = 0; i < g_data->num_machines; i++) {
        free(g_data->schedule[i]);
    }
    free(g_data->processing_times);
    free(g_data->machine_sequence);
    free(g_data->op_map_to_machine);
    free(g_data->schedule);
    free(g_data);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%ld\n", g_data->final_makespan);
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}

static void shuffle_array(int *array, int n) {
    for (int i = n - 1; i > 0; i--) {
        int j = mt_rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

static long calculate_makespan() {
    int nj = g_data->num_jobs;
    int nm = g_data->num_machines;
    int total_ops = nj * nm;

    long *completion_times = (long *)calloc(total_ops, sizeof(long));
    int *in_degree = (int *)calloc(total_ops, sizeof(int));
    Operation *queue = (Operation *)malloc(total_ops * sizeof(Operation));
    int queue_head = 0, queue_tail = 0;

    // Build graph and find initial ready nodes (in_degree == 0)
    for (int j = 0; j < nj; ++j) {
        for (int op = 0; op < nm; ++op) {
            int op_id = j * nm + op;
            int machine = g_data->machine_sequence[j][op];

            // Job predecessor
            if (op > 0) in_degree[op_id]++;

            // Machine predecessor
            int job_order_idx = -1;
            for (int k = 0; k < nj; ++k) {
                 if (g_data->schedule[machine][k] == j) {
                    job_order_idx = k;
                    break;
                 }
            }
            if (job_order_idx > 0) in_degree[op_id]++;

            if (in_degree[op_id] == 0) {
                queue[queue_tail].job_id = j;
                queue[queue_tail].op_idx = op;
                queue_tail++;
            }
        }
    }

    // Process nodes in topological order
    while (queue_head < queue_tail) {
        Operation u = queue[queue_head++];
        int j = u.job_id;
        int op = u.op_idx;
        int op_id = j * nm + op;
        int machine = g_data->machine_sequence[j][op];

        long job_pred_c = (op > 0) ? completion_times[op_id - 1] : 0;

        long machine_pred_c = 0;
        int job_order_idx = -1;
        for (int k = 0; k < nj; ++k) if (g_data->schedule[machine][k] == j) {job_order_idx = k; break;}
        if (job_order_idx > 0) {
            int pred_j = g_data->schedule[machine][job_order_idx - 1];
            int pred_op = g_data->op_map_to_machine[pred_j][machine];
            machine_pred_c = completion_times[pred_j * nm + pred_op];
        }

        completion_times[op_id] = MAX(job_pred_c, machine_pred_c) + g_data->processing_times[j][op];

        // Update successors
        // Job successor
        if (op < nm - 1) {
            int next_op_id = op_id + 1;
            if (--in_degree[next_op_id] == 0) {
                queue[queue_tail++] = (Operation){j, op + 1};
            }
        }
        // Machine successor
        if (job_order_idx < nj - 1) {
            int next_j = g_data->schedule[machine][job_order_idx + 1];
            int next_op = g_data->op_map_to_machine[next_j][machine];
            int next_op_id = next_j * nm + next_op;
            if (--in_degree[next_op_id] == 0) {
                queue[queue_tail++] = (Operation){next_j, next_op};
            }
        }
    }

    long makespan = 0;
    for (int i = 0; i < total_ops; i++) {
        makespan = MAX(makespan, completion_times[i]);
    }

    free(completion_times);
    free(in_degree);
    free(queue);

    return makespan;
}
