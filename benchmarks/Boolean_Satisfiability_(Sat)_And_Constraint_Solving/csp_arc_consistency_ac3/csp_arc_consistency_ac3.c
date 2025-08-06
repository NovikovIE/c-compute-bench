#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>

// --- BEGIN: Mersenne Twister (MT19937) --- Do Not Modify ---
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
// --- END: Mersenne Twister (MT19937) ---

// --- Benchmark-specific data structures and globals ---
typedef struct {
    int modulus;
    int remainder;
} Constraint;

// Parameters
int num_variables;
int domain_size;
int num_constraints;

// CSP data structures
char** domains; // Using char as boolean (0/1)
int* domain_current_sizes;
Constraint*** constraint_matrix; // Pointer matrix for efficient lookup

// Constraint graph for neighbor finding
int** neighbors;
int* neighbor_counts;

// Arc queue for AC-3
int (*arc_queue)[2];
long long queue_head = 0;
long long queue_tail = 0;
long long queue_capacity;

// Final result to prevent dead code elimination
int final_result;

// Forward declaration for helper function
bool revise(int xi_idx, int xj_idx);

// --- Benchmark Functions: setup, computation, cleanup ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_variables domain_size_per_variable num_constraints seed\n", argv[0]);
        exit(1);
    }

    num_variables = atoi(argv[1]);
    domain_size = atoi(argv[2]);
    num_constraints = atoi(argv[3]);
    mt_seed(atoi(argv[4]));

    // Allocate domains and initialize
    domains = (char**)malloc(num_variables * sizeof(char*));
    domain_current_sizes = (int*)malloc(num_variables * sizeof(int));
    for (int i = 0; i < num_variables; ++i) {
        domains[i] = (char*)malloc(domain_size * sizeof(char));
        domain_current_sizes[i] = domain_size;
        for (int j = 0; j < domain_size; ++j) {
            domains[i][j] = 1; // 1 represents 'true' (in domain)
        }
    }

    // Allocate and initialize constraint matrix for O(1) lookup
    constraint_matrix = (Constraint***)malloc(num_variables * sizeof(Constraint**));
    for (int i = 0; i < num_variables; ++i) {
        constraint_matrix[i] = (Constraint**)calloc(num_variables, sizeof(Constraint*));
    }

    // Generate constraints
    int created_constraints = 0;
    neighbor_counts = (int*)calloc(num_variables, sizeof(int));
    while (created_constraints < num_constraints && num_variables > 1) {
        int u = mt_rand() % num_variables;
        int v = mt_rand() % num_variables;

        if (u == v || constraint_matrix[u][v] != NULL) {
            continue;
        }

        Constraint* c = (Constraint*)malloc(sizeof(Constraint));
        c->modulus = (mt_rand() % (domain_size / 2)) + 2; // Avoid 0, 1
        c->remainder = mt_rand() % c->modulus;

        constraint_matrix[u][v] = c;
        constraint_matrix[v][u] = c;

        neighbor_counts[u]++;
        neighbor_counts[v]++;
        created_constraints++;
    }

    // Build neighbor lists from the constraint matrix
    neighbors = (int**)malloc(num_variables * sizeof(int*));
    int* current_neighbor_idx = (int*)calloc(num_variables, sizeof(int));
    for(int i = 0; i < num_variables; ++i) {
        neighbors[i] = (int*)malloc(neighbor_counts[i] * sizeof(int));
    }
    for (int i = 0; i < num_variables; ++i) {
        for (int j = i + 1; j < num_variables; ++j) {
            if (constraint_matrix[i][j] != NULL) {
                neighbors[i][current_neighbor_idx[i]++] = j;
                neighbors[j][current_neighbor_idx[j]++] = i;
            }
        }
    }
    free(current_neighbor_idx);

    // Initialize the arc queue
    // A loose but safe upper bound on insertions is 2*C*(D+1)
    queue_capacity = 2LL * created_constraints * (domain_size + 1);
    arc_queue = malloc(queue_capacity * sizeof(*arc_queue));
    for (int i = 0; i < num_variables; ++i) {
        for (int j = 0; j < num_variables; ++j) {
            if (i != j && constraint_matrix[i][j] != NULL) {
                arc_queue[queue_tail][0] = i;
                arc_queue[queue_tail][1] = j;
                queue_tail++;
            }
        }
    }
}

bool revise(int xi_idx, int xj_idx) {
    bool removed = false;
    Constraint* c = constraint_matrix[xi_idx][xj_idx];

    for (int i = 0; i < domain_size; ++i) {
        if (domains[xi_idx][i]) { // For each value x in D(xi)
            bool has_support = false;
            for (int j = 0; j < domain_size; ++j) {
                if (domains[xj_idx][j]) { // For each value y in D(xj)
                    // Check if (i, j) satisfies the constraint
                    if (((i + j) % c->modulus) != c->remainder) {
                        has_support = true;
                        break; 
                    }
                }
            }
            if (!has_support) {
                domains[xi_idx][i] = 0; // Remove x from D(xi)
                domain_current_sizes[xi_idx]--;
                removed = true;
            }
        }
    }
    return removed;
}

void run_computation() {
    while (queue_head < queue_tail) {
        int xi = arc_queue[queue_head][0];
        int xj = arc_queue[queue_head][1];
        queue_head++;

        if (revise(xi, xj)) {
            if (domain_current_sizes[xi] == 0) {
                final_result = 0; // Inconsistent
                return;
            }
            for (int k = 0; k < neighbor_counts[xi]; ++k) {
                int xk = neighbors[xi][k];
                if (xk != xj) {
                    if (queue_tail >= queue_capacity) {
                        // This should not happen with proper capacity planning
                        fprintf(stderr, "FATAL: Queue overflow!\n"); exit(1);
                    }
                    arc_queue[queue_tail][0] = xk;
                    arc_queue[queue_tail][1] = xi;
                    queue_tail++;
                }
            }
        }
    }

    // If consistent, calculate a result based on remaining domain values
    int sum = 0;
    for (int i = 0; i < num_variables; ++i) {
        sum += domain_current_sizes[i];
    }
    final_result = sum;
}

void cleanup() {
    for (int i = 0; i < num_variables; ++i) {
        for (int j = i + 1; j < num_variables; ++j) {
            if (constraint_matrix[i][j] != NULL) {
                free(constraint_matrix[i][j]);
            }
        }
        free(constraint_matrix[i]);
        free(domains[i]);
        free(neighbors[i]);
    }
    free(constraint_matrix);
    free(domains);
    free(domain_current_sizes);
    free(neighbors);
    free(neighbor_counts);
    free(arc_queue);
}

// --- Main function with timing ---

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
