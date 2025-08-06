#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
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
// --- END MERSENNE TWISTER ---

// Benchmark data and parameters
typedef struct {
    int num_taxa;
    int sequence_length;
    int num_mcmc_generations;
    char **sequences;         // The DNA sequence alignment
    double *branch_lengths;   // Represents the current tree's branch lengths
    int num_branches;
    double final_result;      // Accumulator to prevent dead code elimination
} BenchmarkData;

static BenchmarkData g_data;

// Function to calculate the log-likelihood of a tree given the data.
// This is a computationally intensive proxy for Felsenstein's pruning algorithm.
double calculate_log_likelihood(const double* branches) {
    double total_log_l = 0.0;

    // Pre-calculate a factor based on all branch lengths.
    // This simulates calculating parts of the transition probability matrix (e.g., P(t) = exp(Qt)).
    double branch_factor = 0.0;
    for (int i = 0; i < g_data.num_branches; ++i) {
        branch_factor += 1.0 - exp(-1.333333 * branches[i]); // Inspired by Jukes-Cantor model
    }
    branch_factor /= g_data.num_branches;

    // Iterate over each site in the sequence alignment
    for (int i = 0; i < g_data.sequence_length; ++i) {
        double site_score = 0.0;
        // Combine information from all taxa for this site
        for (int j = 0; j < g_data.num_taxa; ++j) {
            char base = g_data.sequences[j][i];
            // A simple but non-trivial calculation involving the base and branch factor
            site_score += log1p((double)base * branch_factor);
        }
        // The total log-likelihood is the sum of log-likelihoods for each site.
        total_log_l += log(site_score > 1e-9 ? site_score : 1e-9); // Add log of site_score, avoid log(0)
    }

    return total_log_l;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_taxa sequence_length num_mcmc_generations seed\n", argv[0]);
        exit(1);
    }

    g_data.num_taxa = atoi(argv[1]);
    g_data.sequence_length = atoi(argv[2]);
    g_data.num_mcmc_generations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    // For an unrooted binary tree with N taxa, there are 2N-3 branches.
    g_data.num_branches = g_data.num_taxa * 2 - 3;
    if (g_data.num_branches < 1) g_data.num_branches = 1;

    // Allocate and generate random DNA sequences
    g_data.sequences = (char **)malloc(g_data.num_taxa * sizeof(char *));
    if (!g_data.sequences) { perror("malloc failed"); exit(1); }
    const char bases[] = {'A', 'C', 'G', 'T'};
    for (int i = 0; i < g_data.num_taxa; ++i) {
        g_data.sequences[i] = (char *)malloc((g_data.sequence_length + 1) * sizeof(char));
        if (!g_data.sequences[i]) { perror("malloc failed"); exit(1); }
        for (int j = 0; j < g_data.sequence_length; ++j) {
            g_data.sequences[i][j] = bases[mt_rand() % 4];
        }
        g_data.sequences[i][g_data.sequence_length] = '\0';
    }

    // Allocate and initialize branch lengths for the starting tree
    g_data.branch_lengths = (double *)malloc(g_data.num_branches * sizeof(double));
    if (!g_data.branch_lengths) { perror("malloc failed"); exit(1); }
    for (int i = 0; i < g_data.num_branches; ++i) {
        g_data.branch_lengths[i] = (mt_rand() / (double)UINT32_MAX) * 0.1;
    }

    g_data.final_result = 0.0;
}

void run_computation() {
    double current_log_l = calculate_log_likelihood(g_data.branch_lengths);
    double accumulated_log_l = 0.0;

    double *proposal_branches = (double*)malloc(g_data.num_branches * sizeof(double));
    if (!proposal_branches) { perror("malloc failed"); exit(1); }

    for (int i = 0; i < g_data.num_mcmc_generations; ++i) {
        // Propose a new state by perturbing a single branch length
        memcpy(proposal_branches, g_data.branch_lengths, g_data.num_branches * sizeof(double));
        int branch_to_change = mt_rand() % g_data.num_branches;
        double perturbation = (mt_rand() / (double)UINT32_MAX - 0.5) * 0.05;
        proposal_branches[branch_to_change] += perturbation;
        if (proposal_branches[branch_to_change] < 0.0) {
            proposal_branches[branch_to_change] = -proposal_branches[branch_to_change];
        }

        // Calculate likelihood of the proposed state
        double proposal_log_l = calculate_log_likelihood(proposal_branches);

        // Metropolis-Hastings acceptance step
        double log_acceptance_ratio = proposal_log_l - current_log_l;
        double rand_val_log = log((mt_rand() + 1.0) / (double)(UINT32_MAX + 2.0)); // log(U(0,1))

        if (log_acceptance_ratio >= 0.0 || rand_val_log < log_acceptance_ratio) {
            // Accept the new state
            memcpy(g_data.branch_lengths, proposal_branches, g_data.num_branches * sizeof(double));
            current_log_l = proposal_log_l;
        }
        // If rejected, the state remains the same (implicitly)

        accumulated_log_l += current_log_l;
    }

    free(proposal_branches);
    g_data.final_result = accumulated_log_l;
}

void cleanup() {
    for (int i = 0; i < g_data.num_taxa; ++i) {
        free(g_data.sequences[i]);
    }
    free(g_data.sequences);
    free(g_data.branch_lengths);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%f\n", g_data.final_result);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
