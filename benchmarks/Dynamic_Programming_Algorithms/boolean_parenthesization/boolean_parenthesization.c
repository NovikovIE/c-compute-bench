#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// --- Globals for Benchmark Data ---
int expression_length;
char *symbols;    // Array of 'T' or 'F'
char *operators;  // Array of '&', '|', '^'
long long **T_dp; // DP table for True counts
long long **F_dp; // DP table for False counts
long long final_result;

// Use a large prime for modulo arithmetic to prevent overflow
const long long MOD = 1000000007;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

// --- Benchmark Implementation ---

/**
 * @brief Parses command-line args, allocates memory, and generates input data.
 */
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <expression_length> <seed>\n", argv[0]);
        exit(1);
    }

    expression_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if (expression_length <= 1) {
       fprintf(stderr, "FATAL: expression_length must be > 1 to have operators.\n");
       exit(1);
    }

    mt_seed(seed);

    // Allocate memory on the heap
    symbols = (char *)malloc(expression_length * sizeof(char));
    operators = (char *)malloc((expression_length - 1) * sizeof(char));
    T_dp = (long long **)malloc(expression_length * sizeof(long long *));
    F_dp = (long long **)malloc(expression_length * sizeof(long long *));

    if (!symbols || !operators || !T_dp || !F_dp) {
        fprintf(stderr, "FATAL: Top-level memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < expression_length; ++i) {
        T_dp[i] = (long long *)malloc(expression_length * sizeof(long long));
        F_dp[i] = (long long *)malloc(expression_length * sizeof(long long));
        if (!T_dp[i] || !F_dp[i]) {
            fprintf(stderr, "FATAL: DP table row allocation failed.\n");
            exit(1);
        }
    }

    // Generate random expression data
    for (int i = 0; i < expression_length; ++i) {
        symbols[i] = (mt_rand() % 2 == 0) ? 'T' : 'F';
    }
    const char ops[] = {'&', '|', '^'};
    for (int i = 0; i < expression_length - 1; ++i) {
        operators[i] = ops[mt_rand() % 3];
    }
}

/**
 * @brief Executes the boolean parenthesization algorithm.
 */
void run_computation() {
    int n = expression_length;

    // Base case: subexpressions of length 1 (single symbol)
    for (int i = 0; i < n; i++) {
        if (symbols[i] == 'T') {
            T_dp[i][i] = 1;
            F_dp[i][i] = 0;
        } else {
            T_dp[i][i] = 0;
            F_dp[i][i] = 1;
        }
    }

    // Fill DP table for larger subexpressions using tabulation
    // L is the length of the subexpression chain
    for (int L = 2; L <= n; L++) {
        // i is the starting index of the subexpression
        for (int i = 0; i <= n - L; i++) {
            int j = i + L - 1; // Ending index
            T_dp[i][j] = 0;
            F_dp[i][j] = 0;

            // k is the split point (position of the operator)
            for (int k = i; k < j; k++) {
                long long lT = T_dp[i][k];
                long long lF = F_dp[i][k];
                long long rT = T_dp[k + 1][j];
                long long rF = F_dp[k + 1][j];

                long long total_left = (lT + lF);
                long long total_right = (rT + rF);
                long long total_combinations = (total_left * total_right) % MOD;

                if (operators[k] == '&') {
                    long long true_ways = (lT * rT) % MOD;
                    long long false_ways = (total_combinations - true_ways + MOD) % MOD;
                    T_dp[i][j] = (T_dp[i][j] + true_ways) % MOD;
                    F_dp[i][j] = (F_dp[i][j] + false_ways) % MOD;
                } else if (operators[k] == '|') {
                    long long false_ways = (lF * rF) % MOD;
                    long long true_ways = (total_combinations - false_ways + MOD) % MOD;
                    T_dp[i][j] = (T_dp[i][j] + true_ways) % MOD;
                    F_dp[i][j] = (F_dp[i][j] + false_ways) % MOD;
                } else if (operators[k] == '^') {
                    long long true_ways = ((lT * rF) % MOD + (lF * rT) % MOD) % MOD;
                    long long false_ways = (total_combinations - true_ways + MOD) % MOD;
                    T_dp[i][j] = (T_dp[i][j] + true_ways) % MOD;
                    F_dp[i][j] = (F_dp[i][j] + false_ways) % MOD;
                }
            }
        }
    }

    final_result = T_dp[0][n - 1];
}

/**
 * @brief Frees all memory allocated in setup_benchmark.
 */
void cleanup() {
    for (int i = 0; i < expression_length; i++) {
        free(T_dp[i]);
        free(F_dp[i]);
    }
    free(T_dp);
    free(F_dp);
    free(symbols);
    free(operators);
}

// --- Main and Timing ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    // Print the final result (number of ways to evaluate to True) to stdout
    printf("%lld\n", final_result);

    // Calculate and print the execution time to stderr
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
