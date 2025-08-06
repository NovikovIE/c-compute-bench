#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Data and Globals ---

// Description: Beat tracking finds the pulse of music. This program models it by
// identifying prominent onsets (start of musical events) and using dynamic programming
// to find the most likely sequence of beats (the "beat path") that aligns with these
// onsets, evaluating a range of tempi (beats per minute, BPM).

struct BenchmarkData {
    int num_onsets;
    double min_bpm;
    double max_bpm;
    double* onsets;         // Array of onset timestamps in seconds
    double* dp_score;       // DP table for cumulative scores
    int* backpointers;      // DP table for path reconstruction
    int final_beat_path_length; // Final result to prevent dead code elimination
};

static struct BenchmarkData g_data;

// Helper for qsort
int compare_doubles(const void *a, const void *b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_onsets min_bpm max_bpm seed\n", argv[0]);
        exit(1);
    }

    g_data.num_onsets = atoi(argv[1]);
    g_data.min_bpm = atof(argv[2]);
    g_data.max_bpm = atof(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    if (g_data.num_onsets <= 0 || g_data.min_bpm <= 0 || g_data.max_bpm <= g_data.min_bpm) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }

    g_data.onsets = (double*)malloc(g_data.num_onsets * sizeof(double));
    g_data.dp_score = (double*)malloc(g_data.num_onsets * sizeof(double));
    g_data.backpointers = (int*)malloc(g_data.num_onsets * sizeof(int));

    if (!g_data.onsets || !g_data.dp_score || !g_data.backpointers) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Generate random onset times within a 5-minute simulated track
    double max_time = 300.0;
    for (int i = 0; i < g_data.num_onsets; i++) {
        g_data.onsets[i] = ((double)mt_rand() / (double)UINT32_MAX) * max_time;
    }

    // Onsets must be sorted by time for the algorithm to work
    qsort(g_data.onsets, g_data.num_onsets, sizeof(double), compare_doubles);

    g_data.final_beat_path_length = 0;
}

void run_computation() {
    // Initialize DP tables
    for (int i = 0; i < g_data.num_onsets; i++) {
        g_data.dp_score[i] = 1.0; // Base salience score for each onset
        g_data.backpointers[i] = -1;
    }

    double mid_bpm = (g_data.min_bpm + g_data.max_bpm) / 2.0;
    double bpm_range_half = g_data.max_bpm - mid_bpm;

    // Dynamic Programming core loop
    for (int i = 1; i < g_data.num_onsets; i++) {
        double max_prev_score = -1.0/0.0; // Negative infinity
        int best_prev_j = -1;

        for (int j = 0; j < i; j++) {
            double interval = g_data.onsets[i] - g_data.onsets[j];

            // Avoid instability with very close onsets
            if (interval < 1e-5) {
                continue;
            }

            double bpm_equivalent = 60.0 / interval;
            double transition_score;

            if (bpm_equivalent < g_data.min_bpm || bpm_equivalent > g_data.max_bpm) {
                transition_score = -100.0; // Heavy penalty for out-of-range tempo
            } else {
                // A parabolic score that peaks at mid_bpm and is 0 at min/max_bpm
                double normalized_dist = (bpm_equivalent - mid_bpm) / bpm_range_half;
                transition_score = 2.0 * (1.0 - normalized_dist * normalized_dist);
            }
            
            double current_path_score = g_data.dp_score[j] + transition_score;

            if (current_path_score > max_prev_score) {
                max_prev_score = current_path_score;
                best_prev_j = j;
            }
        }

        if (best_prev_j != -1) {
            g_data.dp_score[i] += max_prev_score;
            g_data.backpointers[i] = best_prev_j;
        }
    }

    // Find the best path ending and backtrack to calculate its length
    double max_final_score = -1.0/0.0; // Negative infinity
    int last_beat_index = -1;
    for (int i = 0; i < g_data.num_onsets; i++) {
        if (g_data.dp_score[i] > max_final_score) {
            max_final_score = g_data.dp_score[i];
            last_beat_index = i;
        }
    }

    int path_len = 0;
    if (last_beat_index != -1) {
        int current_index = last_beat_index;
        while (current_index != -1) {
            path_len++;
            current_index = g_data.backpointers[current_index];
        }
    }
    g_data.final_beat_path_length = path_len;
}

void cleanup() {
    free(g_data.onsets);
    free(g_data.dp_score);
    free(g_data.backpointers);
    g_data.onsets = NULL;
    g_data.dp_score = NULL;
    g_data.backpointers = NULL;
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
    printf("%d\n", g_data.final_beat_path_length);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
