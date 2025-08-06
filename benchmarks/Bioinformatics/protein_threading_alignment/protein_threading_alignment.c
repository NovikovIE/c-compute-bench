#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

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

// Benchmark-specific definitions
#define NUM_AMINO_ACIDS 20
#define max(a,b) (((a) > (b)) ? (a) : (b))

// Global struct for benchmark data
struct {
    int protein_sequence_length;
    int num_structural_templates;
    
    // Input data
    int *protein_sequence;            // Sequence of amino acids (0-19) for our protein
    float *all_template_profiles;     // Scoring profiles for all templates
    
    // Output data
    double final_result;
} benchmark_data;

// Setup function: allocates memory and initializes data
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <protein_sequence_length> <num_structural_templates> <seed>\n", argv[0]);
        exit(1);
    }
    
    benchmark_data.protein_sequence_length = atoi(argv[1]);
    benchmark_data.num_structural_templates = atoi(argv[2]);
    uint32_t seed = (uint32_t)strtoul(argv[3], NULL, 10);
    
    mt_seed(seed);
    
    int L = benchmark_data.protein_sequence_length;
    int N = benchmark_data.num_structural_templates;
    
    // Allocate memory for the protein sequence
    benchmark_data.protein_sequence = (int*)malloc(L * sizeof(int));
    if (!benchmark_data.protein_sequence) {
        fprintf(stderr, "Failed to allocate memory for protein sequence.\n");
        exit(1);
    }
    
    // Allocate memory for all template profiles in a single block
    size_t profiles_size = (size_t)N * L * NUM_AMINO_ACIDS;
    benchmark_data.all_template_profiles = (float*)malloc(profiles_size * sizeof(float));
    if (!benchmark_data.all_template_profiles) {
        fprintf(stderr, "Failed to allocate memory for template profiles.\n");
        exit(1);
    }
    
    // Initialize the protein sequence with random amino acids
    for (int i = 0; i < L; i++) {
        benchmark_data.protein_sequence[i] = mt_rand() % NUM_AMINO_ACIDS;
    }
    
    // Initialize the template profiles with random scoring values (potentials)
    // Values are between -1.0 and 1.0 to simulate favorable/unfavorable alignments
    for (size_t i = 0; i < profiles_size; i++) {
        benchmark_data.all_template_profiles[i] = (mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }
    
    benchmark_data.final_result = 0.0;
}

// Computation function: performs the core protein threading alignment
void run_computation() {
    int L = benchmark_data.protein_sequence_length;
    int N = benchmark_data.num_structural_templates;
    const float gap_penalty = -0.5f;
    
    // Allocate space for the DP table (only two rows needed for space optimization)
    float *dp_prev_row = (float*)malloc((L + 1) * sizeof(float));
    float *dp_curr_row = (float*)malloc((L + 1) * sizeof(float));
    if (!dp_prev_row || !dp_curr_row) {
        fprintf(stderr, "Failed to allocate memory for DP table rows.\n");
        exit(1);
    }

    double max_overall_score = -1.0e30; // Initialize with a very small number

    // Iterate over each structural template
    for (int t = 0; t < N; t++) {
        // Initialize the first row of the DP table with gap penalties
        dp_prev_row[0] = 0.0f;
        for (int j = 1; j <= L; j++) {
            dp_prev_row[j] = dp_prev_row[j-1] + gap_penalty;
        }

        // Perform dynamic programming to align the sequence to the template profile
        for (int i = 1; i <= L; i++) { // i-th residue of the protein sequence
            dp_curr_row[0] = dp_prev_row[0] + gap_penalty;
            int amino_acid = benchmark_data.protein_sequence[i - 1];

            for (int j = 1; j <= L; j++) { // j-th position of the template
                // Get the substitution score from the template profile
                size_t profile_idx = (size_t)t * L * NUM_AMINO_ACIDS + (size_t)(j - 1) * NUM_AMINO_ACIDS + amino_acid;
                float match_score = benchmark_data.all_template_profiles[profile_idx];

                // Calculate scores from the three possible previous states
                float score_diag = dp_prev_row[j - 1] + match_score;
                float score_up = dp_prev_row[j] + gap_penalty;
                float score_left = dp_curr_row[j - 1] + gap_penalty;

                // The new score is the maximum of the three
                dp_curr_row[j] = max(score_diag, max(score_up, score_left));
            }
            
            // Swap pointers for the next iteration
            float *temp = dp_prev_row;
            dp_prev_row = dp_curr_row;
            dp_curr_row = temp;
        }

        // The score for this template is the final value in the DP table
        float current_template_score = dp_prev_row[L];
        
        // Update the maximum score found so far
        if (current_template_score > max_overall_score) {
            max_overall_score = current_template_score;
        }
    }

    // Free the DP table rows
    free(dp_prev_row);
    free(dp_curr_row);

    // Store the final result in the global struct
    benchmark_data.final_result = max_overall_score;
}

// Cleanup function: frees all allocated memory
void cleanup() {
    free(benchmark_data.protein_sequence);
    free(benchmark_data.all_template_profiles);
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
    printf("%.2f\n", benchmark_data.final_result);
    
    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
