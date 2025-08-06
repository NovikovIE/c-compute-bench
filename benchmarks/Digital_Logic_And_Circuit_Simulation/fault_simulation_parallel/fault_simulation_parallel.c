#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>

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

// --- BENCHMARK DATA STRUCTURES ---
typedef enum {
    INPUT, AND, OR, NOT, XOR
} GateType;

typedef struct {
    GateType type;
    int num_inputs;
    int input1;
    int input2;
} Gate;

typedef struct {
    int gate_index;
    int stuck_at_value; // 0 or 1
} Fault;

// Global variables for benchmark data
static int num_gates;
static int num_faults;
static int num_test_vectors;
static int num_primary_inputs;
static int final_result;

static Gate* circuit;
static Fault* faults;
static uint8_t** test_vectors;
static uint64_t* gate_values; // For simulation

// Utility for popcount
static int countSetBits(uint64_t n) {
    int count = 0;
    while (n > 0) {
        n &= (n - 1);
        count++;
    }
    return count;
}

// --- BENCHMARK FUNCTIONS ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_gates num_faults num_test_vectors num_primary_inputs seed\n", argv[0]);
        exit(1);
    }

    num_gates = atoi(argv[1]);
    num_faults = atoi(argv[2]);
    num_test_vectors = atoi(argv[3]);
    num_primary_inputs = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if (num_primary_inputs >= num_gates) {
        fprintf(stderr, "Error: num_primary_inputs must be less than num_gates.\n");
        exit(1);
    }
    
    mt_seed(seed);

    // Allocate memory
    circuit = (Gate*)malloc(num_gates * sizeof(Gate));
    faults = (Fault*)malloc(num_faults * sizeof(Fault));
    test_vectors = (uint8_t**)malloc(num_test_vectors * sizeof(uint8_t*));
    for (int i = 0; i < num_test_vectors; i++) {
        test_vectors[i] = (uint8_t*)malloc(num_primary_inputs * sizeof(uint8_t));
    }
    gate_values = (uint64_t*)malloc(num_gates * sizeof(uint64_t));

    // Generate circuit (topologically sorted by construction)
    for (int i = 0; i < num_gates; i++) {
        if (i < num_primary_inputs) {
            circuit[i].type = INPUT;
            circuit[i].num_inputs = 0;
        } else {
            circuit[i].type = (mt_rand() % 4) + 1; // AND, OR, NOT, XOR
            if (circuit[i].type == NOT) {
                circuit[i].num_inputs = 1;
                circuit[i].input1 = mt_rand() % i;
            } else {
                circuit[i].num_inputs = 2;
                circuit[i].input1 = mt_rand() % i;
                circuit[i].input2 = mt_rand() % i;
            }
        }
    }

    // Generate faults (stuck-at-0 or stuck-at-1)
    for (int i = 0; i < num_faults; i++) {
        faults[i].gate_index = num_primary_inputs + (mt_rand() % (num_gates - num_primary_inputs));
        faults[i].stuck_at_value = mt_rand() % 2;
    }

    // Generate random test vectors
    for (int i = 0; i < num_test_vectors; i++) {
        for (int j = 0; j < num_primary_inputs; j++) {
            test_vectors[i][j] = mt_rand() % 2;
        }
    }
}

void run_computation() {
    int total_detected_count = 0;
    const int faults_per_pack = 63;

    // Process faults in packs of 63 (to fit in a 64-bit word with the good circuit)
    for (int pack_start = 0; pack_start < num_faults; pack_start += faults_per_pack) {
        int num_in_pack = (num_faults - pack_start < faults_per_pack) ? (num_faults - pack_start) : faults_per_pack;

        // Pre-calculate fault masks for the current pack for all gates
        uint64_t* stuck_at_1_masks = (uint64_t*)calloc(num_gates, sizeof(uint64_t));
        uint64_t* stuck_at_0_masks = (uint64_t*)calloc(num_gates, sizeof(uint64_t));

        for (int k = 0; k < num_in_pack; k++) {
            const Fault* f = &faults[pack_start + k];
            uint64_t mask = 1ULL << (k + 1);
            if (f->stuck_at_value == 1) {
                stuck_at_1_masks[f->gate_index] |= mask;
            } else {
                stuck_at_0_masks[f->gate_index] |= mask;
            }
        }

        uint64_t pack_detected_mask = 0;

        // Simulate all test vectors for the current fault pack
        for (int i = 0; i < num_test_vectors; i++) {
            // Simulate circuit for one test vector
            for (int g = 0; g < num_gates; g++) {
                uint64_t out_val;
                if (circuit[g].type == INPUT) {
                    out_val = test_vectors[i][g] ? ~0ULL : 0ULL;
                } else {
                    uint64_t in1 = gate_values[circuit[g].input1];
                    if (circuit[g].type == NOT) {
                        out_val = ~in1;
                    } else {
                        uint64_t in2 = gate_values[circuit[g].input2];
                        switch (circuit[g].type) {
                            case AND: out_val = in1 & in2; break;
                            case OR:  out_val = in1 | in2; break;
                            case XOR: out_val = in1 ^ in2; break;
                            default:  out_val = 0; // Should not happen
                        }
                    }
                }
                // Apply faults for this gate
                gate_values[g] = (out_val & ~stuck_at_0_masks[g]) | stuck_at_1_masks[g];
            }

            // Check for detection at primary outputs
            for (int g = num_gates - num_primary_inputs; g < num_gates; g++) {
                uint64_t output_word = gate_values[g];
                uint64_t good_output_bit = output_word & 1ULL;
                uint64_t broadcast_good = good_output_bit ? ~0ULL : 0ULL;
                uint64_t differences = output_word ^ broadcast_good;
                pack_detected_mask |= differences;
            }
        }

        free(stuck_at_1_masks);
        free(stuck_at_0_masks);

        // Count detected faults in this pack (bits 1 to 63)
        total_detected_count += countSetBits(pack_detected_mask >> 1);
    }

    final_result = total_detected_count;
}

void cleanup() {
    free(circuit);
    free(faults);
    for (int i = 0; i < num_test_vectors; i++) {
        free(test_vectors[i]);
    }
    free(test_vectors);
    free(gate_values);
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
