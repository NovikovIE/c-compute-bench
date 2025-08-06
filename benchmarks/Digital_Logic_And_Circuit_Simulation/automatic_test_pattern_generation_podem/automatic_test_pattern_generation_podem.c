#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

// --- Mersenne Twister (MT19937) aken verbatim ---
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

// --- Benchmark Data Structures ---

typedef enum {
    AND, OR, NOT, NAND, NOR, XOR, XNOR, BUF, INPUT
} GateType;

typedef int LogicValue;
#define UNKNOWN -1
#define ZERO 0
#define ONE 1

typedef struct {
    GateType type;
    int input1;        // Index of the first input gate, -1 if not used
    int input2;        // Index of the second input gate, -1 if not used
    bool is_po;        // Is this a primary output?
    LogicValue value;  // Current logic value for simulation
} Gate;

typedef struct {
    int gate_index;       // Which gate's output is faulty
    LogicValue fault_type; // Stuck-at-0 or Stuck-at-1
} Fault;

// --- Global Benchmark Data ---
static int NUM_GATES;
static int NUM_TARGET_FAULTS;
static int NUM_PRIMARY_INPUTS;
static Gate* circuit = NULL;
static Fault* faults = NULL;
static LogicValue* po_good_values = NULL;
static LogicValue* po_faulty_values = NULL;

static volatile long long final_result = 0;

// --- Function Prototypes ---
LogicValue eval_gate(const Gate* g);
void simulate_circuit(bool inject_fault, const Fault* fault);

// --- Core Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_gates> <num_target_faults> <seed>\n", argv[0]);
        exit(1);
    }
    NUM_GATES = atoi(argv[1]);
    NUM_TARGET_FAULTS = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    mt_seed(seed);

    if (NUM_GATES <= 0 || NUM_TARGET_FAULTS < 0) {
        fprintf(stderr, "Parameters must be positive integers.\n");
        exit(1);
    }

    circuit = (Gate*)malloc(NUM_GATES * sizeof(Gate));
    faults = (Fault*)malloc(NUM_TARGET_FAULTS * sizeof(Fault));
    if (!circuit || !faults) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    NUM_PRIMARY_INPUTS = NUM_GATES > 10 ? NUM_GATES / 10 : 1;

    // Generate primary input gates
    for (int i = 0; i < NUM_PRIMARY_INPUTS; ++i) {
        circuit[i].type = INPUT;
        circuit[i].input1 = -1;
        circuit[i].input2 = -1;
        circuit[i].is_po = false;
    }

    // Generate combinational logic gates
    for (int i = NUM_PRIMARY_INPUTS; i < NUM_GATES; ++i) {
        circuit[i].type = (GateType)(mt_rand() % 8); // AND to BUF
        circuit[i].input1 = mt_rand() % i;
        if (circuit[i].type == NOT || circuit[i].type == BUF) {
            circuit[i].input2 = -1;
        } else {
            circuit[i].input2 = mt_rand() % i;
        }
        circuit[i].is_po = false;
    }

    // Designate some gates as primary outputs
    int num_pos = NUM_GATES > 20 ? NUM_GATES / 20 : 1;
    for (int i = 0; i < num_pos; ++i) {
        int po_idx = mt_rand() % (NUM_GATES - NUM_PRIMARY_INPUTS) + NUM_PRIMARY_INPUTS;
        circuit[po_idx].is_po = true;
    }
    po_good_values = (LogicValue*)malloc(num_pos * sizeof(LogicValue));
    po_faulty_values = (LogicValue*)malloc(num_pos * sizeof(LogicValue));
     if (!po_good_values || !po_faulty_values ) {
        fprintf(stderr, "Memory allocation failed for PO buffers\n");
        exit(1);
    }

    // Generate target faults
    for (int i = 0; i < NUM_TARGET_FAULTS; ++i) {
        faults[i].gate_index = mt_rand() % NUM_GATES;
        faults[i].fault_type = mt_rand() % 2; // 0 for Stuck-at-0, 1 for Stuck-at-1
    }
}

LogicValue eval_gate(const Gate* g) {
    LogicValue v1 = (g->input1 != -1) ? circuit[g->input1].value : UNKNOWN;
    LogicValue v2 = (g->input2 != -1) ? circuit[g->input2].value : UNKNOWN;

    switch (g->type) {
        case AND:  return (v1 == ZERO || v2 == ZERO) ? ZERO : ((v1 == ONE && v2 == ONE) ? ONE : UNKNOWN);
        case OR:   return (v1 == ONE || v2 == ONE) ? ONE : ((v1 == ZERO && v2 == ZERO) ? ZERO : UNKNOWN);
        case NAND: { LogicValue r = (v1 == ZERO || v2 == ZERO) ? ZERO : ((v1 == ONE && v2 == ONE) ? ONE : UNKNOWN); return (r == ONE) ? ZERO : ((r == ZERO) ? ONE : UNKNOWN); }
        case NOR:  { LogicValue r = (v1 == ONE || v2 == ONE) ? ONE : ((v1 == ZERO && v2 == ZERO) ? ZERO : UNKNOWN); return (r == ONE) ? ZERO : ((r == ZERO) ? ONE : UNKNOWN); }
        case XOR:  return (v1 == UNKNOWN || v2 == UNKNOWN) ? UNKNOWN : ((v1 != v2) ? ONE : ZERO);
        case XNOR: return (v1 == UNKNOWN || v2 == UNKNOWN) ? UNKNOWN : ((v1 == v2) ? ONE : ZERO);
        case NOT:  return (v1 == ONE) ? ZERO : ((v1 == ZERO) ? ONE : UNKNOWN);
        case BUF:  return v1;
        case INPUT:return g->value;
        default: return UNKNOWN;
    }
}

void simulate_circuit(bool inject_fault, const Fault* fault) {
    for (int i = 0; i < NUM_PRIMARY_INPUTS; ++i) {
        circuit[i].value = mt_rand() % 2;
    }
    for (int i = NUM_PRIMARY_INPUTS; i < NUM_GATES; ++i) {
        circuit[i].value = UNKNOWN;
    }

    for (int i = NUM_PRIMARY_INPUTS; i < NUM_GATES; ++i) {
        circuit[i].value = eval_gate(&circuit[i]);
        if (inject_fault && i == fault->gate_index) {
            circuit[i].value = fault->fault_type;
        }
    }
}

void run_computation() {
    long long detected_fault_sum = 0;

    for (int i = 0; i < NUM_TARGET_FAULTS; ++i) {
        Fault* current_fault = &faults[i];

        // Simulate backtrace by doing a random walk to a PI
        int current_gate_idx = current_fault->gate_index;
        while (circuit[current_gate_idx].type != INPUT && circuit[current_gate_idx].input1 != -1) {
            detected_fault_sum += current_gate_idx;
            if (circuit[current_gate_idx].input2 != -1 && (mt_rand() % 2)) {
                 current_gate_idx = circuit[current_gate_idx].input2;
            } else {
                 current_gate_idx = circuit[current_gate_idx].input1;
            }
        }

        // Simulate good and faulty circuit behavior to see if fault is detected
        simulate_circuit(false, NULL);
        int po_count = 0;
        for(int j=0; j<NUM_GATES; ++j) if(circuit[j].is_po) po_good_values[po_count++] = circuit[j].value;

        simulate_circuit(true, current_fault);
        po_count=0;
        for(int j=0; j<NUM_GATES; ++j) if(circuit[j].is_po) po_faulty_values[po_count++] = circuit[j].value;

        for(int j=0; j<po_count; ++j) {
            if (po_good_values[j] != po_faulty_values[j] && po_good_values[j] != UNKNOWN && po_faulty_values[j] != UNKNOWN) {
                detected_fault_sum++;
                break; // Fault detected
            }
        }
    }
    final_result = detected_fault_sum;
}

void cleanup() {
    free(circuit);
    free(faults);
    free(po_good_values);
    free(po_faulty_values);
    circuit = NULL;
    faults = NULL;
    po_good_values = NULL;
    po_faulty_values = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}