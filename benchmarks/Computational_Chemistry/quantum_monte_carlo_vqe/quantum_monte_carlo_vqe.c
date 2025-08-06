#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- Mersenne Twister (MT19937) Generator (verbatim) ---
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

// --- Benchmark Data Structures ---
/* 
 * This program simulates a Variational Quantum Eigensolver (VQE) calculation,
 * a hybrid quantum-classical algorithm. The goal is to find the ground state energy
 * of a given Hamiltonian. The simulation is a simplified, abstract model focusing
 * on computational kernels rather than a physically accurate quantum simulation.
 * - The Quantum Circuit is a series of parameterized gates that prepares a trial state.
 * - The Hamiltonian defines the energy function of the system.
 * - The core computation involves repeatedly preparing the state and measuring its energy,
 *   averaging the results over many "shots".
 */

// A simplified "quantum gate"
typedef struct {
    int target_qubit;
    double theta; // Parameter for the gate (e.g., a rotation angle)
} Gate;

// A term in the Hamiltonian (e.g., c * Z_i * Z_j)
typedef struct {
    double coefficient;
    int* qubits; // Array of qubit indices this term acts on
    int num_term_qubits; // Number of qubits in this term
} HamiltonianTerm;

// --- Global Pointers for Benchmark Data ---
int NUM_QUBITS;
int CIRCUIT_DEPTH;
int NUM_SHOTS;
int NUM_HAMILTONIAN_TERMS;

Gate** circuit = NULL;
HamiltonianTerm* hamiltonian = NULL;
double final_energy = 0.0; // The final computational result

// Helper to generate a random double in [0, 1)
double rand_double() {
    return (double)mt_rand() / (UINT32_MAX + 1.0);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_qubits circuit_depth num_shots num_hamiltonian_terms seed\n", argv[0]);
        exit(1);
    }

    NUM_QUBITS = atoi(argv[1]);
    CIRCUIT_DEPTH = atoi(argv[2]);
    NUM_SHOTS = atoi(argv[3]);
    NUM_HAMILTONIAN_TERMS = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if (NUM_QUBITS > 32) {
        fprintf(stderr, "Error: num_qubits cannot exceed 32 for this implementation.\n");
        exit(1);
    }

    mt_seed(seed);

    // 1. Allocate and initialize the quantum circuit
    circuit = (Gate**)malloc(CIRCUIT_DEPTH * sizeof(Gate*));
    for (int d = 0; d < CIRCUIT_DEPTH; ++d) {
        circuit[d] = (Gate*)malloc(NUM_QUBITS * sizeof(Gate));
        for (int q = 0; q < NUM_QUBITS; ++q) {
            circuit[d][q].target_qubit = q;
            // Random parameter (rotation angle) between 0 and 2*PI
            circuit[d][q].theta = rand_double() * 2.0 * M_PI;
        }
    }

    // 2. Allocate and initialize the Hamiltonian
    hamiltonian = (HamiltonianTerm*)malloc(NUM_HAMILTONIAN_TERMS * sizeof(HamiltonianTerm));
    for (int i = 0; i < NUM_HAMILTONIAN_TERMS; ++i) {
        // Random coefficient between -1.0 and 1.0
        hamiltonian[i].coefficient = rand_double() * 2.0 - 1.0;

        // Each term acts on 1 or 2 qubits (a common simplification for Pauli terms)
        hamiltonian[i].num_term_qubits = (mt_rand() % 2) + 1;
        hamiltonian[i].qubits = (int*)malloc(hamiltonian[i].num_term_qubits * sizeof(int));
        
        if (hamiltonian[i].num_term_qubits == 1) {
            hamiltonian[i].qubits[0] = mt_rand() % NUM_QUBITS;
        } else { // num_term_qubits == 2
            int q1 = mt_rand() % NUM_QUBITS;
            int q2;
            do {
                q2 = mt_rand() % NUM_QUBITS;
            } while (q1 == q2); // Ensure distinct qubits
            hamiltonian[i].qubits[0] = q1;
            hamiltonian[i].qubits[1] = q2;
        }
    }
}

void run_computation() {
    double total_energy = 0.0;
    uint32_t state_mask = (NUM_QUBITS == 32) ? UINT32_MAX : (1U << NUM_QUBITS) - 1;

    for (int shot = 0; shot < NUM_SHOTS; ++shot) {
        // Start in the |00...0> state for each shot
        uint32_t state = 0;

        // 1. Apply the quantum circuit (simplified, abstract simulation)
        // This part is designed to be computationally intensive. The operations are
        // arbitrary but deterministic, using sin() to burn CPU cycles.
        for (int d = 0; d < CIRCUIT_DEPTH; ++d) {
            for (int q = 0; q < NUM_QUBITS; ++q) {
                Gate g = circuit[d][q];
                uint32_t rotation_effect = (uint32_t)(sin(g.theta + state * 0.001) * 1e6);
                state = (state ^ rotation_effect) & state_mask;
            }
        }
        
        // 2. Measure energy of the final state against the Hamiltonian
        double shot_energy = 0.0;
        for (int i = 0; i < NUM_HAMILTONIAN_TERMS; ++i) {
            HamiltonianTerm term = hamiltonian[i];
            int parity = 0;
            
            // For a product of Pauli-Z operators, calculate the parity of the
            // relevant qubits. state is a bitmask of the measurement outcome.
            for (int k = 0; k < term.num_term_qubits; ++k) {
                int qubit_idx = term.qubits[k];
                if ((state >> qubit_idx) & 1) {
                    parity++;
                }
            }
            
            // Eigenvalue is +1 for even parity, -1 for odd parity.
            double eigenvalue = 1.0 - 2.0 * (double)(parity % 2);
            shot_energy += term.coefficient * eigenvalue;
        }
        
        total_energy += shot_energy;
    }
    
    // The final result is the average energy over all shots
    if (NUM_SHOTS > 0) {
        final_energy = total_energy / (double)NUM_SHOTS;
    }
}

void cleanup() {
    if (circuit) {
        for (int d = 0; d < CIRCUIT_DEPTH; ++d) {
            free(circuit[d]);
        }
        free(circuit);
    }
    
    if (hamiltonian) {
        for (int i = 0; i < NUM_HAMILTONIAN_TERMS; ++i) {
            free(hamiltonian[i].qubits);
        }
        free(hamiltonian);
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final calculated energy to stdout
    printf("%f\n", final_energy);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
