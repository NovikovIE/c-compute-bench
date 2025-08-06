#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// Program: crack_propagation_simulation
// Description: Simulates crack propagation in a 1D material mesh.
// The simulation applies an increasing load over a number of steps.
// For each step, it calculates stress and checks if elements at crack tips
// have reached their fracture toughness. If so, they "fail" and the crack
// propagates to their neighbors.

// Mersenne Twister (MT19937) PRNG
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
// End of MT19937

// --- Benchmark Data Structures ---

typedef struct {
    double youngs_modulus;      // Material stiffness
    double fracture_toughness;  // Resistance to fracture
    double stress;              // Current stress on the element
    double damage;              // Damage state (0.0 = intact, 1.0 = fully cracked)
    int is_crack_tip;           // Flag: is this element at the tip of a crack?
} Element;

// Global parameters
int MESH_ELEMENTS;
int NUM_LOAD_STEPS;

// Global data arrays
Element *mesh;
int *crack_tip_indices;
int *next_crack_tip_indices;
int num_crack_tips;

// Final result accumulator
long long final_damage_metric = 0;

// Helper to generate a random double
double rand_double(double min, double max) {
    return min + (double)mt_rand() / (double)UINT32_MAX * (max - min);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <mesh_elements> <num_load_steps> <seed>\n", argv[0]);
        exit(1);
    }
    MESH_ELEMENTS = atoi(argv[1]);
    NUM_LOAD_STEPS = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    mesh = (Element *)malloc(MESH_ELEMENTS * sizeof(Element));
    crack_tip_indices = (int *)malloc(MESH_ELEMENTS * sizeof(int));
    next_crack_tip_indices = (int *)malloc(MESH_ELEMENTS * sizeof(int));

    if (!mesh || !crack_tip_indices || !next_crack_tip_indices) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize material properties and state for each element
    for (int i = 0; i < MESH_ELEMENTS; i++) {
        mesh[i].youngs_modulus = rand_double(50.0, 250.0); // GPa
        mesh[i].fracture_toughness = rand_double(15.0, 45.0); // Fictional units
        mesh[i].stress = 0.0;
        mesh[i].damage = 0.0;
        mesh[i].is_crack_tip = 0;
    }

    // Create an initial crack at the center of the mesh
    int initial_crack_pos = MESH_ELEMENTS / 2;
    if (initial_crack_pos > 0 && initial_crack_pos < MESH_ELEMENTS - 1) {
        mesh[initial_crack_pos].damage = 1.0; // This element has failed
        // Neighbors are the initial crack tips
        mesh[initial_crack_pos - 1].is_crack_tip = 1;
        mesh[initial_crack_pos + 1].is_crack_tip = 1;
        crack_tip_indices[0] = initial_crack_pos - 1;
        crack_tip_indices[1] = initial_crack_pos + 1;
        num_crack_tips = 2;
    } else {
        num_crack_tips = 0; // No crack if mesh is too small
    }
}

void run_computation() {
    const double STRESS_INTENSITY_FACTOR_GEOMETRY = 1.12; // Simplified constant
    const double STRAIN_PER_STEP = 0.00001; // Small strain applied per load step

    for (int step = 0; step < NUM_LOAD_STEPS; ++step) {
        if (num_crack_tips == 0) {
            break; // Simulation ends if there are no more cracks to propagate
        }

        // 1. Apply stress to all elements based on the current load step
        for (int i = 0; i < MESH_ELEMENTS; ++i) {
            if (mesh[i].damage < 1.0) {
                // Stress = Young's Modulus * Strain
                mesh[i].stress += mesh[i].youngs_modulus * STRAIN_PER_STEP;
            }
        }

        int next_tips_count = 0;

        // 2. Iterate current tips, check for fracture, and identify potential new tips
        for (int i = 0; i < num_crack_tips; ++i) {
            int tip_idx = crack_tip_indices[i];
            
            double stress_intensity = mesh[tip_idx].stress * STRESS_INTENSITY_FACTOR_GEOMETRY;
            
            if (stress_intensity > mesh[tip_idx].fracture_toughness) {
                // Element fractures
                mesh[tip_idx].damage = 1.0;
                final_damage_metric += tip_idx; // Accumulate result

                // Add neighbors to the list of next potential tips
                if (tip_idx > 0 && mesh[tip_idx - 1].damage < 1.0) {
                    next_crack_tip_indices[next_tips_count++] = tip_idx - 1;
                }
                if (tip_idx < MESH_ELEMENTS - 1 && mesh[tip_idx + 1].damage < 1.0) {
                    next_crack_tip_indices[next_tips_count++] = tip_idx + 1;
                }
            } else {
                // Element does not fracture, it remains a tip for the next step
                next_crack_tip_indices[next_tips_count++] = tip_idx;
            }
        }

        // 3. Reset the `is_crack_tip` flag for all current tips
        for (int i = 0; i < num_crack_tips; ++i) {
            mesh[crack_tip_indices[i]].is_crack_tip = 0;
        }

        // 4. Compact the list of next tips, de-duplicating on-the-fly
        num_crack_tips = 0;
        for (int i = 0; i < next_tips_count; ++i) {
            int new_tip_idx = next_crack_tip_indices[i];
            if (!mesh[new_tip_idx].is_crack_tip) {
                mesh[new_tip_idx].is_crack_tip = 1;
                crack_tip_indices[num_crack_tips++] = new_tip_idx;
            }
        }
    }
}

void cleanup() {
    free(mesh);
    free(crack_tip_indices);
    free(next_crack_tip_indices);
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
    printf("%lld\n", final_damage_metric);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
