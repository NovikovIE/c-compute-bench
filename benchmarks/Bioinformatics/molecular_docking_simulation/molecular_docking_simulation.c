#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <float.h>

// --- Mersenne Twister (verbatim) ---
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
// --- end Mersenne Twister ---

// Helper to generate a random float in [min, max]
float rand_float(float min, float max) {
    return min + ((float)mt_rand() / (float)UINT32_MAX) * (max - min);
}

// Represents an atom in 3D space with a partial charge
typedef struct {
    float x, y, z;
    float charge;
} Atom;

// Global pointers and parameters for the simulation
static Atom *ligand_atoms_ptr;
static Atom *receptor_atoms_ptr;
static int num_ligand_atoms;
static int num_receptor_atoms;
static int search_space_size;

// Global variable to store the final result
static double final_best_energy;

// Setup: parse arguments, allocate and initialize data structures
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <ligand_atoms> <receptor_atoms> <search_space_size> <seed>\n", argv[0]);
        exit(1);
    }

    num_ligand_atoms = atoi(argv[1]);
    num_receptor_atoms = atoi(argv[2]);
    search_space_size = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    if (num_ligand_atoms <= 0 || num_receptor_atoms <= 0 || search_space_size <= 0) {
        fprintf(stderr, "Error: All parameters must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    ligand_atoms_ptr = (Atom*)malloc(num_ligand_atoms * sizeof(Atom));
    receptor_atoms_ptr = (Atom*)malloc(num_receptor_atoms * sizeof(Atom));

    if (!ligand_atoms_ptr || !receptor_atoms_ptr) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }
    
    // Initialize ligand atoms in a small volume
    for (int i = 0; i < num_ligand_atoms; i++) {
        ligand_atoms_ptr[i].x = rand_float(-5.0f, 5.0f);
        ligand_atoms_ptr[i].y = rand_float(-5.0f, 5.0f);
        ligand_atoms_ptr[i].z = rand_float(-5.0f, 5.0f);
        ligand_atoms_ptr[i].charge = rand_float(-1.0f, 1.0f);
    }
    
    // Initialize receptor atoms in a larger volume
    for (int i = 0; i < num_receptor_atoms; i++) {
        receptor_atoms_ptr[i].x = rand_float(-50.0f, 50.0f);
        receptor_atoms_ptr[i].y = rand_float(-50.0f, 50.0f);
        receptor_atoms_ptr[i].z = rand_float(-50.0f, 50.0f);
        receptor_atoms_ptr[i].charge = rand_float(-1.0f, 1.0f);
    }
}

// Computation: simulate docking by searching for the minimum energy configuration
void run_computation() {
    double best_energy = DBL_MAX;
    const float search_range = 50.0f;

    for (int i = 0; i < search_space_size; ++i) {
        // Generate a random translation (pose) for the ligand
        float trans_x = rand_float(-search_range, search_range);
        float trans_y = rand_float(-search_range, search_range);
        float trans_z = rand_float(-search_range, search_range);
        
        double current_pose_energy = 0.0;
        
        // Calculate interaction energy for the current pose
        for (int j = 0; j < num_ligand_atoms; ++j) {
            float lx = ligand_atoms_ptr[j].x + trans_x;
            float ly = ligand_atoms_ptr[j].y + trans_y;
            float lz = ligand_atoms_ptr[j].z + trans_z;
            float l_charge = ligand_atoms_ptr[j].charge;

            for (int k = 0; k < num_receptor_atoms; ++k) {
                float dx = lx - receptor_atoms_ptr[k].x;
                float dy = ly - receptor_atoms_ptr[k].y;
                float dz = lz - receptor_atoms_ptr[k].z;
                
                // Use squared distance to avoid costly sqrt in the inner loop
                float dist_sq = dx * dx + dy * dy + dz * dz;
                // Add a small epsilon to avoid division by zero
                if (dist_sq < 1e-6f) {
                    dist_sq = 1e-6f;
                }
                
                // Simplified energy model (e.g., electrostatic-like interaction)
                // A real simulation would use complex potentials like Lennard-Jones.
                // This E = q1*q2/r^2 model serves as a good computational workload.
                current_pose_energy += (l_charge * receptor_atoms_ptr[k].charge) / dist_sq;
            }
        }
        
        if (current_pose_energy < best_energy) {
            best_energy = current_pose_energy;
        }
    }
    
    final_best_energy = best_energy;
}

// Cleanup: free all allocated memory
void cleanup() {
    free(ligand_atoms_ptr);
    free(receptor_atoms_ptr);
    ligand_atoms_ptr = NULL;
    receptor_atoms_ptr = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    
    setup_benchmark(argc, argv);
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print the final calculated result to stdout
    printf("%f\n", final_best_energy);
    
    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);
    
    return 0;
}
