#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

//
// --- Mersenne Twister (MT19937) Generator ---
//
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

// --- Benchmark Globals ---
int num_atoms;
int num_rotatable_bonds;
int num_mc_steps;

// Current and temporary atomic coordinates
double* coords_x;
double* coords_y;
double* coords_z;
double* temp_coords_x;
double* temp_coords_y;
double* temp_coords_z;

// Current angles of rotatable bonds
double* bond_angles;

// Final result: the lowest energy conformation found
double lowest_energy;

// --- Function Prototypes ---
double calculate_energy(const double* x, const double* y, const double* z);
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

// --- Benchmark Implementation ---

/**
 * @brief A simplified potential energy function.
 * It calculates the sum of squared distances between all pairs of atoms.
 * This O(N^2) calculation is a proxy for complex non-bonded interactions
 * like Lennard-Jones or electrostatic forces.
 */
double calculate_energy(const double* x, const double* y, const double* z) {
    double energy = 0.0;
    for (int i = 0; i < num_atoms; ++i) {
        for (int j = i + 1; j < num_atoms; ++j) {
            double dx = x[i] - x[j];
            double dy = y[i] - y[j];
            double dz = z[i] - z[j];
            energy += dx*dx + dy*dy + dz*dz;
        }
    }
    return energy;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_atoms num_rotatable_bonds num_mc_steps seed\n", argv[0]);
        exit(1);
    }

    num_atoms = atoi(argv[1]);
    num_rotatable_bonds = atoi(argv[2]);
    num_mc_steps = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    coords_x = (double*)malloc(num_atoms * sizeof(double));
    coords_y = (double*)malloc(num_atoms * sizeof(double));
    coords_z = (double*)malloc(num_atoms * sizeof(double));
    temp_coords_x = (double*)malloc(num_atoms * sizeof(double));
    temp_coords_y = (double*)malloc(num_atoms * sizeof(double));
    temp_coords_z = (double*)malloc(num_atoms * sizeof(double));
    bond_angles = (double*)malloc(num_rotatable_bonds * sizeof(double));

    if (!coords_x || !coords_y || !coords_z || !temp_coords_x || !temp_coords_y || !temp_coords_z || !bond_angles) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize atoms in a random configuration (e.g., in a 20x20x20 box)
    for (int i = 0; i < num_atoms; ++i) {
        coords_x[i] = ((double)mt_rand() / UINT32_MAX) * 20.0;
        coords_y[i] = ((double)mt_rand() / UINT32_MAX) * 20.0;
        coords_z[i] = ((double)mt_rand() / UINT32_MAX) * 20.0;
    }

    // Initialize bond angles randomly from 0 to 2*PI
    #ifndef M_PI
    #define M_PI 3.14159265358979323846
    #endif
    for (int i = 0; i < num_rotatable_bonds; ++i) {
        bond_angles[i] = ((double)mt_rand() / UINT32_MAX) * 2.0 * M_PI;
    }

    lowest_energy = 1e30; // Initialize with a very large value
}

void run_computation() {
    double current_energy = calculate_energy(coords_x, coords_y, coords_z);
    lowest_energy = current_energy;

    // A constant representing the temperature factor (kT) in the Metropolis criterion
    const double boltzmann_temp = 1.0;

    for (int step = 0; step < num_mc_steps; ++step) {
        // 1. Pick a random bond and generate a new random angle for it
        int bond_idx = mt_rand() % num_rotatable_bonds;
        double old_angle = bond_angles[bond_idx];
        double new_angle = ((double)mt_rand() / UINT32_MAX) * 2.0 * M_PI;

        // 2. Create a trial conformation by applying a transformation
        // This simplified transformation rotates all atoms around the Z-axis.
        // It serves as a computationally-intensive proxy for a real bond rotation.
        memcpy(temp_coords_x, coords_x, num_atoms * sizeof(double));
        memcpy(temp_coords_y, coords_y, num_atoms * sizeof(double));
        // Z coordinates do not change in this simplified rotation, so no need to copy `temp_coords_z`.

        double delta_angle = new_angle - old_angle;
        double cos_delta = cos(delta_angle);
        double sin_delta = sin(delta_angle);
        
        for (int i = 0; i < num_atoms; ++i) {
             temp_coords_x[i] = coords_x[i] * cos_delta - coords_y[i] * sin_delta;
             temp_coords_y[i] = coords_x[i] * sin_delta + coords_y[i] * cos_delta;
        }
        // Copy the unchanged Z coordinates to the temporary storage for the energy calculation.
        memcpy(temp_coords_z, coords_z, num_atoms * sizeof(double));

        // 3. Calculate the energy of the new conformation
        double new_energy = calculate_energy(temp_coords_x, temp_coords_y, temp_coords_z);
        
        // 4. Metropolis-Hastings acceptance criterion
        double energy_diff = new_energy - current_energy;
        if (energy_diff < 0.0 || ((double)mt_rand() / UINT32_MAX) < exp(-energy_diff / boltzmann_temp)) {
            // Acceptance: Update the system state
            current_energy = new_energy;
            bond_angles[bond_idx] = new_angle;
            memcpy(coords_x, temp_coords_x, num_atoms * sizeof(double));
            memcpy(coords_y, temp_coords_y, num_atoms * sizeof(double));

            if (current_energy < lowest_energy) {
                lowest_energy = current_energy;
            }
        }
        // Rejection: The state remains unchanged, and we proceed to the next step.
    }
}

void cleanup() {
    free(coords_x);
    free(coords_y);
    free(coords_z);
    free(temp_coords_x);
    free(temp_coords_y);
    free(temp_coords_z);
    free(bond_angles);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final result (lowest energy found) to stdout
    printf("%f\n", lowest_energy);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
