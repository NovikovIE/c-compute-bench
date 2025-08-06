#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

// --- MERSENNE TWISTER (MT19937) --- (DO NOT MODIFY)
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

// Benchmark-specific data structures
typedef struct {
    double x, y, z;
} Vec3D;

typedef struct {
    Vec3D pos;     // Position
    Vec3D vel;     // Velocity
    int type;      // 0: hydrophobic, 1: hydrophilic
} AminoAcid;

// Global state struct to hold all data
static struct {
    int num_amino_acids;
    int num_simulation_steps;
    AminoAcid *chain;
    Vec3D *forces;
    double result_checksum;
} g_data;

// Generates a random double between -1.0 and 1.0
double rand_double() {
    return ((double)mt_rand() / (double)UINT32_MAX) * 2.0 - 1.0;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_amino_acids> <num_simulation_steps> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_amino_acids = atoi(argv[1]);
    g_data.num_simulation_steps = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    g_data.chain = (AminoAcid*)malloc(g_data.num_amino_acids * sizeof(AminoAcid));
    g_data.forces = (Vec3D*)malloc(g_data.num_amino_acids * sizeof(Vec3D));
    if (!g_data.chain || !g_data.forces) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize amino acids in a chain-like manner
    for (int i = 0; i < g_data.num_amino_acids; ++i) {
        // Stagger initial positions to form a loose initial chain
        g_data.chain[i].pos.x = (double)i * 0.8;
        g_data.chain[i].pos.y = rand_double();
        g_data.chain[i].pos.z = rand_double();

        g_data.chain[i].vel.x = 0.0;
        g_data.chain[i].vel.y = 0.0;
        g_data.chain[i].vel.z = 0.0;

        g_data.chain[i].type = mt_rand() % 2;
    }
    g_data.result_checksum = 0.0;
}

void run_computation() {
    const double DT = 1e-5;               // Simulation time step
    const double DAMPING = 0.98;          // Velocity damping factor
    const double BOND_K = 20.0;           // Bond spring constant
    const double BOND_TARGET_DIST_SQ = 1.0; // Target bond length squared
    const double DIST_EPSILON = 1e-6;     // To prevent division by zero

    for (int step = 0; step < g_data.num_simulation_steps; ++step) {
        memset(g_data.forces, 0, g_data.num_amino_acids * sizeof(Vec3D));

        // 1. Calculate non-bonded forces (all-pairs interaction)
        for (int i = 0; i < g_data.num_amino_acids; ++i) {
            for (int j = i + 1; j < g_data.num_amino_acids; ++j) {
                double dx = g_data.chain[j].pos.x - g_data.chain[i].pos.x;
                double dy = g_data.chain[j].pos.y - g_data.chain[i].pos.y;
                double dz = g_data.chain[j].pos.z - g_data.chain[i].pos.z;
                double dist_sq = dx * dx + dy * dy + dz * dz + DIST_EPSILON;

                double attractive_coeff, repulsive_coeff;
                if (g_data.chain[i].type == 0 && g_data.chain[j].type == 0) { // Hydrophobic-Hydrophobic
                    attractive_coeff = 2.0; repulsive_coeff = 1.0;
                } else if (g_data.chain[i].type == 1 && g_data.chain[j].type == 1) { // Hydrophilic-Hydrophilic
                    attractive_coeff = 0.2; repulsive_coeff = 1.5;
                } else { // Mixed
                    attractive_coeff = 0.5; repulsive_coeff = 0.5;
                }
                
                double inv_dist_sq = 1.0 / dist_sq;
                double force_scalar = inv_dist_sq * (attractive_coeff - repulsive_coeff * inv_dist_sq);

                g_data.forces[i].x += dx * force_scalar;
                g_data.forces[i].y += dy * force_scalar;
                g_data.forces[i].z += dz * force_scalar;
                g_data.forces[j].x -= dx * force_scalar;
                g_data.forces[j].y -= dy * force_scalar;
                g_data.forces[j].z -= dz * force_scalar;
            }
        }

        // 2. Calculate bonded forces (maintain chain structure)
        for (int i = 0; i < g_data.num_amino_acids - 1; ++i) {
            double dx = g_data.chain[i+1].pos.x - g_data.chain[i].pos.x;
            double dy = g_data.chain[i+1].pos.y - g_data.chain[i].pos.y;
            double dz = g_data.chain[i+1].pos.z - g_data.chain[i].pos.z;
            double dist_sq = dx * dx + dy * dy + dz * dz + DIST_EPSILON;
            
            double error = dist_sq - BOND_TARGET_DIST_SQ;
            double force_scalar = BOND_K * error;

            g_data.forces[i].x += dx * force_scalar;
            g_data.forces[i].y += dy * force_scalar;
            g_data.forces[i].z += dz * force_scalar;
            g_data.forces[i+1].x -= dx * force_scalar;
            g_data.forces[i+1].y -= dy * force_scalar;
            g_data.forces[i+1].z -= dz * force_scalar;
        }

        // 3. Update positions and velocities using simple Euler integration
        for (int i = 0; i < g_data.num_amino_acids; ++i) {
            g_data.chain[i].vel.x = (g_data.chain[i].vel.x + g_data.forces[i].x * DT) * DAMPING;
            g_data.chain[i].vel.y = (g_data.chain[i].vel.y + g_data.forces[i].y * DT) * DAMPING;
            g_data.chain[i].vel.z = (g_data.chain[i].vel.z + g_data.forces[i].z * DT) * DAMPING;

            g_data.chain[i].pos.x += g_data.chain[i].vel.x * DT;
            g_data.chain[i].pos.y += g_data.chain[i].vel.y * DT;
            g_data.chain[i].pos.z += g_data.chain[i].vel.z * DT;
        }
    }
    
    // Calculate a final checksum to prevent dead code elimination
    double checksum = 0.0;
    for (int i = 0; i < g_data.num_amino_acids; ++i) {
        checksum += g_data.chain[i].pos.x + g_data.chain[i].pos.y + g_data.chain[i].pos.z;
    }
    g_data.result_checksum = checksum;
}

void cleanup() {
    free(g_data.chain);
    free(g_data.forces);
}

// --- MAIN FUNCTION ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", g_data.result_checksum);

    cleanup();

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
