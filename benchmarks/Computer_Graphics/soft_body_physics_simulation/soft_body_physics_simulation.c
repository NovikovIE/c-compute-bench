#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

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

// --- Benchmark Data Structures ---
typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    int nodeA, nodeB;
    float rest_length;
} Spring;

// Parameters
static int NUM_NODES;
static int NUM_SPRINGS;
static int NUM_SIMULATION_STEPS;
static int SOLVER_ITERATIONS;

// Simulation data
static Vec3* positions;
static Vec3* prev_positions;
static Spring* springs;

// Result
static float final_result = 0.0f;

// Helper function
float rand_float(float min, float max) {
    return min + (max - min) * ((float)mt_rand() / (float)UINT32_MAX);
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_nodes num_springs num_simulation_steps solver_iterations seed\n", argv[0]);
        exit(1);
    }

    NUM_NODES = atoi(argv[1]);
    NUM_SPRINGS = atoi(argv[2]);
    NUM_SIMULATION_STEPS = atoi(argv[3]);
    SOLVER_ITERATIONS = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    // Allocate memory
    positions = (Vec3*)malloc(NUM_NODES * sizeof(Vec3));
    prev_positions = (Vec3*)malloc(NUM_NODES * sizeof(Vec3));
    springs = (Spring*)malloc(NUM_SPRINGS * sizeof(Spring));

    if (!positions || !prev_positions || !springs) {
        fprintf(stderr, "Failed to allocate memory\n");
        exit(1);
    }
    
    // Initialize nodes
    for (int i = 0; i < NUM_NODES; i++) {
        positions[i].x = rand_float(-1.0f, 1.0f);
        positions[i].y = rand_float(1.0f, 3.0f);
        positions[i].z = rand_float(-1.0f, 1.0f);
        prev_positions[i] = positions[i];
    }

    // Initialize springs
    for (int i = 0; i < NUM_SPRINGS; i++) {
        int nodeA, nodeB;
        do {
            nodeA = mt_rand() % NUM_NODES;
            nodeB = mt_rand() % NUM_NODES;
        } while (nodeA == nodeB);

        springs[i].nodeA = nodeA;
        springs[i].nodeB = nodeB;

        Vec3 pA = positions[nodeA];
        Vec3 pB = positions[nodeB];
        float dx = pA.x - pB.x;
        float dy = pA.y - pB.y;
        float dz = pA.z - pB.z;
        springs[i].rest_length = sqrtf(dx*dx + dy*dy + dz*dz);
    }
}

void run_computation() {
    const float dt = 0.016f;
    const float dt_sq = dt * dt;
    const Vec3 gravity = {0.0f, -9.81f, 0.0f};
    const float damping = 0.995f;

    for (int step = 0; step < NUM_SIMULATION_STEPS; step++) {
        // Predict positions (Verlet integration)
        #pragma omp parallel for
        for (int i = 0; i < NUM_NODES; i++) {
            Vec3 velocity = {
                (positions[i].x - prev_positions[i].x) * damping,
                (positions[i].y - prev_positions[i].y) * damping,
                (positions[i].z - prev_positions[i].z) * damping
            };
            
            prev_positions[i] = positions[i];
            
            positions[i].x += velocity.x + gravity.x * dt_sq;
            positions[i].y += velocity.y + gravity.y * dt_sq;
            positions[i].z += velocity.z + gravity.z * dt_sq;
        }

        // Solve constraints
        for (int iter = 0; iter < SOLVER_ITERATIONS; iter++) {
            // Spring constraints
            #pragma omp parallel for
            for (int i = 0; i < NUM_SPRINGS; i++) {
                int idxA = springs[i].nodeA;
                int idxB = springs[i].nodeB;

                Vec3 pA = positions[idxA];
                Vec3 pB = positions[idxB];
                
                Vec3 delta = {pA.x - pB.x, pA.y - pB.y, pA.z - pB.z};
                
                float dist = sqrtf(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);
                
                if (dist > 1e-6f) {
                    float diff = (dist - springs[i].rest_length) / dist;
                    float correction_scalar = diff * 0.5f; // 0.5 because we move both points
                    
                    Vec3 correction_vec = {
                        delta.x * correction_scalar,
                        delta.y * correction_scalar,
                        delta.z * correction_scalar
                    };

                    positions[idxA].x -= correction_vec.x;
                    positions[idxA].y -= correction_vec.y;
                    positions[idxA].z -= correction_vec.z;

                    positions[idxB].x += correction_vec.x;
                    positions[idxB].y += correction_vec.y;
                    positions[idxB].z += correction_vec.z;
                }
            }
            
            // Boundary constraints (ground plane)
            #pragma omp parallel for
            for (int i = 0; i < NUM_NODES; i++) {
                if (positions[i].y < 0.0f) {
                    positions[i].y = 0.0f;
                }
            }
        }
    }

    // Accumulate a final result to prevent dead code elimination
    float total_x_pos = 0.0f;
    for (int i = 0; i < NUM_NODES; i++) {
        total_x_pos += positions[i].x;
    }
    final_result = total_x_pos;
}

void cleanup() {
    free(positions);
    free(prev_positions);
    free(springs);
    positions = NULL;
    prev_positions = NULL;
    springs = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();
    
    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
