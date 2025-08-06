#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (MT19937) ---
// Provided verbatim as per requirements
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
// --- end MT19937 ---

// --- Benchmark Data Structures ---

typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    Vec3 position;
    Vec3 velocity;
    float radius;
    float mass;
} RigidBody;

// Global struct to hold all benchmark data
static struct {
    int num_objects;
    int num_simulation_steps;
    int solver_iterations;

    RigidBody *objects;
    Vec3 gravity;
    float domain_size;
    
    // Accumulator for the final result to prevent dead-code elimination
    double result_accumulator;
} g_state;

// --- Helper Functions ---

// Generates a random float between min and max
float rand_float(float min, float max) {
    return min + (float)mt_rand() / ((float)UINT32_MAX / (max - min));
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_objects> <num_simulation_steps> <solver_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    g_state.num_objects = atoi(argv[1]);
    g_state.num_simulation_steps = atoi(argv[2]);
    g_state.solver_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    g_state.objects = (RigidBody *)malloc(g_state.num_objects * sizeof(RigidBody));
    if (!g_state.objects) {
        fprintf(stderr, "Failed to allocate memory for objects\n");
        exit(1);
    }
    
    g_state.gravity = (Vec3){0.0f, -9.81f, 0.0f};
    g_state.domain_size = 100.0f;
    const float domain_half_size = g_state.domain_size / 2.0f;

    for (int i = 0; i < g_state.num_objects; ++i) {
        g_state.objects[i].radius = rand_float(0.5f, 1.5f);
        g_state.objects[i].mass = g_state.objects[i].radius * g_state.objects[i].radius;

        float pos_range = domain_half_size - g_state.objects[i].radius;
        g_state.objects[i].position = (Vec3){
            rand_float(-pos_range, pos_range),
            rand_float(-pos_range, pos_range),
            rand_float(-pos_range, pos_range)
        };
        
        g_state.objects[i].velocity = (Vec3){
            rand_float(-2.0f, 2.0f),
            rand_float(-2.0f, 2.0f),
            rand_float(-2.0f, 2.0f)
        };
    }
    
    g_state.result_accumulator = 0.0;
}

void run_computation() {
    const float dt = 0.016f; // Timestep for a ~60Hz simulation
    const float domain_half_size = g_state.domain_size / 2.0f;
    const float restitution = 0.6f; // Coefficient of restitution (bounciness)

    for (int step = 0; step < g_state.num_simulation_steps; ++step) {
        // 1. Update velocities based on gravity
        for (int i = 0; i < g_state.num_objects; ++i) {
            g_state.objects[i].velocity.x += g_state.gravity.x * dt;
            g_state.objects[i].velocity.y += g_state.gravity.y * dt;
            g_state.objects[i].velocity.z += g_state.gravity.z * dt;
        }

        // 2. Iteratively solve constraints (collisions with domain boundaries)
        for (int k = 0; k < g_state.solver_iterations; ++k) {
            for (int i = 0; i < g_state.num_objects; ++i) {
                RigidBody *obj = &g_state.objects[i];

                if (obj->position.x - obj->radius < -domain_half_size) {
                    obj->position.x = -domain_half_size + obj->radius;
                    obj->velocity.x *= -restitution;
                } else if (obj->position.x + obj->radius > domain_half_size) {
                    obj->position.x = domain_half_size - obj->radius;
                    obj->velocity.x *= -restitution;
                }

                if (obj->position.y - obj->radius < -domain_half_size) {
                    obj->position.y = -domain_half_size + obj->radius;
                    obj->velocity.y *= -restitution;
                } else if (obj->position.y + obj->radius > domain_half_size) {
                    obj->position.y = domain_half_size - obj->radius;
                    obj->velocity.y *= -restitution;
                }

                if (obj->position.z - obj->radius < -domain_half_size) {
                    obj->position.z = -domain_half_size + obj->radius;
                    obj->velocity.z *= -restitution;
                } else if (obj->position.z + obj->radius > domain_half_size) {
                    obj->position.z = domain_half_size - obj->radius;
                    obj->velocity.z *= -restitution;
                }
            }
        }
        
        // 3. Update positions based on new velocities
        for (int i = 0; i < g_state.num_objects; ++i) {
            g_state.objects[i].position.x += g_state.objects[i].velocity.x * dt;
            g_state.objects[i].position.y += g_state.objects[i].velocity.y * dt;
            g_state.objects[i].position.z += g_state.objects[i].velocity.z * dt;
        }
    }

    // Accumulate a final result to ensure computation is not optimized away
    double sum_x = 0.0;
    for (int i = 0; i < g_state.num_objects; ++i) {
        sum_x += g_state.objects[i].position.x;
    }
    g_state.result_accumulator = sum_x;
}

void cleanup() {
    free(g_state.objects);
    g_state.objects = NULL;
}


// --- Main ---

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);
    
    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    // Print the final result to stdout
    printf("%f\n", g_state.result_accumulator);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
