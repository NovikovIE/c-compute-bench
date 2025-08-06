#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator - Do Not Modify ---
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
// --- End of MT19937 ---

// --- Benchmark Data Structures ---
typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    Vec3 position;
    Vec3 old_position;
    Vec3 acceleration;
    int is_movable;
} Particle;

// --- Global Variables ---
int grid_width;
int grid_height;
int num_particles;
int num_simulation_steps;
int constraint_iterations;

Particle *particles;
float rest_length;

double final_result = 0.0;

// --- Vector Math Helpers ---
static inline Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vec3 vec3_scale(Vec3 v, float s) {
    return (Vec3){v.x * s, v.y * s, v.z * s};
}

static inline float vec3_length(Vec3 v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <grid_width> <grid_height> <num_steps> <constraint_iter> <seed>\n", argv[0]);
        exit(1);
    }

    grid_width = atoi(argv[1]);
    grid_height = atoi(argv[2]);
    num_simulation_steps = atoi(argv[3]);
    constraint_iterations = atoi(argv[4]);
    uint32_t seed = atoi(argv[5]);
    mt_seed(seed);

    num_particles = grid_width * grid_height;
    particles = (Particle *)malloc(num_particles * sizeof(Particle));
    if (!particles) {
        fprintf(stderr, "Failed to allocate memory for particles.\n");
        exit(1);
    }

    float spacing = 0.5f;
    rest_length = spacing;
    Vec3 gravity = {0.0f, -9.81f, 0.0f};

    for (int y = 0; y < grid_height; ++y) {
        for (int x = 0; x < grid_width; ++x) {
            int i = y * grid_width + x;
            particles[i].position = (Vec3){x * spacing, 0.0f, y * spacing};
            particles[i].old_position = particles[i].position;
            particles[i].acceleration = gravity;
            if (y == 0) { // Pin the top row
                particles[i].is_movable = 0;
            } else {
                particles[i].is_movable = 1;
            }
        }
    }
}

static inline void satisfy_constraint(Particle* p1, Particle* p2, float rl) {
    Vec3 delta = vec3_sub(p2->position, p1->position);
    float current_dist = vec3_length(delta);
    // Avoid division by zero
    if (current_dist < 1e-6f) return;
    float diff = (current_dist - rl) / current_dist;
    Vec3 correction = vec3_scale(delta, 0.5f * diff);

    if (p1->is_movable) {
        p1->position = vec3_add(p1->position, correction);
    }
    if (p2->is_movable) {
        p2->position = vec3_sub(p2->position, correction);
    }
}

void run_computation() {
    float dt = 0.016f;
    float dt_sq = dt * dt;
    float rl_structural = rest_length;
    float rl_shear = rest_length * sqrtf(2.0f);

    for (int step = 0; step < num_simulation_steps; ++step) {
        // Verlet integration
        for (int i = 0; i < num_particles; ++i) {
            if (particles[i].is_movable) {
                Vec3 temp = particles[i].position;
                Vec3 displacement = vec3_sub(particles[i].position, particles[i].old_position);
                Vec3 new_pos = vec3_add(particles[i].position, displacement);
                new_pos = vec3_add(new_pos, vec3_scale(particles[i].acceleration, dt_sq));
                particles[i].position = new_pos;
                particles[i].old_position = temp;
            }
        }

        // Satisfy constraints
        for (int iter = 0; iter < constraint_iterations; ++iter) {
            // Structural constraints (horizontal and vertical)
            for (int y = 0; y < grid_height; ++y) {
                for (int x = 0; x < grid_width; ++x) {
                    Particle* p1 = &particles[y * grid_width + x];
                    if (x < grid_width - 1) {
                        satisfy_constraint(p1, &particles[y * grid_width + x + 1], rl_structural);
                    }
                    if (y < grid_height - 1) {
                        satisfy_constraint(p1, &particles[(y + 1) * grid_width + x], rl_structural);
                    }
                }
            }
            // Shear constraints (diagonal)
            for (int y = 0; y < grid_height - 1; ++y) {
                for (int x = 0; x < grid_width - 1; ++x) {
                    Particle* p1 = &particles[y * grid_width + x];
                    satisfy_constraint(p1, &particles[(y + 1) * grid_width + x + 1], rl_shear);
                    Particle* p2 = &particles[y * grid_width + x + 1];
                    satisfy_constraint(p2, &particles[(y + 1) * grid_width + x], rl_shear);
                }
            }
        }
    }

    double sum_y = 0.0;
    for (int i = 0; i < num_particles; ++i) {
        sum_y += particles[i].position.y;
    }
    final_result = sum_y;
}

void cleanup() {
    free(particles);
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
