#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (DO NOT MODIFY) ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Globals ---
typedef struct { float x, y, z; } Vec3;
typedef struct { Vec3 pos; Vec3 vel; } Particle;
typedef struct { Vec3 pos; float strength; float radius_sq; } Influencer;

int num_particles;
int num_simulation_steps;
int num_influencers;

Particle* particles;
Influencer* influencers;
double final_result = 0.0;

// Helper to generate a random float between 0 and 1
float rand_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <num_particles> <num_simulation_steps> <num_influencers> <seed>\n", argv[0]);
        exit(1);
    }

    num_particles = atoi(argv[1]);
    num_simulation_steps = atoi(argv[2]);
    num_influencers = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    particles = (Particle*)malloc(num_particles * sizeof(Particle));
    if (!particles) {
        fprintf(stderr, "Failed to allocate memory for particles.\n");
        exit(1);
    }

    influencers = (Influencer*)malloc(num_influencers * sizeof(Influencer));
    if (!influencers) {
        fprintf(stderr, "Failed to allocate memory for influencers.\n");
        free(particles);
        exit(1);
    }

    // Initialize particles in a plane at the bottom
    for (int i = 0; i < num_particles; ++i) {
        particles[i].pos.x = rand_float() * 100.0f - 50.0f;
        particles[i].pos.y = rand_float() * 10.0f; // Start in a thin layer
        particles[i].pos.z = rand_float() * 100.0f - 50.0f;
        particles[i].vel.x = 0.0f;
        particles[i].vel.y = 0.0f;
        particles[i].vel.z = 0.0f;
    }

    // Initialize influencers randomly in the volume
    for (int i = 0; i < num_influencers; ++i) {
        influencers[i].pos.x = rand_float() * 100.0f - 50.0f;
        influencers[i].pos.y = rand_float() * 80.0f + 20.0f; // Positioned above the particles
        influencers[i].pos.z = rand_float() * 100.0f - 50.0f;
        influencers[i].strength = rand_float() * 5.0f + 1.0f;
        float radius = rand_float() * 15.0f + 10.0f;
        influencers[i].radius_sq = radius * radius;
    }
}

void run_computation() {
    const float dt = 0.1f;
    const float damping = 0.99f;
    const float base_buoyancy = 0.05f;

    for (int step = 0; step < num_simulation_steps; ++step) {
        for (int i = 0; i < num_particles; ++i) {
            Particle* p = &particles[i];

            Vec3 total_force = {0.0f, base_buoyancy, 0.0f};

            for (int j = 0; j < num_influencers; ++j) {
                Influencer* inf = &influencers[j];

                float dx = p->pos.x - inf->pos.x;
                float dy = p->pos.y - inf->pos.y;
                float dz = p->pos.z - inf->pos.z;
                float dist_sq = dx * dx + dy * dy + dz * dz;

                if (dist_sq < inf->radius_sq && dist_sq > 1e-6) {
                    float dist = sqrtf(dist_sq);
                    float inv_dist = 1.0f / dist;
                    
                    // Repulsive force from influencer center to push particles away
                    float repulsion = inf->strength * (1.0f - dist / sqrtf(inf->radius_sq));
                    total_force.x += dx * inv_dist * repulsion;
                    total_force.y += dy * inv_dist * repulsion;
                    total_force.z += dz * inv_dist * repulsion;

                    // Add a swirl/vortex component
                    float swirl = inf->strength * 0.5f;
                    total_force.x += -dz * inv_dist * swirl;
                    total_force.z +=  dx * inv_dist * swirl;
                }
            }

            // Update velocity
            p->vel.x += total_force.x * dt;
            p->vel.y += total_force.y * dt;
            p->vel.z += total_force.z * dt;

            // Apply damping
            p->vel.x *= damping;
            p->vel.y *= damping;
            p->vel.z *= damping;

            // Update position
            p->pos.x += p->vel.x * dt;
            p->pos.y += p->vel.y * dt;
            p->pos.z += p->vel.z * dt;
        }
    }

    // Calculate a final result to prevent dead code elimination
    double checksum = 0.0;
    for (int i = 0; i < num_particles; ++i) {
        checksum += particles[i].pos.x + particles[i].pos.y + particles[i].pos.z;
    }
    final_result = checksum;
}

void cleanup() {
    free(particles);
    free(influencers);
}


int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
