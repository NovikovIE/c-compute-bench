#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Mersenne Twister (verbatim)
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

// Utility function for random doubles in [0, 1)
double rand_double() {
    return (double)mt_rand() / (UINT32_MAX + 1.0);
}

// --- Benchmark Globals ---

// Parameters
int NUM_PARTICLES;
int MAP_WIDTH;
int MAP_HEIGHT;
int NUM_UPDATES;
unsigned int SEED;

// Data structures
typedef struct {
    double x;
    double y;
    double theta;
    double weight;
} Particle;

Particle *particles;
char *map;
double final_result = 0.0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_particles map_width map_height num_updates seed\n", argv[0]);
        exit(1);
    }
    NUM_PARTICLES = atoi(argv[1]);
    MAP_WIDTH = atoi(argv[2]);
    MAP_HEIGHT = atoi(argv[3]);
    NUM_UPDATES = atoi(argv[4]);
    SEED = (unsigned int)strtoul(argv[5], NULL, 10);

    mt_seed(SEED);

    // Allocate map
    map = (char*)malloc((size_t)MAP_WIDTH * MAP_HEIGHT * sizeof(char));
    if (!map) {
        perror("Failed to allocate map");
        exit(1);
    }

    // Generate map with obstacles (approx 20% density)
    for (int y = 0; y < MAP_HEIGHT; y++) {
        for (int x = 0; x < MAP_WIDTH; x++) {
            if (y == 0 || y == MAP_HEIGHT - 1 || x == 0 || x == MAP_WIDTH - 1 || rand_double() < 0.2) {
                map[y * MAP_WIDTH + x] = 1; // Obstacle
            } else {
                map[y * MAP_WIDTH + x] = 0; // Free space
            }
        }
    }

    // Allocate particles
    particles = (Particle*)malloc((size_t)NUM_PARTICLES * sizeof(Particle));
    if (!particles) {
        perror("Failed to allocate particles");
        exit(1);
    }

    // Initialize particles randomly in free space
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int x, y;
        do {
            x = (int)(rand_double() * (MAP_WIDTH - 2)) + 1;
            y = (int)(rand_double() * (MAP_HEIGHT - 2)) + 1;
        } while (map[y * MAP_WIDTH + x] == 1); // Ensure not starting in an obstacle
        
        particles[i].x = (double)x;
        particles[i].y = (double)y;
        particles[i].theta = rand_double() * 2.0 * M_PI;
        particles[i].weight = 1.0 / NUM_PARTICLES;
    }
}

// Sensor model: returns sum of distances to obstacles in 4 directions
double sense(double x, double y) {
    if (x < 0 || x >= MAP_WIDTH || y < 0 || y >= MAP_HEIGHT) return 0.0;

    double distances[4] = {0};
    const int max_dist = 100;
    int ix = (int)x;
    int iy = (int)y;

    // North
    for (int i = 1; i < max_dist; i++) { int ny = iy - i; if (ny < 0 || map[ny * MAP_WIDTH + ix] == 1) { distances[0] = i; break; } }
    // East
    for (int i = 1; i < max_dist; i++) { int nx = ix + i; if (nx >= MAP_WIDTH || map[iy * MAP_WIDTH + nx] == 1) { distances[1] = i; break; } }
    // South
    for (int i = 1; i < max_dist; i++) { int ny = iy + i; if (ny >= MAP_HEIGHT || map[ny * MAP_WIDTH + ix] == 1) { distances[2] = i; break; } }
    // West
    for (int i = 1; i < max_dist; i++) { int nx = ix - i; if (nx < 0 || map[iy * MAP_WIDTH + nx] == 1) { distances[3] = i; break; } }
    
    return distances[0] + distances[1] + distances[2] + distances[3];
}

void run_computation() {
    const double move_dist = 2.0;
    const double turn_angle = 0.1;
    const double move_noise = 0.5;
    const double turn_noise = 0.05;
    const double sensor_noise_stddev_sq = 25.0; // 5.0 * 5.0

    Particle* new_particles = (Particle*)malloc((size_t)NUM_PARTICLES * sizeof(Particle));

    double true_robot_x = MAP_WIDTH / 2.0;
    double true_robot_y = MAP_HEIGHT / 2.0;

    for (int i = 0; i < NUM_UPDATES; i++) {
        // 1. Prediction (Motion Update)
        for (int p = 0; p < NUM_PARTICLES; p++) {
            double noisy_dist = move_dist + (rand_double() - 0.5) * 2.0 * move_noise;
            double noisy_turn = turn_angle + (rand_double() - 0.5) * 2.0 * turn_noise;
            
            particles[p].theta = fmod(particles[p].theta + noisy_turn, 2.0 * M_PI);
            particles[p].x += cos(particles[p].theta) * noisy_dist;
            particles[p].y += sin(particles[p].theta) * noisy_dist;

            if (particles[p].x < 0) particles[p].x = 0; else if (particles[p].x >= MAP_WIDTH) particles[p].x = MAP_WIDTH - 1;
            if (particles[p].y < 0) particles[p].y = 0; else if (particles[p].y >= MAP_HEIGHT) particles[p].y = MAP_HEIGHT - 1;
        }

        // 2. Correction (Measurement Update)
        double true_measurement = sense(true_robot_x, true_robot_y);
        double total_weight = 0.0;
        for (int p = 0; p < NUM_PARTICLES; p++) {
            double particle_measurement = sense(particles[p].x, particles[p].y);
            double diff = true_measurement - particle_measurement;
            double weight = exp(-(diff * diff) / (2.0 * sensor_noise_stddev_sq));
            particles[p].weight = weight;
            total_weight += weight;
        }
        
        // 3. Normalize Weights
        if (total_weight > 1e-12) {
            for (int p = 0; p < NUM_PARTICLES; p++) { particles[p].weight /= total_weight; }
        } else {
            for (int p = 0; p < NUM_PARTICLES; p++) { particles[p].weight = 1.0 / NUM_PARTICLES; }
        }

        // 4. Resampling (Low Variance Sampler)
        double r = rand_double() / NUM_PARTICLES;
        double c = particles[0].weight;
        int j = 0;
        for (int m = 0; m < NUM_PARTICLES; m++) {
            double u = r + (double)m / NUM_PARTICLES;
            while (u > c) {
                j++;
                c += particles[j].weight;
            }
            new_particles[m] = particles[j];
        }
        memcpy(particles, new_particles, (size_t)NUM_PARTICLES * sizeof(Particle));
        
        true_robot_x += cos(i * 0.05) * 2.0;
        true_robot_y += sin(i * 0.05) * 2.0;
    }

    free(new_particles);

    double sum_x = 0.0, sum_y = 0.0;
    for (int p = 0; p < NUM_PARTICLES; p++) {
        sum_x += particles[p].x;
        sum_y += particles[p].y;
    }
    final_result = sum_x + sum_y;
}

void cleanup() {
    free(particles);
    free(map);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", (int)final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
