#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

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
// --- End of Mersenne Twister ---

// --- Benchmark Data Structures ---
typedef struct {
    double x, y, z;      // Position
    double vx, vy, vz;   // Velocity
    double mass;
} Star;

typedef struct {
    double fx, fy, fz;    // Force
} Force;

// --- Global Benchmark Parameters and Data ---
static int NUM_STARS_GALAXY1;
static int NUM_STARS_GALAXY2;
static int NUM_TIME_STEPS;
static int TOTAL_STARS;

static Star* stars = NULL;
static Force* forces = NULL;

static double final_result = 0.0;

// --- Helper function for random number generation ---
double rand_double() {
    // Returns a random double in [0, 1)
    return (double)mt_rand() / 4294967296.0;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_stars_galaxy1 num_stars_galaxy2 num_time_steps seed\n", argv[0]);
        exit(1);
    }

    NUM_STARS_GALAXY1 = atoi(argv[1]);
    NUM_STARS_GALAXY2 = atoi(argv[2]);
    NUM_TIME_STEPS = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    TOTAL_STARS = NUM_STARS_GALAXY1 + NUM_STARS_GALAXY2;

    stars = (Star*)malloc(TOTAL_STARS * sizeof(Star));
    forces = (Force*)malloc(TOTAL_STARS * sizeof(Force));

    if (stars == NULL || forces == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    // Initialize Galaxy 1
    const double g1_radius = 5.0;
    const double g1_center_x = -15.0, g1_center_y = 0.0, g1_center_z = 0.0;
    const double g1_vel_x = 2.0, g1_vel_y = 0.0, g1_vel_z = 0.0;
    for (int i = 0; i < NUM_STARS_GALAXY1; i++) {
        double px, py, pz;
        do {
            px = (rand_double() * 2.0 - 1.0) * g1_radius;
            py = (rand_double() * 2.0 - 1.0) * g1_radius;
            pz = (rand_double() * 2.0 - 1.0) * g1_radius;
        } while (px*px + py*py + pz*pz > g1_radius*g1_radius);

        stars[i].x = g1_center_x + px;
        stars[i].y = g1_center_y + py;
        stars[i].z = g1_center_z + pz;
        stars[i].vx = g1_vel_x + (rand_double() - 0.5);
        stars[i].vy = g1_vel_y + (rand_double() - 0.5);
        stars[i].vz = g1_vel_z + (rand_double() - 0.5);
        stars[i].mass = 1.0 + rand_double() * 9.0; // Mass between 1.0 and 10.0
    }

    // Initialize Galaxy 2
    const double g2_radius = 3.0;
    const double g2_center_x = 15.0, g2_center_y = 0.0, g2_center_z = 0.0;
    const double g2_vel_x = -2.0, g2_vel_y = 0.0, g2_vel_z = 0.0;
    for (int i = NUM_STARS_GALAXY1; i < TOTAL_STARS; i++) {
        double px, py, pz;
        do {
            px = (rand_double() * 2.0 - 1.0) * g2_radius;
            py = (rand_double() * 2.0 - 1.0) * g2_radius;
            pz = (rand_double() * 2.0 - 1.0) * g2_radius;
        } while (px*px + py*py + pz*pz > g2_radius*g2_radius);

        stars[i].x = g2_center_x + px;
        stars[i].y = g2_center_y + py;
        stars[i].z = g2_center_z + pz;
        stars[i].vx = g2_vel_x + (rand_double() - 0.5);
        stars[i].vy = g2_vel_y + (rand_double() - 0.5);
        stars[i].vz = g2_vel_z + (rand_double() - 0.5);
        stars[i].mass = 0.5 + rand_double() * 4.5; // Mass between 0.5 and 5.0
    }
}

void run_computation() {
    const double G = 1.0;            // Scaled gravitational constant
    const double dt = 0.01;          // Time step
    const double softening_sq = 0.01; // Softening factor squared (0.1*0.1)

    for (int t = 0; t < NUM_TIME_STEPS; t++) {
        memset(forces, 0, TOTAL_STARS * sizeof(Force));

        for (int i = 0; i < TOTAL_STARS; i++) {
            for (int j = 0; j < TOTAL_STARS; j++) {
                if (i == j) continue;

                double dx = stars[j].x - stars[i].x;
                double dy = stars[j].y - stars[i].y;
                double dz = stars[j].z - stars[i].z;

                double dist_sq = dx*dx + dy*dy + dz*dz + softening_sq;
                double inv_dist = 1.0 / sqrt(dist_sq);
                double inv_dist_cubed = inv_dist * inv_dist * inv_dist;

                double force_magnitude = G * stars[i].mass * stars[j].mass * inv_dist_cubed;

                forces[i].fx += force_magnitude * dx;
                forces[i].fy += force_magnitude * dy;
                forces[i].fz += force_magnitude * dz;
            }
        }

        for (int i = 0; i < TOTAL_STARS; i++) {
            double ax = forces[i].fx / stars[i].mass;
            double ay = forces[i].fy / stars[i].mass;
            double az = forces[i].fz / stars[i].mass;

            stars[i].vx += ax * dt;
            stars[i].vy += ay * dt;
            stars[i].vz += az * dt;

            stars[i].x += stars[i].vx * dt;
            stars[i].y += stars[i].vy * dt;
            stars[i].z += stars[i].vz * dt;
        }
    }

    double total_kinetic_energy = 0.0;
    for (int i = 0; i < TOTAL_STARS; i++) {
        double v_sq = stars[i].vx * stars[i].vx + stars[i].vy * stars[i].vy + stars[i].vz * stars[i].vz;
        total_kinetic_energy += 0.5 * stars[i].mass * v_sq;
    }
    final_result = total_kinetic_energy;
}

void cleanup() {
    if (stars) {
        free(stars);
        stars = NULL;
    }
    if (forces) {
        free(forces);
        forces = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    cleanup();

    printf("%f\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
