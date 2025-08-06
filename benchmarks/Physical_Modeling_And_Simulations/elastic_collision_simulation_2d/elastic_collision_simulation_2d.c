#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (DO NOT MODIFY) ---
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
// --- End Mersenne Twister ---


// --- Benchmark Globals ---
typedef struct {
    double x, y;
} Vector2D;

typedef struct {
    Vector2D pos;
    Vector2D vel;
    double radius;
    double mass;
} Ball;

int NUM_BALLS;
int NUM_TIME_STEPS;
Ball* balls;
double final_result;

const double BOX_WIDTH = 1000.0;
const double BOX_HEIGHT = 1000.0;
const double MAX_INITIAL_VEL = 2.0;
const double PI = 3.14159265358979323846;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_balls> <num_time_steps> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_BALLS = atoi(argv[1]);
    NUM_TIME_STEPS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if(NUM_BALLS <= 0 || NUM_TIME_STEPS <= 0) {
        fprintf(stderr, "Error: num_balls and num_time_steps must be positive integers.\n");
        exit(1);
    }
    
    mt_seed(seed);

    balls = (Ball*)malloc(NUM_BALLS * sizeof(Ball));
    if (!balls) {
        fprintf(stderr, "Error: Memory allocation failed for balls.\n");
        exit(1);
    }

    for (int i = 0; i < NUM_BALLS; ++i) {
        balls[i].radius = (mt_rand() / (double)UINT32_MAX) * 5.0 + 5.0; // Radius [5, 10]
        balls[i].mass = balls[i].radius * balls[i].radius * PI; // Mass proportional to area
        balls[i].pos.x = (mt_rand() / (double)UINT32_MAX) * (BOX_WIDTH - 2*balls[i].radius) + balls[i].radius;
        balls[i].pos.y = (mt_rand() / (double)UINT32_MAX) * (BOX_HEIGHT - 2*balls[i].radius) + balls[i].radius;
        balls[i].vel.x = (mt_rand() / (double)UINT32_MAX) * 2.0 * MAX_INITIAL_VEL - MAX_INITIAL_VEL;
        balls[i].vel.y = (mt_rand() / (double)UINT32_MAX) * 2.0 * MAX_INITIAL_VEL - MAX_INITIAL_VEL;
    }
}


void run_computation() {
    for (int t = 0; t < NUM_TIME_STEPS; ++t) {
        // Ball-Ball Collisions
        for (int i = 0; i < NUM_BALLS; ++i) {
            for (int j = i + 1; j < NUM_BALLS; ++j) {
                double dx = balls[j].pos.x - balls[i].pos.x;
                double dy = balls[j].pos.y - balls[i].pos.y;
                double dist_sq = dx * dx + dy * dy;
                double r_sum = balls[i].radius + balls[j].radius;

                if (dist_sq < r_sum * r_sum && dist_sq > 1e-9) {
                    double dist = sqrt(dist_sq);

                    // Normal vector
                    double nx = dx / dist;
                    double ny = dy / dist;

                    // Relative velocity vector
                    double dvx = balls[j].vel.x - balls[i].vel.x;
                    double dvy = balls[j].vel.y - balls[i].vel.y;

                    // Velocity along normal
                    double vel_normal = dvx * nx + dvy * ny;

                    // Only resolve if balls are approaching
                    if (vel_normal < 0) {
                        double m1 = balls[i].mass;
                        double m2 = balls[j].mass;
                        double impulse_magnitude = (2.0 * vel_normal) / (m1 + m2);

                        balls[i].vel.x += impulse_magnitude * m2 * nx;
                        balls[i].vel.y += impulse_magnitude * m2 * ny;
                        balls[j].vel.x -= impulse_magnitude * m1 * nx;
                        balls[j].vel.y -= impulse_magnitude * m1 * ny;
                        
                        // Static resolution to prevent balls sinking into each other
                        double overlap = 0.5 * (r_sum - dist);
                        balls[i].pos.x -= overlap * nx;
                        balls[i].pos.y -= overlap * ny;
                        balls[j].pos.x += overlap * nx;
                        balls[j].pos.y += overlap * ny;
                    }
                }
            }
        }

        // Update positions and handle boundary collisions
        for (int i = 0; i < NUM_BALLS; ++i) {
            balls[i].pos.x += balls[i].vel.x;
            balls[i].pos.y += balls[i].vel.y;

            if (balls[i].pos.x - balls[i].radius < 0) {
                balls[i].vel.x *= -1;
                balls[i].pos.x = balls[i].radius;
            } else if (balls[i].pos.x + balls[i].radius > BOX_WIDTH) {
                balls[i].vel.x *= -1;
                balls[i].pos.x = BOX_WIDTH - balls[i].radius;
            }

            if (balls[i].pos.y - balls[i].radius < 0) {
                balls[i].vel.y *= -1;
                balls[i].pos.y = balls[i].radius;
            } else if (balls[i].pos.y + balls[i].radius > BOX_HEIGHT) {
                balls[i].vel.y *= -1;
                balls[i].pos.y = BOX_HEIGHT - balls[i].radius;
            }
        }
    }

    // Calculate final result to prevent dead code elimination
    double total_x_pos = 0.0;
    for (int i = 0; i < NUM_BALLS; ++i) {
        total_x_pos += balls[i].pos.x;
    }
    final_result = total_x_pos;
}


void cleanup() {
    free(balls);
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
