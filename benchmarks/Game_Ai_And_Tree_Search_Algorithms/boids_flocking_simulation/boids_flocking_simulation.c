#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// Mersenne Twister PRNG (DO NOT MODIFY)
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

// Boid structure
typedef struct {
    float px, py; // position
    float vx, vy; // velocity
} Boid;

// --- Boid Simulation Constants ---
const float PI = 3.1415926535f;
const float VISUAL_RANGE_SQR = 5625.0f;    // 75*75
const float PROTECTED_RANGE_SQR = 64.0f;  // 8*8
const float CENTERING_FACTOR = 0.0005f;   // Cohesion strength
const float AVOID_FACTOR = 0.05f;         // Separation strength
const float MATCHING_FACTOR = 0.05f;      // Alignment strength
const float MAX_SPEED = 6.0f;
const float MIN_SPEED = 3.0f;
const float WORLD_WIDTH = 800.0f;
const float WORLD_HEIGHT = 600.0f;

// Global benchmark parameters
int NUM_AGENTS;
int NUM_SIMULATION_STEPS;

// Global data structures
Boid* boids = NULL;
float result_checksum = 0.0f;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_agents> <num_simulation_steps> <seed>\n", argv[0]);
        exit(1);
    }
    NUM_AGENTS = atoi(argv[1]);
    NUM_SIMULATION_STEPS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_AGENTS <= 0 || NUM_SIMULATION_STEPS <= 0) {
        fprintf(stderr, "Parameters must be positive integers.\n");
        exit(1);
    }
    
    mt_seed(seed);

    boids = (Boid*)malloc(NUM_AGENTS * sizeof(Boid));
    if (boids == NULL) {
        fprintf(stderr, "Failed to allocate memory for boids.\n");
        exit(1);
    }

    for (int i = 0; i < NUM_AGENTS; ++i) {
        boids[i].px = (mt_rand() / (float)UINT32_MAX) * WORLD_WIDTH;
        boids[i].py = (mt_rand() / (float)UINT32_MAX) * WORLD_HEIGHT;
        float angle = (mt_rand() / (float)UINT32_MAX) * 2.0f * PI;
        float speed = (mt_rand() / (float)UINT32_MAX) * (MAX_SPEED - MIN_SPEED) + MIN_SPEED;
        boids[i].vx = cosf(angle) * speed;
        boids[i].vy = sinf(angle) * speed;
    }
}

void run_computation() {
    for (int step = 0; step < NUM_SIMULATION_STEPS; ++step) {
        for (int i = 0; i < NUM_AGENTS; ++i) {
            float avg_vx = 0, avg_vy = 0;   // For Alignment
            float center_x = 0, center_y = 0; // For Cohesion
            float close_dx = 0, close_dy = 0; // For Separation
            int neighbor_count = 0;

            // For each boid, find its neighbors and calculate forces
            for (int j = 0; j < NUM_AGENTS; ++j) {
                if (i == j) continue;

                float dx = boids[j].px - boids[i].px;
                float dy = boids[j].py - boids[i].py;
                float dist_sqr = dx * dx + dy * dy;

                if (dist_sqr < VISUAL_RANGE_SQR) {
                    neighbor_count++;
                    avg_vx += boids[j].vx;
                    avg_vy += boids[j].vy;
                    center_x += boids[j].px;
                    center_y += boids[j].py;

                    if (dist_sqr < PROTECTED_RANGE_SQR) {
                        // Inverse weighting can be complex, a simple sum is good for a benchmark
                        close_dx -= dx;
                        close_dy -= dy;
                    }
                }
            }
            
            // Apply Separation force
            boids[i].vx += close_dx * AVOID_FACTOR;
            boids[i].vy += close_dy * AVOID_FACTOR;

            if (neighbor_count > 0) {
                // Apply Cohesion force
                center_x /= neighbor_count;
                center_y /= neighbor_count;
                boids[i].vx += (center_x - boids[i].px) * CENTERING_FACTOR;
                boids[i].vy += (center_y - boids[i].py) * CENTERING_FACTOR;
                 
                // Apply Alignment force
                avg_vx /= neighbor_count;
                avg_vy /= neighbor_count;
                boids[i].vx += (avg_vx - boids[i].vx) * MATCHING_FACTOR;
                boids[i].vy += (avg_vy - boids[i].vy) * MATCHING_FACTOR;
            }
            
            // Limit speed
            float speed = sqrtf(boids[i].vx * boids[i].vx + boids[i].vy * boids[i].vy);
            if (speed > MAX_SPEED) {
                boids[i].vx = (boids[i].vx / speed) * MAX_SPEED;
                boids[i].vy = (boids[i].vy / speed) * MAX_SPEED;
            }
            if (speed < MIN_SPEED) {
                boids[i].vx = (boids[i].vx / speed) * MIN_SPEED;
                boids[i].vy = (boids[i].vy / speed) * MIN_SPEED;
            }
            
            // Update position
            boids[i].px += boids[i].vx;
            boids[i].py += boids[i].vy;

            // World wrapping
            if (boids[i].px < 0) boids[i].px += WORLD_WIDTH;
            if (boids[i].px > WORLD_WIDTH) boids[i].px -= WORLD_WIDTH;
            if (boids[i].py < 0) boids[i].py += WORLD_HEIGHT;
            if (boids[i].py > WORLD_HEIGHT) boids[i].py -= WORLD_HEIGHT;
        }
    }

    float checksum = 0.0f;
    for (int i = 0; i < NUM_AGENTS; ++i) {
        checksum += boids[i].px + boids[i].py + boids[i].vx + boids[i].vy;
    }
    result_checksum = checksum;
}

void cleanup() {
    free(boids);
    boids = NULL;
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
    printf("%f\n", result_checksum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
