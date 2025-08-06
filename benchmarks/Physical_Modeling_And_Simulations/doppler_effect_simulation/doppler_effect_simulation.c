#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator --- DO NOT MODIFY ---
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

// --- Benchmark Globals ---
typedef struct {
    double *pos_x, *pos_y;
    double *vel_x, *vel_y;
    double *base_freq;
} SoundSources;

SoundSources sources;
int NUM_SOURCES;
int NUM_TIME_STEPS;
double final_result;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_sources> <num_time_steps> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_SOURCES = atoi(argv[1]);
    NUM_TIME_STEPS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    sources.pos_x = (double*)malloc(NUM_SOURCES * sizeof(double));
    sources.pos_y = (double*)malloc(NUM_SOURCES * sizeof(double));
    sources.vel_x = (double*)malloc(NUM_SOURCES * sizeof(double));
    sources.vel_y = (double*)malloc(NUM_SOURCES * sizeof(double));
    sources.base_freq = (double*)malloc(NUM_SOURCES * sizeof(double));

    if (!sources.pos_x || !sources.pos_y || !sources.vel_x || !sources.vel_y || !sources.base_freq) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < NUM_SOURCES; ++i) {
        // Initial positions in a 2000x2000 meter box centered at origin
        sources.pos_x[i] = (mt_rand() / (double)UINT32_MAX) * 2000.0 - 1000.0;
        sources.pos_y[i] = (mt_rand() / (double)UINT32_MAX) * 2000.0 - 1000.0;

        // Initial velocities up to 50 m/s in any direction
        sources.vel_x[i] = (mt_rand() / (double)UINT32_MAX) * 100.0 - 50.0;
        sources.vel_y[i] = (mt_rand() / (double)UINT32_MAX) * 100.0 - 50.0;

        // Base frequencies from 200 Hz to 1200 Hz
        sources.base_freq[i] = (mt_rand() / (double)UINT32_MAX) * 1000.0 + 200.0;
    }

    final_result = 0.0;
}

void run_computation() {
    const double SPEED_OF_SOUND = 343.0; // m/s
    const double TIME_DELTA = 0.01;     // s, simulation time step
    double total_perceived_freq = 0.0;

    for (int t = 0; t < NUM_TIME_STEPS; ++t) {
        for (int i = 0; i < NUM_SOURCES; ++i) {
            // 1. Update source position
            sources.pos_x[i] += sources.vel_x[i] * TIME_DELTA;
            sources.pos_y[i] += sources.vel_y[i] * TIME_DELTA;

            // 2. Calculate Doppler effect for a stationary observer at (0,0)
            double dist_sq = sources.pos_x[i] * sources.pos_x[i] + sources.pos_y[i] * sources.pos_y[i];

            // Skip if source is at the observer's location to avoid division by zero
            if (dist_sq < 1e-9) continue;
            double dist = sqrt(dist_sq);

            // 3. Calculate radial velocity (component of velocity along the line of sight)
            // This is the dot product of the velocity vector and the unit position vector.
            // A positive value means the source is moving away from the observer.
            double v_radial_away = (sources.vel_x[i] * sources.pos_x[i] + sources.vel_y[i] * sources.pos_y[i]) / dist;

            // 4. Calculate perceived frequency using the Doppler effect formula for a moving source:
            // f_observed = f_source * (v_sound / (v_sound + v_radial_away))
            if (fabs(SPEED_OF_SOUND + v_radial_away) < 1e-9) continue; // Avoid division by zero
            double perceived_freq = sources.base_freq[i] * (SPEED_OF_SOUND / (SPEED_OF_SOUND + v_radial_away));
            
            // 5. Accumulate a value to prevent dead-code elimination
            total_perceived_freq += perceived_freq;
        }
    }
    final_result = total_perceived_freq;
}

void cleanup() {
    free(sources.pos_x);
    free(sources.pos_y);
    free(sources.vel_x);
    free(sources.vel_y);
    free(sources.base_freq);
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the final accumulated result to stdout
    printf("%f\n", final_result);

    // Print the timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
