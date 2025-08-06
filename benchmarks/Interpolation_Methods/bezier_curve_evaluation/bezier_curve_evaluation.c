#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (MT19937) Generator - DO NOT MODIFY ---
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

// --- Benchmark-specific data structures and globals ---
typedef struct {
    double x;
    double y;
} Point;

int NUM_CONTROL_POINTS;
int NUM_EVAL_POINTS;

Point *control_points = NULL;
double *eval_timesteps = NULL;
Point *evaluated_points = NULL;
Point *temp_points = NULL; // Scratch space for De Casteljau's algorithm

double final_result = 0.0; // Accumulated result to prevent dead code elimination

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_control_points num_eval_points seed\n", argv[0]);
        exit(1);
    }

    NUM_CONTROL_POINTS = atoi(argv[1]);
    NUM_EVAL_POINTS = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    control_points = (Point *)malloc(NUM_CONTROL_POINTS * sizeof(Point));
    eval_timesteps = (double *)malloc(NUM_EVAL_POINTS * sizeof(double));
    evaluated_points = (Point *)malloc(NUM_EVAL_POINTS * sizeof(Point));
    temp_points = (Point *)malloc(NUM_CONTROL_POINTS * sizeof(Point));

    if (!control_points || !eval_timesteps || !evaluated_points || !temp_points) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < NUM_CONTROL_POINTS; ++i) {
        control_points[i].x = (double)mt_rand() / (double)UINT32_MAX * 100.0;
        control_points[i].y = (double)mt_rand() / (double)UINT32_MAX * 100.0;
    }

    for (int i = 0; i < NUM_EVAL_POINTS; ++i) {
        eval_timesteps[i] = (double)i / (double)(NUM_EVAL_POINTS - 1);
    }
}

void run_computation() {
    double total_sum = 0.0;
    for (int i = 0; i < NUM_EVAL_POINTS; ++i) {
        double t = eval_timesteps[i];

        // Initialize temp points with control points for this evaluation
        for (int j = 0; j < NUM_CONTROL_POINTS; ++j) {
            temp_points[j] = control_points[j];
        }

        // De Casteljau's algorithm
        for (int k = 1; k < NUM_CONTROL_POINTS; ++k) {
            for (int j = 0; j < NUM_CONTROL_POINTS - k; ++j) {
                temp_points[j].x = (1.0 - t) * temp_points[j].x + t * temp_points[j + 1].x;
                temp_points[j].y = (1.0 - t) * temp_points[j].y + t * temp_points[j + 1].y;
            }
        }

        evaluated_points[i] = temp_points[0];
        total_sum += evaluated_points[i].x + evaluated_points[i].y;
    }
    final_result = total_sum;
}

void cleanup() {
    free(control_points);
    free(eval_timesteps);
    free(evaluated_points);
    free(temp_points);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout to prevent optimizer from removing the work
    printf("%f\n", final_result);

    // Print timing information to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
