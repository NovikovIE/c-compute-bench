#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// --- Benchmark Data Structures and Globals ---
typedef struct {
    double x, y, theta;
} Pose;

typedef struct {
    int from;
    int to;
    Pose measurement;
} Edge;

int NUM_POSES;
int NUM_LOOP_CLOSURES;
int TOTAL_EDGES;

Pose* poses; // Array of pose estimates
Edge* edges; // Array of constraints (odometry and loop closures)
double final_result = 0.0;

// --- Helper Functions ---
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

double rand_noise(double scale) {
    return (rand_double() * 2.0 - 1.0) * scale;
}

// T_ac = T_ab * T_bc
static inline Pose pose_compose(Pose a, Pose b) {
    Pose result;
    double cos_a = cos(a.theta);
    double sin_a = sin(a.theta);
    result.x = a.x + cos_a * b.x - sin_a * b.y;
    result.y = a.y + sin_a * b.x + cos_a * b.y;
    result.theta = a.theta + b.theta;
    return result;
}

// T_ba = inv(T_ab)
static inline Pose pose_invert(Pose a) {
    Pose result;
    double cos_a = cos(a.theta);
    double sin_a = sin(a.theta);
    result.x = -a.x * cos_a - a.y * sin_a;
    result.y = -a.y * cos_a + a.x * sin_a;
    result.theta = -a.theta;
    return result;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_poses num_loop_closures seed\n", argv[0]);
        exit(1);
    }

    NUM_POSES = atoi(argv[1]);
    NUM_LOOP_CLOSURES = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);
    mt_seed(seed);

    TOTAL_EDGES = (NUM_POSES - 1) + NUM_LOOP_CLOSURES;

    poses = (Pose*)malloc(NUM_POSES * sizeof(Pose));
    edges = (Edge*)malloc(TOTAL_EDGES * sizeof(Edge));
    Pose* ground_truth_poses = (Pose*)malloc(NUM_POSES * sizeof(Pose));

    if (!poses || !edges || !ground_truth_poses) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // 1. Generate a ground truth trajectory (e.g., a curving path)
    ground_truth_poses[0] = (Pose){0.0, 0.0, 0.0};
    for (int i = 1; i < NUM_POSES; ++i) {
        Pose step = {0.1, 0.0, 0.02};
        ground_truth_poses[i] = pose_compose(ground_truth_poses[i - 1], step);
    }

    // 2. Generate odometry edges with noise
    double odom_pos_noise = 0.02;
    double odom_rot_noise = 0.005;
    for (int i = 0; i < NUM_POSES - 1; ++i) {
        Pose p1_inv = pose_invert(ground_truth_poses[i]);
        Pose relative_pose = pose_compose(p1_inv, ground_truth_poses[i+1]);
        Pose noise = {rand_noise(odom_pos_noise), rand_noise(odom_pos_noise), rand_noise(odom_rot_noise)};
        edges[i].from = i;
        edges[i].to = i + 1;
        edges[i].measurement = pose_compose(relative_pose, noise);
    }

    // 3. Generate loop closure edges with noise
    double loop_pos_noise = 0.1;
    double loop_rot_noise = 0.02;
    for (int i = 0; i < NUM_LOOP_CLOSURES; ++i) {
        int from_id, to_id;
        do {
            from_id = mt_rand() % (NUM_POSES - 20);
            to_id = from_id + (mt_rand() % 19) + 1; // Ensure to > from and not consecutive
        } while (from_id >= to_id);

        Pose p_from_inv = pose_invert(ground_truth_poses[from_id]);
        Pose relative_pose = pose_compose(p_from_inv, ground_truth_poses[to_id]);
        Pose noise = {rand_noise(loop_pos_noise), rand_noise(loop_rot_noise), rand_noise(loop_rot_noise)};
        
        edges[NUM_POSES - 1 + i].from = from_id;
        edges[NUM_POSES - 1 + i].to = to_id;
        edges[NUM_POSES - 1 + i].measurement = pose_compose(relative_pose, noise);
    }

    // 4. Initialize pose estimates with noise (except for the first pose, which is fixed)
    poses[0] = ground_truth_poses[0]; // Anchor the graph
    for (int i = 1; i < NUM_POSES; ++i) {
        Pose noise = {rand_noise(0.5), rand_noise(0.5), rand_noise(0.1)};
        poses[i] = pose_compose(ground_truth_poses[i], noise);
    }

    free(ground_truth_poses);
}

void run_computation() {
    const int NUM_ITERATIONS = 5;
    const double DAMPING = 0.5;
    
    for (int iter = 0; iter < NUM_ITERATIONS; ++iter) {
        for (int i = 0; i < TOTAL_EDGES; ++i) {
            Edge* e = &edges[i];
            if (e->from == 0) continue; // Skip updates originating from the anchor pose

            Pose p_from = poses[e->from];
            Pose p_to = poses[e->to];

            // Calculate predicted relative pose from current estimates
            Pose p_from_inv = pose_invert(p_from);
            Pose predicted = pose_compose(p_from_inv, p_to);

            // Calculate error: error = inv(measurement) * predicted
            Pose measurement_inv = pose_invert(e->measurement);
            Pose error = pose_compose(measurement_inv, predicted);
            error.theta = fmod(error.theta + M_PI, 2.0 * M_PI) - M_PI; // Normalize angle

            // Distribute a fraction of the error to the poses.
            // This is a simplified relaxation, not a full nonlinear solve,
            // but provides a representative computational load.
            Pose error_half = {error.x * 0.5 * DAMPING, error.y * 0.5 * DAMPING, error.theta * 0.5 * DAMPING};
            Pose error_half_inv = pose_invert(error_half);

            poses[e->from] = pose_compose(poses[e->from], error_half_inv);
            poses[e->to]   = pose_compose(pose_invert(error_half), poses[e->to]);
        }
    }

    // Calculate final graph error to prevent dead code elimination.
    double total_error_sq = 0.0;
    for (int i = 0; i < TOTAL_EDGES; ++i) {
        Edge e = edges[i];
        Pose p_from_inv = pose_invert(poses[e.from]);
        Pose predicted = pose_compose(p_from_inv, poses[e.to]);
        Pose measurement_inv = pose_invert(e.measurement);
        Pose error = pose_compose(measurement_inv, predicted);
        total_error_sq += error.x * error.x + error.y * error.y + error.theta * error.theta;
    }
    final_result = total_error_sq;
}

void cleanup() {
    free(poses);
    free(edges);
}

int main(int argc, char* argv[]) {
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
