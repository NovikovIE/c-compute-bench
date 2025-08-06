#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// M_PI is not part of the C standard, so define it if not available (e.g., on MSVC)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- START: Mersenne Twister (MT19937) implementation ---
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
// --- END: Mersenne Twister (MT19937) implementation ---

// --- BENCHMARK DATA AND GLOBALS ---
typedef struct {
    int num_links;
    int max_iterations;

    // Robot Arm state
    double* link_lengths;     // size: num_links
    double* joint_angles;     // size: num_links

    // Target position
    double target_x;
    double target_y;

    // Pre-allocated workspace arrays
    double* cumulative_angles; // size: num_links
    double* jacobian_flat;     // size: 2 * num_links (stores a 2xN matrix)

    // Result accumulator to prevent dead code elimination
    int total_successful_solves;

} BenchmarkData;

static BenchmarkData g_data;

// --- BENCHMARK FUNCTIONS ---

// Helper to generate a random double in [min, max)
static double rand_double(double min, double max) {
    return min + ((double)mt_rand() / (double)UINT32_MAX) * (max - min);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_links> <max_iterations> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_links = atoi(argv[1]);
    g_data.max_iterations = atoi(argv[2]);
    uint32_t seed = (uint32_t)atol(argv[3]);

    if (g_data.num_links <= 0 || g_data.max_iterations <= 0) {
        fprintf(stderr, "FATAL: num_links and max_iterations must be positive integers.\n");
        exit(1);
    }

    mt_seed(seed);

    // Allocate memory
    g_data.link_lengths = (double*)malloc(g_data.num_links * sizeof(double));
    g_data.joint_angles = (double*)malloc(g_data.num_links * sizeof(double));
    g_data.cumulative_angles = (double*)malloc(g_data.num_links * sizeof(double));
    g_data.jacobian_flat = (double*)malloc(2 * g_data.num_links * sizeof(double));

    if (!g_data.link_lengths || !g_data.joint_angles || !g_data.cumulative_angles || !g_data.jacobian_flat) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize robot arm with random properties
    double max_reach = 0.0;
    for (int i = 0; i < g_data.num_links; ++i) {
        g_data.link_lengths[i] = rand_double(0.5, 1.5);
        g_data.joint_angles[i] = rand_double(-M_PI, M_PI);
        max_reach += g_data.link_lengths[i];
    }

    // Set a random, reachable target
    double target_dist = rand_double(0.1, 0.8) * max_reach;
    double target_angle = rand_double(-M_PI, M_PI);
    g_data.target_x = target_dist * cos(target_angle);
    g_data.target_y = target_dist * sin(target_angle);

    g_data.total_successful_solves = 0;
}

void run_computation() {
    const double alpha = 0.1;           // Learning rate
    const double tolerance = 1e-4;      // Convergence tolerance
    const double tolerance_sq = tolerance * tolerance;
    const int N = g_data.num_links;

    // Main solver loop
    for (int iter = 0; iter < g_data.max_iterations; ++iter) {
        // 1. Forward Kinematics: Calculate current end-effector position
        // First, compute cumulative joint angles
        g_data.cumulative_angles[0] = g_data.joint_angles[0];
        for (int i = 1; i < N; ++i) {
            g_data.cumulative_angles[i] = g_data.cumulative_angles[i - 1] + g_data.joint_angles[i];
        }

        // Then, compute end-effector position
        double current_x = 0.0, current_y = 0.0;
        for (int i = 0; i < N; ++i) {
            current_x += g_data.link_lengths[i] * cos(g_data.cumulative_angles[i]);
            current_y += g_data.link_lengths[i] * sin(g_data.cumulative_angles[i]);
        }

        // 2. Check for convergence
        double dx = g_data.target_x - current_x;
        double dy = g_data.target_y - current_y;
        if ((dx * dx + dy * dy) < tolerance_sq) {
            g_data.total_successful_solves++;
            // Reset joint angles to a new random configuration to continue workload
            for (int i = 0; i < N; ++i) {
                g_data.joint_angles[i] = rand_double(-M_PI, M_PI);
            }
            continue; // Continue to the next iteration
        }

        // 3. Calculate Jacobian matrix (J)
        for (int k = 0; k < N; ++k) {
            double J_xk = 0.0, J_yk = 0.0;
            for (int i = k; i < N; ++i) {
                J_xk -= g_data.link_lengths[i] * sin(g_data.cumulative_angles[i]);
                J_yk += g_data.link_lengths[i] * cos(g_data.cumulative_angles[i]);
            }
            g_data.jacobian_flat[k] = J_xk;          // J[0][k]
            g_data.jacobian_flat[N + k] = J_yk;    // J[1][k]
        }

        // 4. Calculate J * J^T (a 2x2 matrix)
        double jjt[2][2] = {{0, 0}, {0, 0}};
        for (int k = 0; k < N; ++k) {
            jjt[0][0] += g_data.jacobian_flat[k] * g_data.jacobian_flat[k];
            jjt[0][1] += g_data.jacobian_flat[k] * g_data.jacobian_flat[N + k];
            jjt[1][1] += g_data.jacobian_flat[N + k] * g_data.jacobian_flat[N + k];
        }
        jjt[1][0] = jjt[0][1];

        // 5. Invert J*J^T using Damped Least Squares to avoid singularities
        double damping = 1e-4;
        double det = (jjt[0][0] + damping) * (jjt[1][1] + damping) - jjt[0][1] * jjt[1][0];
        if (fabs(det) < 1e-10) { det = (det > 0 ? 1e-10 : -1e-10); }
        double inv_det = 1.0 / det;

        double jjt_inv[2][2];
        jjt_inv[0][0] = (jjt[1][1] + damping) * inv_det;
        jjt_inv[0][1] = -jjt[0][1] * inv_det;
        jjt_inv[1][0] = -jjt[1][0] * inv_det;
        jjt_inv[1][1] = (jjt[0][0] + damping) * inv_det;

        // 6. Calculate delta_p * (J*J^T)^-1
        double temp_vec_x = jjt_inv[0][0] * dx + jjt_inv[0][1] * dy;
        double temp_vec_y = jjt_inv[1][0] * dx + jjt_inv[1][1] * dy;

        // 7. Calculate delta_theta = J^T * temp_vec and update angles
        for (int k = 0; k < N; ++k) {
            double delta_theta_k = g_data.jacobian_flat[k] * temp_vec_x + g_data.jacobian_flat[N + k] * temp_vec_y;
            g_data.joint_angles[k] += alpha * delta_theta_k;
        }
    }
}

void cleanup() {
    free(g_data.link_lengths);
    free(g_data.joint_angles);
    free(g_data.cumulative_angles);
    free(g_data.jacobian_flat);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print the accumulated result to stdout to prevent dead code elimination
    printf("%d\n", g_data.total_successful_solves);

    // Print the time taken to stderr
    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
