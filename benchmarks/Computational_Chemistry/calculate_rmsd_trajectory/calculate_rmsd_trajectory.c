#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- START: Mersenne Twister (MT19937) --- Do Not Modify ---
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

// Data structure for 3D coordinates
typedef struct {
    double x, y, z;
} Point3D;

// --- Global Benchmark data ---
int NUM_FRAMES;
int NUM_ATOMS_TO_ALIGN;

// A reference structure and a trajectory of other structures
Point3D *reference_structure;
Point3D **trajectory;

// Variable to store the final result, preventing dead code elimination
double total_rmsd_accumulator = 0.0;

// Generates a random coordinate between -10.0 and 10.0
double rand_coord() {
    return ((double)mt_rand() / (double)UINT32_MAX) * 20.0 - 10.0;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_frames> <num_atoms_to_align> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_FRAMES = atoi(argv[1]);
    NUM_ATOMS_TO_ALIGN = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    // Allocate memory for the reference structure
    reference_structure = (Point3D*)malloc(NUM_ATOMS_TO_ALIGN * sizeof(Point3D));
    if (!reference_structure) {
        fprintf(stderr, "Failed to allocate memory for reference_structure\n");
        exit(1);
    }

    // Allocate memory for the trajectory (array of pointers to frames)
    trajectory = (Point3D**)malloc(NUM_FRAMES * sizeof(Point3D*));
    if (!trajectory) {
        fprintf(stderr, "Failed to allocate memory for trajectory\n");
        exit(1);
    }

    // Populate reference structure with random coordinates
    for (int i = 0; i < NUM_ATOMS_TO_ALIGN; ++i) {
        reference_structure[i] = (Point3D){rand_coord(), rand_coord(), rand_coord()};
    }

    // Allocate and populate each frame in the trajectory
    for (int i = 0; i < NUM_FRAMES; ++i) {
        trajectory[i] = (Point3D*)malloc(NUM_ATOMS_TO_ALIGN * sizeof(Point3D));
        if (!trajectory[i]) {
            fprintf(stderr, "Failed to allocate memory for trajectory frame %d\n", i);
            exit(1);
        }
        for (int j = 0; j < NUM_ATOMS_TO_ALIGN; ++j) {
            trajectory[i][j] = (Point3D){rand_coord(), rand_coord(), rand_coord()};
        }
    }
}

void run_computation() {
    double accumulator = 0.0;
    for (int i = 0; i < NUM_FRAMES; ++i) {
        double sum_sq_dist = 0.0;
        for (int j = 0; j < NUM_ATOMS_TO_ALIGN; ++j) {
            double dx = reference_structure[j].x - trajectory[i][j].x;
            double dy = reference_structure[j].y - trajectory[i][j].y;
            double dz = reference_structure[j].z - trajectory[i][j].z;
            sum_sq_dist += dx*dx + dy*dy + dz*dz;
        }
        // Calculate RMSD for the current frame and add to the accumulator
        // The actual RMSD calculation is sqrt(sum_sq_dist / N),
        // but we accumulate the sum of RMSDs for all frames.
        double frame_rmsd = sqrt(sum_sq_dist / NUM_ATOMS_TO_ALIGN);
        accumulator += frame_rmsd;
    }
    total_rmsd_accumulator = accumulator;
}

void cleanup() {
    free(reference_structure);
    for (int i = 0; i < NUM_FRAMES; ++i) {
        free(trajectory[i]);
    }
    free(trajectory);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", total_rmsd_accumulator);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
