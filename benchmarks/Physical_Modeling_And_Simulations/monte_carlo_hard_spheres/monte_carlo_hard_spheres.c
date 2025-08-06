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

// --- Benchmark Configuration ---
typedef struct {
    double x, y, z;
} Particle;

// Parameters
static int num_particles;
static long long num_mc_cycles;
static double box_size;
static double sphere_radius;
static double max_displacement;

// Data
static Particle *particles;
static double sphere_diameter_sq;

// Result
static long long total_accepted_moves = 0;

// Helper to generate a random double between 0 and 1
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_particles num_mc_cycles box_size sphere_radius max_displacement seed\n", argv[0]);
        exit(1);
    }
    num_particles = atoi(argv[1]);
    num_mc_cycles = atoll(argv[2]);
    box_size = atof(argv[3]);
    sphere_radius = atof(argv[4]);
    max_displacement = atof(argv[5]);
    uint32_t seed = (uint32_t)atoi(argv[6]);

    mt_seed(seed);

    particles = (Particle*)malloc(num_particles * sizeof(Particle));
    if (particles == NULL) {
        fprintf(stderr, "Failed to allocate memory for particles.\n");
        exit(1);
    }
    
    sphere_diameter_sq = (2.0 * sphere_radius) * (2.0 * sphere_radius);

    // Initial particle placement on a regular grid to avoid initial overlaps
    int particles_per_dim = (int)ceil(cbrt((double)num_particles));
    double spacing = box_size / particles_per_dim;

    if (spacing < 2.0 * sphere_radius) {
        fprintf(stderr, "Particles do not fit in the box with the given parameters.\n");
        free(particles);
        exit(1);
    }

    int p_idx = 0;
    for (int i = 0; i < particles_per_dim && p_idx < num_particles; ++i) {
        for (int j = 0; j < particles_per_dim && p_idx < num_particles; ++j) {
            for (int k = 0; k < particles_per_dim && p_idx < num_particles; ++k) {
                particles[p_idx].x = i * spacing;
                particles[p_idx].y = j * spacing;
                particles[p_idx].z = k * spacing;
                p_idx++;
            }
        }
    }
}

void run_computation() {
    long long accepted_count = 0;
    
    for (long long cycle = 0; cycle < num_mc_cycles; ++cycle) {
        // One Monte Carlo sweep attempts to move each particle once on average
        for (int i = 0; i < num_particles; ++i) {
            int p_idx = mt_rand() % num_particles;
            
            // Store old position
            Particle old_pos = particles[p_idx];

            // Propose a new random position
            double dx = (rand_double() * 2.0 - 1.0) * max_displacement;
            double dy = (rand_double() * 2.0 - 1.0) * max_displacement;
            double dz = (rand_double() * 2.0 - 1.0) * max_displacement;
            
            Particle new_pos;
            new_pos.x = old_pos.x + dx;
            new_pos.y = old_pos.y + dy;
            new_pos.z = old_pos.z + dz;

            // Apply periodic boundary conditions
            new_pos.x -= floor(new_pos.x / box_size) * box_size;
            new_pos.y -= floor(new_pos.y / box_size) * box_size;
            new_pos.z -= floor(new_pos.z / box_size) * box_size;
            
            // Check for overlaps with other particles
            int has_overlap = 0;
            for (int j = 0; j < num_particles; ++j) {
                if (p_idx == j) continue;
                
                // Minimum image convention for periodic boundary conditions
                double dist_x = new_pos.x - particles[j].x;
                double dist_y = new_pos.y - particles[j].y;
                double dist_z = new_pos.z - particles[j].z;

                dist_x -= box_size * round(dist_x / box_size);
                dist_y -= box_size * round(dist_y / box_size);
                dist_z -= box_size * round(dist_z / box_size);
                
                double dist_sq = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
                
                if (dist_sq < sphere_diameter_sq) {
                    has_overlap = 1;
                    break;
                }
            }
            
            // If no overlap, accept the move
            if (!has_overlap) {
                particles[p_idx] = new_pos;
                accepted_count++;
            }
        }
    }
    total_accepted_moves = accepted_count;
}

void cleanup() {
    free(particles);
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
    
    // Print result to stdout
    printf("%lld\n", total_accepted_moves);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
