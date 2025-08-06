#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- DO NOT MODIFY --- MERSENNE TWISTER GENERATOR -------------------------
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
// --- END OF MERSENNE TWISTER -----------------------------------------------

// Simulation constants for cell states
#define EMPTY 0
#define PREY 1
#define PREDATOR 2

// Global struct to hold all benchmark data
struct BenchmarkData {
    int num_steps;
    int world_size;
    int* world;
    int* next_world;
    long long final_checksum;
};
struct BenchmarkData g_data;

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_steps> <initial_population> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_steps = atoi(argv[1]);
    int initial_population = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (g_data.num_steps <= 0 || initial_population <= 0) {
        fprintf(stderr, "ERROR: num_steps and initial_population must be positive integers.\n");
        exit(1);
    }
    
    mt_seed(seed);

    // World size is a multiple of the initial population to provide empty space
    g_data.world_size = initial_population * 2;
    int num_predators = initial_population / 20;
    if (num_predators == 0) num_predators = 1;
    
    // Ensure world_size is large enough to hold the initial population
    if (g_data.world_size < initial_population + num_predators) {
        g_data.world_size = initial_population + num_predators;
    }

    g_data.world = (int*)malloc(g_data.world_size * sizeof(int));
    g_data.next_world = (int*)malloc(g_data.world_size * sizeof(int));
    if (!g_data.world || !g_data.next_world) {
        fprintf(stderr, "ERROR: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize world to be empty
    for (int i = 0; i < g_data.world_size; ++i) {
        g_data.world[i] = EMPTY;
    }

    // Place initial prey population randomly
    int prey_placed = 0;
    while (prey_placed < initial_population) {
        int idx = mt_rand() % g_data.world_size;
        if (g_data.world[idx] == EMPTY) {
            g_data.world[idx] = PREY;
            prey_placed++;
        }
    }

    // Place initial predator population randomly
    int predators_placed = 0;
    while (predators_placed < num_predators) {
        int idx = mt_rand() % g_data.world_size;
        if (g_data.world[idx] == EMPTY) {
            g_data.world[idx] = PREDATOR;
            predators_placed++;
        }
    }
    
    g_data.final_checksum = 0;
}

void run_computation() {
    int* current_world = g_data.world;
    int* next_gen_world = g_data.next_world;

    for (int step = 0; step < g_data.num_steps; ++step) {
        for (int i = 0; i < g_data.world_size; ++i) {
            // Get neighbors with wrap-around (toroidal world)
            int left = (i == 0) ? g_data.world_size - 1 : i - 1;
            int right = (i == g_data.world_size - 1) ? 0 : i + 1;

            int num_prey_neighbors = (current_world[left] == PREY) + (current_world[right] == PREY);
            int num_predator_neighbors = (current_world[left] == PREDATOR) + (current_world[right] == PREDATOR);

            int current_cell_state = current_world[i];
            int next_cell_state = current_cell_state;

            switch (current_cell_state) {
                case EMPTY:
                    // Reproduce if there is exactly one prey neighbor and no predator neighbors
                    if (num_prey_neighbors == 1 && num_predator_neighbors == 0) {
                        next_cell_state = PREY;
                    }
                    break;
                case PREY:
                    // Eaten if any predator neighbors
                    if (num_predator_neighbors > 0) {
                        next_cell_state = EMPTY;
                    }
                    break;
                case PREDATOR:
                    // Starve if no prey neighbors
                    if (num_prey_neighbors == 0) {
                        next_cell_state = EMPTY;
                    }
                    break;
            }
            next_gen_world[i] = next_cell_state;
        }

        // Swap world buffers for the next iteration
        int* temp = current_world;
        current_world = next_gen_world;
        next_gen_world = temp;
    }
    
    // Ensure the final state is in the original 'world' buffer for consistent cleanup
    if (current_world != g_data.world) {
        for(int i = 0; i < g_data.world_size; ++i) {
            g_data.world[i] = current_world[i];
        }
    }
    
    // Calculate a checksum of the final world state to prevent dead code elimination
    long long checksum = 0;
    for (int i = 0; i < g_data.world_size; i++) {
        checksum += g_data.world[i] * (long long)(i + 1);
    }
    g_data.final_checksum = checksum;
}

void cleanup() {
    free(g_data.world);
    free(g_data.next_world);
    g_data.world = NULL;
    g_data.next_world = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", g_data.final_checksum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
