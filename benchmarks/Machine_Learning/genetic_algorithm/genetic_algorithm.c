#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

/************************************************************************
 *                  Mersenne Twister (MT19937)                          *
 *                 (Do Not Modify - Included Verbatim)                  *
 ************************************************************************/

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

/************************************************************************
 *                  BENCHMARK IMPLEMENTATION                          *
 ************************************************************************/

// --- Benchmark Globals & Constants ---

// Encapsulates all benchmark data
typedef struct {
    int population_size;
    int num_genes;
    int num_generations;

    // The target "perfect" individual
    uint8_t *target_genes;

    // Two buffers for populations, to swap between generations
    uint8_t **population_a;
    uint8_t **population_b;

    // Pointer to the current active population
    uint8_t **current_population;
    // Pointer to the next generation's population buffer
    uint8_t **next_population;

    // Fitness scores for the current population
    int *fitness;

    // Final result: best fitness found across all generations
    int best_fitness_overall;
} BenchmarkData;

static BenchmarkData g_data;

// --- Helper Functions in support of Core Algorithm ---

// Helper to allocate a 2D population array using a contiguous memory block
static uint8_t** allocate_population_buffer(int pop_size, int gene_count) {
    uint8_t** population = (uint8_t**)malloc(pop_size * sizeof(uint8_t*));
    if (!population) return NULL;

    size_t genes_size = (size_t)pop_size * gene_count * sizeof(uint8_t);
    uint8_t* genes_block = (uint8_t*)malloc(genes_size);
    if (!genes_block) {
        free(population);
        return NULL;
    }

    for (int i = 0; i < pop_size; ++i) {
        population[i] = &genes_block[(size_t)i * gene_count];
    }

    return population;
}

// Calculate fitness for the entire population
static void calculate_fitness() {
    for (int i = 0; i < g_data.population_size; ++i) {
        int score = 0;
        for (int j = 0; j < g_data.num_genes; ++j) {
            if (g_data.current_population[i][j] == g_data.target_genes[j]) {
                score++;
            }
        }
        g_data.fitness[i] = score;
    }
}

// Select a parent using tournament selection
static int select_parent() {
    const int TOURNAMENT_SIZE = 5;
    int best_participant_index = -1;
    int best_fitness = -1;

    for (int i = 0; i < TOURNAMENT_SIZE; ++i) {
        int participant_index = mt_rand() % g_data.population_size;
        if (g_data.fitness[participant_index] > best_fitness) {
            best_fitness = g_data.fitness[participant_index];
            best_participant_index = participant_index;
        }
    }
    return best_participant_index;
}

// --- Core Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <population_size> <num_genes> <num_generations> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.population_size = atoi(argv[1]);
    g_data.num_genes = atoi(argv[2]);
    g_data.num_generations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    g_data.target_genes = (uint8_t*)malloc(g_data.num_genes * sizeof(uint8_t));
    g_data.fitness = (int*)malloc(g_data.population_size * sizeof(int));
    g_data.population_a = allocate_population_buffer(g_data.population_size, g_data.num_genes);
    g_data.population_b = allocate_population_buffer(g_data.population_size, g_data.num_genes);

    if (!g_data.target_genes || !g_data.fitness || !g_data.population_a || !g_data.population_b) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    g_data.current_population = g_data.population_a;
    g_data.next_population = g_data.population_b;
    g_data.best_fitness_overall = 0;

    // Initialize target solution with random bits
    for (int i = 0; i < g_data.num_genes; i++) {
        g_data.target_genes[i] = mt_rand() % 2;
    }

    // Initialize the starting population with random individuals
    for (int i = 0; i < g_data.population_size; i++) {
        for (int j = 0; j < g_data.num_genes; j++) {
            g_data.current_population[i][j] = mt_rand() % 2;
        }
    }
}

void run_computation() {
    const float MUTATION_RATE = 0.01f;
    const uint32_t MUTATION_THRESHOLD = (uint32_t)(MUTATION_RATE * (float)UINT32_MAX);

    for (int gen = 0; gen < g_data.num_generations; ++gen) {
        calculate_fitness();

        int best_current_fitness = -1;
        int elite_index = -1;
        for (int i = 0; i < g_data.population_size; ++i) {
            if (g_data.fitness[i] > best_current_fitness) {
                best_current_fitness = g_data.fitness[i];
                elite_index = i;
            }
        }

        if (best_current_fitness > g_data.best_fitness_overall) {
            g_data.best_fitness_overall = best_current_fitness;
        }

        // Elitism: copy the best individual from current gen to the next
        memcpy(g_data.next_population[0], g_data.current_population[elite_index], g_data.num_genes * sizeof(uint8_t));

        // Create the rest of the new generation
        for (int i = 1; i < g_data.population_size; ++i) {
            int p1_idx = select_parent();
            int p2_idx = select_parent();
            uint8_t* parent1 = g_data.current_population[p1_idx];
            uint8_t* parent2 = g_data.current_population[p2_idx];
            uint8_t* offspring = g_data.next_population[i];

            // Crossover
            int crossover_point = mt_rand() % (g_data.num_genes);
            memcpy(offspring, parent1, crossover_point * sizeof(uint8_t));
            memcpy(offspring + crossover_point, parent2 + crossover_point, (g_data.num_genes - crossover_point) * sizeof(uint8_t));

            // Mutation
            for (int j = 0; j < g_data.num_genes; ++j) {
                if (mt_rand() < MUTATION_THRESHOLD) {
                    offspring[j] = 1 - offspring[j]; // Flip the bit
                }
            }
        }

        // Swap population buffers for the next iteration
        uint8_t **temp_pop = g_data.current_population;
        g_data.current_population = g_data.next_population;
        g_data.next_population = temp_pop;
    }

    // Final fitness check on the last generation
    calculate_fitness();
    for (int i = 0; i < g_data.population_size; ++i) {
        if (g_data.fitness[i] > g_data.best_fitness_overall) {
            g_data.best_fitness_overall = g_data.fitness[i];
        }
    }
}

void cleanup() {
    free(g_data.population_a[0]); // Free the contiguous block of genes
    free(g_data.population_a);    // Free the array of pointers

    free(g_data.population_b[0]);
    free(g_data.population_b);

    free(g_data.target_genes);
    free(g_data.fitness);
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
    printf("%d\n", g_data.best_fitness_overall);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
