#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA STRUCTURES ---
typedef struct {
    double cash;
    int* portfolio; // size: num_assets
    double risk_aversion; // Tendency to avoid risk (0.0=none, 1.0=high)
    double optimism;      // Tendency to buy (0.0=pessimist, 1.0=optimist)
} Agent;

typedef struct {
    double price;
    double volatility;    // Price fluctuation factor
} Asset;

// --- GLOBAL PARAMETERS AND DATA ---
int num_agents;
int num_simulation_steps;
int num_assets;

Agent* agents;
Asset* assets;

double final_market_value;

// Utility function to get a random double between 0.0 and 1.0
double rand_double() {
    return mt_rand() / (double)UINT32_MAX;
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_agents num_simulation_steps num_assets seed\n", argv[0]);
        exit(1);
    }

    num_agents = atoi(argv[1]);
    num_simulation_steps = atoi(argv[2]);
    num_assets = atoi(argv[3]);
    uint32_t seed = atoi(argv[4]);

    mt_seed(seed);

    // Allocate memory
    assets = (Asset*)malloc(num_assets * sizeof(Asset));
    agents = (Agent*)malloc(num_agents * sizeof(Agent));
    if (!assets || !agents) {
        fprintf(stderr, "Memory allocation failed for assets or agents.\n");
        exit(1);
    }

    // Initialize assets
    for (int i = 0; i < num_assets; i++) {
        assets[i].price = 50.0 + rand_double() * 100.0; // Price between 50 and 150
        assets[i].volatility = 0.005 + rand_double() * 0.02; // Volatility between 0.5% and 2.5%
    }

    // Initialize agents
    for (int i = 0; i < num_agents; i++) {
        agents[i].cash = 10000.0 + rand_double() * 40000.0; // Cash between 10k and 50k
        agents[i].risk_aversion = rand_double();
        agents[i].optimism = rand_double();
        agents[i].portfolio = (int*)malloc(num_assets * sizeof(int));
        if (!agents[i].portfolio) {
             fprintf(stderr, "Memory allocation failed for agent portfolio.\n");
             exit(1);
        }
        for (int j = 0; j < num_assets; j++) {
            agents[i].portfolio[j] = 0; // Start with no assets
        }
    }
}

void run_computation() {
    for (int step = 0; step < num_simulation_steps; ++step) {
        // 1. Update asset prices based on a simple random walk model
        for (int j = 0; j < num_assets; ++j) {
            double price_change_factor = 1.0 + assets[j].volatility * (rand_double() - 0.5) * 2.0;
            assets[j].price *= price_change_factor;
            if (assets[j].price < 0.01) assets[j].price = 0.01; // Prevent price from becoming zero/negative
        }

        // 2. Agents make buy/sell decisions
        for (int i = 0; i < num_agents; ++i) {
            for (int j = 0; j < num_assets; ++j) {
                // Create a trading signal based on agent traits and a random factor
                double random_impulse = rand_double() - 0.5;
                double signal = agents[i].optimism - agents[i].risk_aversion + random_impulse;

                if (signal > 0.15) { // Threshold to buy
                    if (agents[i].cash >= assets[j].price) {
                        agents[i].portfolio[j]++;
                        agents[i].cash -= assets[j].price;
                    }
                } else if (signal < -0.15) { // Threshold to sell
                    if (agents[i].portfolio[j] > 0) {
                        agents[i].portfolio[j]--;
                        agents[i].cash += assets[j].price;
                    }
                }
            }
        }
    }

    // 3. Calculate final total market value as the result
    double total_value = 0.0;
    for (int i = 0; i < num_agents; ++i) {
        total_value += agents[i].cash;
        for (int j = 0; j < num_assets; ++j) {
            total_value += agents[i].portfolio[j] * assets[j].price;
        }
    }
    final_market_value = total_value;
}

void cleanup() {
    for (int i = 0; i < num_agents; ++i) {
        free(agents[i].portfolio);
    }
    free(agents);
    free(assets);
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
    printf("%f\n", final_market_value);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
