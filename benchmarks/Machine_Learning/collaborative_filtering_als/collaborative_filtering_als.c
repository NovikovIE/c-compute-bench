#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// --- BENCHMARK DATA STRUCTURES ---

typedef struct {
    int user_id;
    int item_id;
    float rating;
} Rating;

int num_users;
int num_items;
long num_ratings;
int num_factors;
int num_iterations;

Rating *ratings;
float **user_factors;
float **item_factors;

double final_result; // To prevent dead code elimination

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 7) {
        fprintf(stderr, "Usage: %s num_users num_items num_ratings num_factors num_iterations seed\n", argv[0]);
        exit(1);
    }

    num_users = atoi(argv[1]);
    num_items = atoi(argv[2]);
    num_ratings = atol(argv[3]);
    num_factors = atoi(argv[4]);
    num_iterations = atoi(argv[5]);
    uint32_t seed = atoi(argv[6]);

    mt_seed(seed);

    // Allocate ratings data
    ratings = (Rating *)malloc(num_ratings * sizeof(Rating));
    if (!ratings) {
        perror("Failed to allocate ratings");
        exit(1);
    }

    // Generate random ratings
    for (long i = 0; i < num_ratings; ++i) {
        ratings[i].user_id = mt_rand() % num_users;
        ratings[i].item_id = mt_rand() % num_items;
        ratings[i].rating = (float)(mt_rand() % 501) / 100.0f; // Ratings from 0.0 to 5.0
    }

    // Allocate and initialize user factors
    user_factors = (float **)malloc(num_users * sizeof(float *));
    for (int i = 0; i < num_users; ++i) {
        user_factors[i] = (float *)malloc(num_factors * sizeof(float));
        for (int j = 0; j < num_factors; ++j) {
            user_factors[i][j] = (float)mt_rand() / (float)UINT32_MAX * 0.1f;
        }
    }

    // Allocate and initialize item factors
    item_factors = (float **)malloc(num_items * sizeof(float *));
    for (int i = 0; i < num_items; ++i) {
        item_factors[i] = (float *)malloc(num_factors * sizeof(float));
        for (int j = 0; j < num_factors; ++j) {
            item_factors[i][j] = (float)mt_rand() / (float)UINT32_MAX * 0.1f;
        }
    }
}

void run_computation() {
    const float learning_rate = 0.005f;
    const float lambda = 0.02f;

    for (int iter = 0; iter < num_iterations; ++iter) {
        // Iterate through all ratings and update factors
        for (long i = 0; i < num_ratings; ++i) {
            int u = ratings[i].user_id;
            int item = ratings[i].item_id;
            float r = ratings[i].rating;

            // Predict rating
            float prediction = 0.0f;
            for (int f = 0; f < num_factors; ++f) {
                prediction += user_factors[u][f] * item_factors[item][f];
            }

            float error = r - prediction;

            // Update user and item factors
            for (int f = 0; f < num_factors; ++f) {
                float uf = user_factors[u][f];
                float itf = item_factors[item][f];
                user_factors[u][f] += learning_rate * (error * itf - lambda * uf);
                item_factors[item][f] += learning_rate * (error * uf - lambda * itf);
            }
        }
    }

    // Calculate a final result to prevent optimization
    double sum = 0.0;
    for (int i = 0; i < num_users; ++i) {
        for (int j = 0; j < num_factors; ++j) {
            sum += user_factors[i][j];
        }
    }
    for (int i = 0; i < num_items; ++i) {
        for (int j = 0; j < num_factors; ++j) {
            sum += item_factors[i][j];
        }
    }
    final_result = sum;
}

void cleanup() {
    for (int i = 0; i < num_users; ++i) {
        free(user_factors[i]);
    }
    free(user_factors);

    for (int i = 0; i < num_items; ++i) {
        free(item_factors[i]);
    }
    free(item_factors);

    free(ratings);
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
    printf("%f\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
