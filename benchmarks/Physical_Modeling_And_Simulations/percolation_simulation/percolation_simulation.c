#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

// --- Mersenne Twister (Do Not Modify) ---
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

// --- Benchmark Data ---
struct Benchmark {
    int grid_size_x;
    int grid_size_y;
    int num_trials;
    int total_sites;

    // Data for computation
    int* site_indices; // For shuffling
    char* grid;        // The grid itself (0=blocked, 1=open)
    int* parent;       // DSU parent array
    int* sz;           // DSU size array

    // Result
    double final_result; // Use double for final average
};

static struct Benchmark benchmark;

// --- DSU Helper Functions ---
// Find root with path compression
static int dsu_find(int i) {
    if (benchmark.parent[i] == i) {
        return i;
    }
    return benchmark.parent[i] = dsu_find(benchmark.parent[i]);
}

// Unite by size
static void dsu_unite(int i, int j) {
    int root_i = dsu_find(i);
    int root_j = dsu_find(j);
    if (root_i != root_j) {
        if (benchmark.sz[root_i] < benchmark.sz[root_j]) {
            int temp = root_i;
            root_i = root_j;
            root_j = temp;
        }
        benchmark.parent[root_j] = root_i;
        benchmark.sz[root_i] += benchmark.sz[root_j];
    }
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s grid_size_x grid_size_y num_trials seed\n", argv[0]);
        exit(1);
    }
    benchmark.grid_size_x = atoi(argv[1]);
    benchmark.grid_size_y = atoi(argv[2]);
    benchmark.num_trials = atoi(argv[3]);
    uint32_t seed = (uint32_t)strtoul(argv[4], NULL, 10);
    
    mt_seed(seed);

    benchmark.total_sites = benchmark.grid_size_x * benchmark.grid_size_y;
    if (benchmark.total_sites <= 0) {
        fprintf(stderr, "FATAL: Invalid grid dimensions.\n");
        exit(1);
    }
    
    // +2 for virtual top and bottom nodes in DSU
    int dsu_size = benchmark.total_sites + 2;

    benchmark.site_indices = (int*)malloc(benchmark.total_sites * sizeof(int));
    benchmark.grid = (char*)malloc(benchmark.total_sites * sizeof(char));
    benchmark.parent = (int*)malloc(dsu_size * sizeof(int));
    benchmark.sz = (int*)malloc(dsu_size * sizeof(int));

    if (!benchmark.site_indices || !benchmark.grid || !benchmark.parent || !benchmark.sz) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }
}

void run_computation() {
    long long total_opened_sites = 0;
    const int gx = benchmark.grid_size_x;
    const int gy = benchmark.grid_size_y;
    const int N = benchmark.total_sites;
    const int v_top = N;
    const int v_bottom = N + 1;

    for (int t = 0; t < benchmark.num_trials; ++t) {
        // --- Reset for new trial ---
        // 1. Reset grid and DSU
        for (int i = 0; i < N; ++i) {
            benchmark.grid[i] = 0; // 0 for blocked
        }
        for (int i = 0; i < N + 2; ++i) {
            benchmark.parent[i] = i;
            benchmark.sz[i] = 1;
        }

        // 2. Create shuffled site order
        for (int i = 0; i < N; ++i) {
            benchmark.site_indices[i] = i;
        }
        // Fisher-Yates shuffle
        for (int i = N - 1; i > 0; --i) {
            int j = mt_rand() % (i + 1);
            int temp = benchmark.site_indices[i];
            benchmark.site_indices[i] = benchmark.site_indices[j];
            benchmark.site_indices[j] = temp;
        }

        // --- Run one trial ---
        for (int k = 0; k < N; ++k) {
            int idx = benchmark.site_indices[k];
            benchmark.grid[idx] = 1; // Open the site

            int r = idx / gx;
            int c = idx % gx;

            // Connect to neighbors
            // Up
            if (r > 0 && benchmark.grid[idx - gx]) dsu_unite(idx, idx - gx);
            // Down
            if (r < gy - 1 && benchmark.grid[idx + gx]) dsu_unite(idx, idx + gx);
            // Left
            if (c > 0 && benchmark.grid[idx - 1]) dsu_unite(idx, idx - 1);
            // Right
            if (c < gx - 1 && benchmark.grid[idx + 1]) dsu_unite(idx, idx + 1);

            // Connect to virtual nodes
            if (r == 0) dsu_unite(idx, v_top);
            if (r == gy - 1) dsu_unite(idx, v_bottom);

            // Check for percolation
            if (dsu_find(v_top) == dsu_find(v_bottom)) {
                total_opened_sites += (k + 1);
                break;
            }
        }
    }
    
    benchmark.final_result = (double)total_opened_sites / benchmark.num_trials;
}

void cleanup() {
    free(benchmark.site_indices);
    free(benchmark.grid);
    free(benchmark.parent);
    free(benchmark.sz);
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
    printf("%.4f\n", benchmark.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
