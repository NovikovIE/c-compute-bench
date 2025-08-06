#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h> // For sqrt
#include <float.h> // For DBL_MAX

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


// --- Global Benchmark Data ---
int g_width;
int g_height;
int g_num_sources;
float* g_cost_surface = NULL;
double* g_distance_surface = NULL;
int* g_source_points_indices = NULL; // Array of flattened 2D indices
double g_final_sum = 0.0;

// -- Priority Queue (Min-Heap) for Dijkstra's Algorithm --
typedef struct {
    int index;
    double distance;
} PQNode;

PQNode* g_pq = NULL;
int g_pq_size = 0;
int g_pq_capacity = 0;

void pq_swap(PQNode* a, PQNode* b) {
    PQNode temp = *a;
    *a = *b;
    *b = temp;
}

void sift_up(int i) {
    while (i > 0) {
        int p = (i - 1) / 2;
        if (g_pq[i].distance < g_pq[p].distance) {
            pq_swap(&g_pq[i], &g_pq[p]);
            i = p;
        } else {
            break;
        }
    }
}

void sift_down(int i) {
    int min_index = i;
    while (1) {
        int l = 2 * i + 1;
        int r = 2 * i + 2;
        if (l < g_pq_size && g_pq[l].distance < g_pq[min_index].distance) {
            min_index = l;
        }
        if (r < g_pq_size && g_pq[r].distance < g_pq[min_index].distance) {
            min_index = r;
        }
        if (i != min_index) {
            pq_swap(&g_pq[i], &g_pq[min_index]);
            i = min_index;
        } else {
            break;
        }
    }
}

void pq_push(int index, double distance) {
    if (g_pq_size >= g_pq_capacity) {
        return;
    }
    g_pq[g_pq_size].index = index;
    g_pq[g_pq_size].distance = distance;
    sift_up(g_pq_size);
    g_pq_size++;
}

PQNode pq_pop() {
    PQNode top = g_pq[0];
    g_pq_size--;
    if (g_pq_size > 0) {
        g_pq[0] = g_pq[g_pq_size];
        sift_down(0);
    }
    return top;
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s cost_surface_width cost_surface_height num_source_points seed\n", argv[0]);
        exit(1);
    }

    g_width = atoi(argv[1]);
    g_height = atoi(argv[2]);
    g_num_sources = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);

    mt_seed(seed);

    long long num_cells = (long long)g_width * g_height;

    // Allocate memory
    g_cost_surface = (float*)malloc(num_cells * sizeof(float));
    g_distance_surface = (double*)malloc(num_cells * sizeof(double));
    g_source_points_indices = (int*)malloc(g_num_sources * sizeof(int));
    g_pq_capacity = num_cells;
    g_pq = (PQNode*)malloc(g_pq_capacity * sizeof(PQNode));

    if (!g_cost_surface || !g_distance_surface || !g_source_points_indices || !g_pq) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize cost surface with random values (1.0 to 100.0)
    for (long long i = 0; i < num_cells; i++) {
        g_cost_surface[i] = (mt_rand() / (float)UINT32_MAX) * 99.0f + 1.0f;
    }
    
    // Initialize distance surface
    for (long long i = 0; i < num_cells; i++) {
        g_distance_surface[i] = DBL_MAX;
    }

    // Generate unique source points
    for (int i = 0; i < g_num_sources; i++) {
        int index = mt_rand() % num_cells;
        
        int is_duplicate = 0;
        for (int j = 0; j < i; j++) {
            if (g_source_points_indices[j] == index) {
                is_duplicate = 1;
                break;
            }
        }
        if (is_duplicate) {
            i--; // Retry
            continue;
        }

        g_source_points_indices[i] = index;
    }
}

void run_computation() {
    // Dijkstra's algorithm to calculate cost distance
    
    // 1. Initialize priority queue with source points
    for (int i = 0; i < g_num_sources; i++) {
        int source_idx = g_source_points_indices[i];
        g_distance_surface[source_idx] = 0.0;
        pq_push(source_idx, 0.0);
    }

    const double SQRT2 = sqrt(2.0);

    // 2. Process nodes from the priority queue
    while (g_pq_size > 0) {
        PQNode current_node = pq_pop();
        int u_idx = current_node.index;
        double u_dist = current_node.distance;

        if (u_dist > g_distance_surface[u_idx]) {
            continue;
        }

        int u_x = u_idx % g_width;
        int u_y = u_idx / g_width;

        // Explore 8 neighbors
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0) continue;

                int v_x = u_x + dx;
                int v_y = u_y + dy;

                if (v_x >= 0 && v_x < g_width && v_y >= 0 && v_y < g_height) {
                    int v_idx = v_y * g_width + v_x;
                    
                    double move_dist = (dx == 0 || dy == 0) ? 1.0 : SQRT2;
                    double travel_cost = (g_cost_surface[u_idx] + g_cost_surface[v_idx]) / 2.0 * move_dist;
                    double new_dist = u_dist + travel_cost;

                    if (new_dist < g_distance_surface[v_idx]) {
                        g_distance_surface[v_idx] = new_dist;
                        pq_push(v_idx, new_dist);
                    }
                }
            }
        }
    }

    // 3. Post-computation: calculate an aggregate checksum
    double sum = 0.0;
    long long num_cells = (long long)g_width * g_height;
    for (long long i = 0; i < num_cells; i++) {
        if (g_distance_surface[i] < DBL_MAX) {
            sum += g_distance_surface[i];
        }
    }
    g_final_sum = sum;
}

void cleanup() {
    free(g_cost_surface);
    free(g_distance_surface);
    free(g_source_points_indices);
    free(g_pq);
    g_cost_surface = NULL;
    g_distance_surface = NULL;
    g_source_points_indices = NULL;
    g_pq = NULL;
}

int main(int argc, char *argv[]) {
    struct timespec start, end;
    double time_taken;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", g_final_sum);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
