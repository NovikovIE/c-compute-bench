#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h> 

// --- START MERSENNE TWISTER (Provided) ---
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
// --- END MERSENNE TWISTER ---

// --- BENCHMARK DATA STRUCTURES ---

// Axis-Aligned Bounding Box
typedef struct {
    float min[3]; // x, y, z
    float max[3]; // x, y, z
    int set_id;   // 0 for set A, 1 for set B
} AABB;

// Endpoint for the Sweep and Prune algorithm
typedef struct {
    float value;      // coordinate value on the sweep axis
    int is_min;       // 1 if min, 0 if max
    AABB* parent_aabb; // pointer to the AABB this endpoint belongs to
} Endpoint;

// --- GLOBAL VARIABLES ---
int num_objects_a;
int num_objects_b;

AABB *objects_a;
AABB *objects_b;

long long collision_pairs_count = 0; // Final result

// --- HELPER FUNCTIONS ---

// Generate a random float in a given range
float rand_float(float min, float max) {
    return min + (max - min) * (mt_rand() / (float)UINT32_MAX);
}

// Comparison function for qsort to sort endpoints
int compare_endpoints(const void* a, const void* b) {
    Endpoint* ep1 = (Endpoint*)a;
    Endpoint* ep2 = (Endpoint*)b;
    if (ep1->value < ep2->value) return -1;
    if (ep1->value > ep2->value) return 1;
    return 0;
}

// Check for 3D AABB intersection
int aabb_intersect(const AABB* a, const AABB* b) {
    return (a->min[0] <= b->max[0] && a->max[0] >= b->min[0]) &&
           (a->min[1] <= b->max[1] && a->max[1] >= b->min[1]) &&
           (a->min[2] <= b->max[2] && a->max[2] >= b->min[2]);
}

// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_objects_a> <num_objects_b> <seed>\n", argv[0]);
        exit(1);
    }

    num_objects_a = atoi(argv[1]);
    num_objects_b = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    mt_seed(seed);

    objects_a = (AABB*)malloc(num_objects_a * sizeof(AABB));
    objects_b = (AABB*)malloc(num_objects_b * sizeof(AABB));

    if (!objects_a || !objects_b) {
        fprintf(stderr, "Failed to allocate memory for objects.\n");
        exit(1);
    }
    
    // Generate objects for set A
    for (int i = 0; i < num_objects_a; i++) {
        float center[3] = {rand_float(-1000.0f, 1000.0f), rand_float(-1000.0f, 1000.0f), rand_float(-1000.0f, 1000.0f)};
        float size[3] = {rand_float(1.0f, 20.0f), rand_float(1.0f, 20.0f), rand_float(1.0f, 20.0f)};
        for (int j = 0; j < 3; j++) {
            objects_a[i].min[j] = center[j] - size[j] / 2.0f;
            objects_a[i].max[j] = center[j] + size[j] / 2.0f;
        }
        objects_a[i].set_id = 0;
    }

    // Generate objects for set B
    for (int i = 0; i < num_objects_b; i++) {
        float center[3] = {rand_float(-1000.0f, 1000.0f), rand_float(-1000.0f, 1000.0f), rand_float(-1000.0f, 1000.0f)};
        float size[3] = {rand_float(1.0f, 20.0f), rand_float(1.0f, 20.0f), rand_float(1.0f, 20.0f)};
        for (int j = 0; j < 3; j++) {
            objects_b[i].min[j] = center[j] - size[j] / 2.0f;
            objects_b[i].max[j] = center[j] + size[j] / 2.0f;
        }
        objects_b[i].set_id = 1;
    }
}

void run_computation() {
    int total_objects = num_objects_a + num_objects_b;
    Endpoint* endpoints_x = (Endpoint*)malloc(2 * total_objects * sizeof(Endpoint));
    if (!endpoints_x) {
        fprintf(stderr, "Failed to allocate memory for endpoints.\n");
        exit(1);
    }
    
    // Populate endpoints list from both sets
    int endpoint_idx = 0;
    for (int i = 0; i < num_objects_a; i++) {
        endpoints_x[endpoint_idx++] = (Endpoint){objects_a[i].min[0], 1, &objects_a[i]};
        endpoints_x[endpoint_idx++] = (Endpoint){objects_a[i].max[0], 0, &objects_a[i]};
    }
    for (int i = 0; i < num_objects_b; i++) {
        endpoints_x[endpoint_idx++] = (Endpoint){objects_b[i].min[0], 1, &objects_b[i]};
        endpoints_x[endpoint_idx++] = (Endpoint){objects_b[i].max[0], 0, &objects_b[i]};
    }

    // Sort endpoints along the X-axis
    qsort(endpoints_x, 2 * total_objects, sizeof(Endpoint), compare_endpoints);

    // Active lists for objects currently overlapping on the sweep line
    AABB** active_list_a = (AABB**)malloc(num_objects_a * sizeof(AABB*));
    AABB** active_list_b = (AABB**)malloc(num_objects_b * sizeof(AABB*));
    if (!active_list_a || !active_list_b) {
        fprintf(stderr, "Failed to allocate memory for active lists.\n");
        exit(1);
    }
    int active_count_a = 0;
    int active_count_b = 0;

    collision_pairs_count = 0;
    
    // Sweep through the sorted endpoints
    for (int i = 0; i < 2 * total_objects; i++) {
        Endpoint* ep = &endpoints_x[i];
        if (ep->is_min) { // start of an interval
            if (ep->parent_aabb->set_id == 0) { // Object from set A
                // Check against all active objects in set B
                for (int j = 0; j < active_count_b; j++) {
                    if (aabb_intersect(ep->parent_aabb, active_list_b[j])) {
                        collision_pairs_count++;
                    }
                }
                active_list_a[active_count_a++] = ep->parent_aabb;
            } else { // Object from set B
                // Check against all active objects in set A
                for (int j = 0; j < active_count_a; j++) {
                    if (aabb_intersect(ep->parent_aabb, active_list_a[j])) {
                        collision_pairs_count++;
                    }
                }
                active_list_b[active_count_b++] = ep->parent_aabb;
            }
        } else { // end of an interval
            if (ep->parent_aabb->set_id == 0) { // Object from set A
                // Remove from active list A (linear scan to find, O(1) to remove)
                for (int j = 0; j < active_count_a; j++) {
                    if (active_list_a[j] == ep->parent_aabb) {
                        active_list_a[j] = active_list_a[active_count_a - 1];
                        active_count_a--;
                        break;
                    }
                }
            } else { // Object from set B
                // Remove from active list B (linear scan to find, O(1) to remove)
                for (int j = 0; j < active_count_b; j++) {
                    if (active_list_b[j] == ep->parent_aabb) {
                        active_list_b[j] = active_list_b[active_count_b - 1];
                        active_count_b--;
                        break;
                    }
                }
            }
        }
    }
    
    free(endpoints_x);
    free(active_list_a);
    free(active_list_b);
}

void cleanup() {
    free(objects_a);
    free(objects_b);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", collision_pairs_count);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
