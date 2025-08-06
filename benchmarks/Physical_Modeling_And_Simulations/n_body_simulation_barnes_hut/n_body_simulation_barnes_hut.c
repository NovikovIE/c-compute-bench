#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

// --- Benchmark Globals ---
int NUM_BODIES;
int NUM_TIME_STEPS;
double TIME_DELTA;
double THETA_PARAMETER;

const double G_CONST = 1.0; // Normalized gravitational constant
const double EPSILON = 1e-3; // Softening parameter to avoid singularities

double final_result = 0.0;

typedef struct {
    double x, y, z;
} Vector3D;

typedef struct {
    Vector3D pos;
    Vector3D vel;
    Vector3D acc;
    double mass;
} Body;

#define MAX_CHILDREN 8
typedef struct Node {
    int is_leaf;
    int num_bodies;
    int body_index; // -1 if internal node

    Vector3D center_of_mass;
    double total_mass;

    Vector3D min_bound;
    Vector3D max_bound;
    double s_width; // width of the cubic cell

    struct Node* children[MAX_CHILDREN];
} Node;

Body* bodies;
Node* node_pool;
size_t next_node_idx;
size_t node_pool_size;

// --- Helper Functions ---
double mt_rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

void reset_node_pool() {
    next_node_idx = 0;
}

Node* alloc_node() {
    if (next_node_idx >= node_pool_size) {
        fprintf(stderr, "FATAL: Node pool exhausted.\n");
        exit(1);
    }
    Node* node = &node_pool[next_node_idx++];
    memset(node, 0, sizeof(Node));
    node->is_leaf = 1;
    node->body_index = -1;
    return node;
}

void get_bounding_box(Vector3D* min_b, Vector3D* max_b) {
    *min_b = (Vector3D){1e10, 1e10, 1e10};
    *max_b = (Vector3D){-1e10, -1e10, -1e10};
    for (int i = 0; i < NUM_BODIES; ++i) {
        if (bodies[i].pos.x < min_b->x) min_b->x = bodies[i].pos.x;
        if (bodies[i].pos.y < min_b->y) min_b->y = bodies[i].pos.y;
        if (bodies[i].pos.z < min_b->z) min_b->z = bodies[i].pos.z;
        if (bodies[i].pos.x > max_b->x) max_b->x = bodies[i].pos.x;
        if (bodies[i].pos.y > max_b->y) max_b->y = bodies[i].pos.y;
        if (bodies[i].pos.z > max_b->z) max_b->z = bodies[i].pos.z;
    }
    double dx = max_b->x - min_b->x;
    double dy = max_b->y - min_b->y;
    double dz = max_b->z - min_b->z;
    double max_dim = fmax(dx, fmax(dy, dz)) + 1e-9; // Add small epsilon
    Vector3D center = {(min_b->x + max_b->x)/2.0, (min_b->y + max_b->y)/2.0, (min_b->z + max_b->z)/2.0};
    *min_b = (Vector3D){center.x - max_dim/2.0, center.y - max_dim/2.0, center.z - max_dim/2.0};
    *max_b = (Vector3D){center.x + max_dim/2.0, center.y + max_dim/2.0, center.z + max_dim/2.0};
}

int get_octant(const Vector3D* pos, const Node* node) {
    Vector3D center = {(node->min_bound.x + node->max_bound.x) / 2.0, (node->min_bound.y + node->max_bound.y) / 2.0, (node->min_bound.z + node->max_bound.z) / 2.0};
    int octant = 0;
    if (pos->x >= center.x) octant |= 1;
    if (pos->y >= center.y) octant |= 2;
    if (pos->z >= center.z) octant |= 4;
    return octant;
}

void get_child_bounds(Vector3D* child_min, Vector3D* child_max, const Node* parent, int octant) {
    Vector3D parent_center = {(parent->min_bound.x + parent->max_bound.x)/2.0, (parent->min_bound.y + parent->max_bound.y)/2.0, (parent->min_bound.z + parent->max_bound.z)/2.0};
    child_min->x = (octant & 1) ? parent_center.x : parent->min_bound.x;
    child_max->x = (octant & 1) ? parent->max_bound.x : parent_center.x;
    child_min->y = (octant & 2) ? parent_center.y : parent->min_bound.y;
    child_max->y = (octant & 2) ? parent->max_bound.y : parent_center.y;
    child_min->z = (octant & 4) ? parent_center.z : parent->min_bound.z;
    child_max->z = (octant & 4) ? parent->max_bound.z : parent_center.z;
}

void insert_body(Node* node, int body_idx) {
    if (node->is_leaf) {
        if (node->num_bodies == 0) {
            node->body_index = body_idx;
        } else {
            node->is_leaf = 0;
            int old_body_idx = node->body_index;
            node->body_index = -1;

            int old_octant = get_octant(&bodies[old_body_idx].pos, node);
            node->children[old_octant] = alloc_node();
            get_child_bounds(&node->children[old_octant]->min_bound, &node->children[old_octant]->max_bound, node, old_octant);
            node->children[old_octant]->s_width = (node->children[old_octant]->max_bound.x - node->children[old_octant]->min_bound.x);
            insert_body(node->children[old_octant], old_body_idx);

            int new_octant = get_octant(&bodies[body_idx].pos, node);
            if (!node->children[new_octant]) {
                node->children[new_octant] = alloc_node();
                get_child_bounds(&node->children[new_octant]->min_bound, &node->children[new_octant]->max_bound, node, new_octant);
                node->children[new_octant]->s_width = (node->children[new_octant]->max_bound.x - node->children[new_octant]->min_bound.x);
            }
            insert_body(node->children[new_octant], body_idx);
        }
    } else {
        int octant = get_octant(&bodies[body_idx].pos, node);
        if (!node->children[octant]) {
            node->children[octant] = alloc_node();
            get_child_bounds(&node->children[octant]->min_bound, &node->children[octant]->max_bound, node, octant);
            node->children[octant]->s_width = (node->children[octant]->max_bound.x - node->children[octant]->min_bound.x);
        }
        insert_body(node->children[octant], body_idx);
    }
    node->num_bodies++;
}

void compute_com(Node *node) {
    if (node->is_leaf) {
        if(node->body_index != -1) {
            node->total_mass = bodies[node->body_index].mass;
            node->center_of_mass = bodies[node->body_index].pos;
        }
    } else {
        node->total_mass = 0.0;
        node->center_of_mass = (Vector3D){0,0,0};
        for (int i = 0; i < MAX_CHILDREN; ++i) {
            if (node->children[i]) {
                compute_com(node->children[i]);
                node->total_mass += node->children[i]->total_mass;
                node->center_of_mass.x += node->children[i]->center_of_mass.x * node->children[i]->total_mass;
                node->center_of_mass.y += node->children[i]->center_of_mass.y * node->children[i]->total_mass;
                node->center_of_mass.z += node->children[i]->center_of_mass.z * node->children[i]->total_mass;
            }
        }
        if (node->total_mass > 1e-9) {
            node->center_of_mass.x /= node->total_mass;
            node->center_of_mass.y /= node->total_mass;
            node->center_of_mass.z /= node->total_mass;
        }
    }
}

void compute_force_on_body(int body_idx, Node* node) {
    if (!node || node->num_bodies == 0) return;

    double dx = node->center_of_mass.x - bodies[body_idx].pos.x;
    double dy = node->center_of_mass.y - bodies[body_idx].pos.y;
    double dz = node->center_of_mass.z - bodies[body_idx].pos.z;
    double dist_sq = dx*dx + dy*dy + dz*dz;

    if (node->is_leaf) {
        if (node->body_index != -1 && node->body_index != body_idx) {
            double dist_sq_soft = dist_sq + EPSILON * EPSILON;
            double force_mag_over_dist = (G_CONST * bodies[body_idx].mass * node->total_mass) / (dist_sq_soft * sqrt(dist_sq_soft));
            bodies[body_idx].acc.x += dx * force_mag_over_dist;
            bodies[body_idx].acc.y += dy * force_mag_over_dist;
            bodies[body_idx].acc.z += dz * force_mag_over_dist;
        }
    } else {
        double s = node->s_width;
        if ((s*s / dist_sq) < (THETA_PARAMETER * THETA_PARAMETER)) {
            double dist_sq_soft = dist_sq + EPSILON * EPSILON;
            double force_mag_over_dist = (G_CONST * bodies[body_idx].mass * node->total_mass) / (dist_sq_soft * sqrt(dist_sq_soft));
            bodies[body_idx].acc.x += dx * force_mag_over_dist;
            bodies[body_idx].acc.y += dy * force_mag_over_dist;
            bodies[body_idx].acc.z += dz * force_mag_over_dist;
        } else {
            for (int i = 0; i < MAX_CHILDREN; ++i) {
                compute_force_on_body(body_idx, node->children[i]);
            }
        }
    }
}

// --- Main Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_bodies num_time_steps time_delta theta_parameter seed\n", argv[0]);
        exit(1);
    }
    NUM_BODIES = atoi(argv[1]);
    NUM_TIME_STEPS = atoi(argv[2]);
    TIME_DELTA = atof(argv[3]);
    THETA_PARAMETER = atof(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    if(NUM_BODIES <= 0 || NUM_TIME_STEPS <= 0 || TIME_DELTA <= 0 || THETA_PARAMETER < 0) {
        fprintf(stderr, "Invalid parameters.\n");
        exit(1);
    }

    mt_seed(seed);

    bodies = (Body*)malloc(NUM_BODIES * sizeof(Body));
    if (!bodies) {
        fprintf(stderr, "Failed to allocate memory for bodies.\n");
        exit(1);
    }
    
    node_pool_size = 10 * NUM_BODIES;
    node_pool = (Node*)malloc(node_pool_size * sizeof(Node));
    if (!node_pool) {
        fprintf(stderr, "Failed to allocate memory for node pool.\n");
        free(bodies);
        exit(1);
    }

    // Initialize bodies in a sphere
    for (int i = 0; i < NUM_BODIES; i++) {
        double r = mt_rand_double();
        double theta = mt_rand_double() * 2.0 * M_PI;
        double phi = acos(2.0 * mt_rand_double() - 1.0);

        bodies[i].pos = (Vector3D){r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)};
        bodies[i].vel = (Vector3D){0.0, 0.0, 0.0};
        bodies[i].acc = (Vector3D){0.0, 0.0, 0.0};
        bodies[i].mass = 1.0 + mt_rand_double() * 9.0;
    }
}

void run_computation() {
    for (int t = 0; t < NUM_TIME_STEPS; ++t) {
        for (int i = 0; i < NUM_BODIES; ++i) {
            bodies[i].acc = (Vector3D){0.0, 0.0, 0.0};
        }

        reset_node_pool();
        Vector3D min_b, max_b;
        get_bounding_box(&min_b, &max_b);
        Node* root = alloc_node();
        root->min_bound = min_b;
        root->max_bound = max_b;
        root->s_width = root->max_bound.x - root->min_bound.x;

        for (int i=0; i<NUM_BODIES; ++i) {
            insert_body(root, i);
        }

        compute_com(root);

        for (int i = 0; i < NUM_BODIES; ++i) {
            compute_force_on_body(i, root);
        }

        for(int i = 0; i < NUM_BODIES; ++i) {
            bodies[i].vel.x += bodies[i].acc.x * TIME_DELTA / bodies[i].mass;
            bodies[i].vel.y += bodies[i].acc.y * TIME_DELTA / bodies[i].mass;
            bodies[i].vel.z += bodies[i].acc.z * TIME_DELTA / bodies[i].mass;
            bodies[i].pos.x += bodies[i].vel.x * TIME_DELTA;
            bodies[i].pos.y += bodies[i].vel.y * TIME_DELTA;
            bodies[i].pos.z += bodies[i].vel.z * TIME_DELTA;
        }
    }

    double total_kinetic_energy = 0.0;
    for (int i = 0; i < NUM_BODIES; ++i) {
        total_kinetic_energy += 0.5 * bodies[i].mass * (bodies[i].vel.x*bodies[i].vel.x + bodies[i].vel.y*bodies[i].vel.y + bodies[i].vel.z*bodies[i].vel.z);
    }
    final_result = total_kinetic_energy;
}

void cleanup() {
    free(bodies);
    free(node_pool);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_result);

    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
