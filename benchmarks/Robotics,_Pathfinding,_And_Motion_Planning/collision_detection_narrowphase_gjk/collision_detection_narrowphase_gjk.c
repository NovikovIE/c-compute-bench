#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) implementation ---
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
// --- End of MT19937 ---

// --- Vector and GJK Data Structures ---
typedef struct {
    float x, y, z;
} Vec3;

// --- Global Benchmark Data ---
int p_num_vertices_a;
int p_num_vertices_b;
const int NUM_TESTS = 1200;
const int GJK_MAX_ITERATIONS = 32;

Vec3* vertices_a;
Vec3* vertices_b;
Vec3* initial_directions;

int final_result; // Accumulated result to prevent dead code elimination

// --- Vector Utility Functions ---
static inline float rand_float() {
    return (float)mt_rand() / (float)UINT32_MAX;
}

static inline Vec3 vec_sub(Vec3 a, Vec3 b) { return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline Vec3 vec_negate(Vec3 a) { return (Vec3){-a.x, -a.y, -a.z}; }
static inline float vec_dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline Vec3 vec_cross(Vec3 a, Vec3 b) {
    return (Vec3){
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

// --- GJK Core Logic ---

// Support function: finds the vertex on a shape furthest in a given direction.
Vec3 support(const Vec3* shape_vertices, int num_v, Vec3 dir) {
    float max_dot = -1e20f;
    Vec3 best_v = shape_vertices[0];
    for (int i = 0; i < num_v; ++i) {
        float dot = vec_dot(shape_vertices[i], dir);
        if (dot > max_dot) {
            max_dot = dot;
            best_v = shape_vertices[i];
        }
    }
    return best_v;
}

// Core GJK collision detection routine.
// It builds a simplex (triangle/tetrahedron) in the Minkowski Difference space
// and checks if it encloses the origin.
int gjk_collide(const Vec3* initial_direction) {
    Vec3 simplex[4];
    int simplex_size = 0;

    // Initial support point
    Vec3 dir = *initial_direction;
    simplex[simplex_size++] = vec_sub(support(vertices_a, p_num_vertices_a, dir), support(vertices_b, p_num_vertices_b, vec_negate(dir)));
    
    // New direction is towards the origin from the first point
    dir = vec_negate(simplex[0]);

    for (int iter = 0; iter < GJK_MAX_ITERATIONS; ++iter) {
        Vec3 a = vec_sub(support(vertices_a, p_num_vertices_a, dir), support(vertices_b, p_num_vertices_b, vec_negate(dir)));

        // If the new point is not past the origin, no collision is possible
        if (vec_dot(a, dir) < 0) {
            return 0; // No collision
        }

        simplex[simplex_size++] = a;

        // Process the simplex to see if it contains the origin
        Vec3 a_pt = simplex[simplex_size - 1];
        Vec3 ao = vec_negate(a_pt);

        if (simplex_size == 2) { // Line case
            Vec3 b_pt = simplex[0];
            Vec3 ab = vec_sub(b_pt, a_pt);
            if (vec_dot(ab, ao) > 0) {
                dir = vec_cross(vec_cross(ab, ao), ab);
            } else {
                simplex_size = 1;
                dir = ao;
            }
        } else if (simplex_size == 3) { // Triangle case
            Vec3 b_pt = simplex[1];
            Vec3 c_pt = simplex[0];
            Vec3 ab = vec_sub(b_pt, a_pt);
            Vec3 ac = vec_sub(c_pt, a_pt);
            Vec3 abc_perp = vec_cross(ab, ac);

            if (vec_dot(vec_cross(abc_perp, ac), ao) > 0) {
                if (vec_dot(ac, ao) > 0) {
                    simplex[0] = c_pt; simplex[1] = a_pt;
                    simplex_size = 2;
                    dir = vec_cross(vec_cross(ac, ao), ac);
                } else {
                    simplex[0] = b_pt; simplex[1] = a_pt;
                    simplex_size = 2;
                    goto line_ab_check;
                }
            } else {
                if (vec_dot(vec_cross(ab, abc_perp), ao) > 0) {
                    line_ab_check:
                    if (vec_dot(ab, ao) > 0) {
                        simplex[0] = b_pt; simplex[1] = a_pt;
                        simplex_size = 2;
                        dir = vec_cross(vec_cross(ab, ao), ab);
                    } else {
                        simplex[0] = a_pt;
                        simplex_size = 1;
                        dir = ao;
                    }
                } else {
                    if (vec_dot(abc_perp, ao) > 0) {
                        dir = abc_perp;
                    } else {
                        dir = vec_negate(abc_perp);
                        // swap b and c to keep winding order consistent for tetrahedron case
                        simplex[0] = b_pt; simplex[1] = c_pt;
                    }
                }
            }
        } else { // Tetrahedron case (simplex_size == 4)
            Vec3 b_pt = simplex[2];
            Vec3 c_pt = simplex[1];
            Vec3 d_pt = simplex[0];
            Vec3 ab = vec_sub(b_pt, a_pt);
            Vec3 ac = vec_sub(c_pt, a_pt);
            Vec3 ad = vec_sub(d_pt, a_pt);

            Vec3 abc_perp = vec_cross(ab, ac);
            if (vec_dot(abc_perp, ao) > 0) {
                simplex_size = 3; // a, b, c
                continue;
            }

            Vec3 acd_perp = vec_cross(ac, ad);
            if (vec_dot(acd_perp, ao) > 0) {
                simplex[2] = d_pt; // a, c, d
                simplex_size = 3;
                continue;
            }

            Vec3 adb_perp = vec_cross(ad, ab);
            if (vec_dot(adb_perp, ao) > 0) {
                simplex[1] = d_pt; // a, d, b
                simplex_size = 3;
                continue;
            }

            return 1; // Origin is inside tetrahedron
        }
    }
    return 0; // No collision found within max iterations
}


// --- Benchmark Setup, Computation, and Cleanup ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s num_vertices_a num_vertices_b seed\n", argv[0]);
        exit(1);
    }

    p_num_vertices_a = atoi(argv[1]);
    p_num_vertices_b = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    vertices_a = (Vec3*)malloc(p_num_vertices_a * sizeof(Vec3));
    vertices_b = (Vec3*)malloc(p_num_vertices_b * sizeof(Vec3));
    initial_directions = (Vec3*)malloc(NUM_TESTS * sizeof(Vec3));

    if (!vertices_a || !vertices_b || !initial_directions) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }

    // Generate vertices for two separated convex-like hulls
    // Shape A around (-2, 0, 0)
    for (int i = 0; i < p_num_vertices_a; i++) {
        vertices_a[i] = (Vec3){-2.0f + rand_float(), rand_float() - 0.5f, rand_float() - 0.5f};
    }
    // Shape B around (2, 0, 0)
    for (int i = 0; i < p_num_vertices_b; i++) {
        vertices_b[i] = (Vec3){2.0f + rand_float(), rand_float() - 0.5f, rand_float() - 0.5f};
    }

    // Generate random initial directions for the tests
    for (int i = 0; i < NUM_TESTS; i++) {
        float x = rand_float() * 2.0f - 1.0f;
        float y = rand_float() * 2.0f - 1.0f;
        float z = rand_float() * 2.0f - 1.0f;
        float mag = sqrtf(x*x + y*y + z*z);
        if (mag < 1e-6f) mag = 1.0f;
        initial_directions[i] = (Vec3){ x/mag, y/mag, z/mag };
    }
}

void run_computation() {
    int collisions = 0;
    for (int i = 0; i < NUM_TESTS; ++i) {
        collisions += gjk_collide(&initial_directions[i]);
    }
    final_result = collisions;
}

void cleanup() {
    free(vertices_a);
    free(vertices_b);
    free(initial_directions);
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
    printf("%d\n", final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
