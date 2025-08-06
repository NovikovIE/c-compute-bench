#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Mersenne Twister Generator (verbatim)
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

double mt_rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// --- Data Structures ---
typedef struct { double x, y, z; } Vec3;
typedef struct { double m[4][4]; } Mat4;
typedef struct { int v1, v2; } Edge;

// --- Global Benchmark Data ---
int initial_vertex_count;
int target_vertex_count;
int triangle_count;
int edge_count;

Vec3* vertices_pos;
Mat4* vertices_q;
int* vertices_active;
int* triangles;
Edge* edges;

double final_checksum;

// --- Helper Functions ---

static Vec3 vec3_sub(Vec3 a, Vec3 b) { return (Vec3){a.x - b.x, a.y - b.y, a.z - b.z}; }
static Vec3 vec3_cross(Vec3 a, Vec3 b) { return (Vec3){a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
static double vec3_dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static double vec3_length(Vec3 v) { return sqrt(vec3_dot(v, v)); }
static Vec3 vec3_normalize(Vec3 v) {
    double len = vec3_length(v);
    if (len > 1e-9) {
        double inv_len = 1.0 / len;
        return (Vec3){v.x * inv_len, v.y * inv_len, v.z * inv_len};
    }
    return (Vec3){0, 0, 0};
}

static void mat4_add(Mat4* out, const Mat4* a, const Mat4* b) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            out->m[i][j] = a->m[i][j] + b->m[i][j];
        }
    }
}

static int invert_3x3(double m_inv[3][3], const Mat4* q) {
    const double (*m)[4] = q->m;
    double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                 m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                 m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if (fabs(det) < 1e-12) {
        return 0; // Not invertible
    }

    double inv_det = 1.0 / det;
    m_inv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * inv_det;
    m_inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det;
    m_inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;
    m_inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det;
    m_inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
    m_inv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * inv_det;
    m_inv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * inv_det;
    m_inv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * inv_det;
    m_inv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * inv_det;
    return 1;
}

static double calculate_error(const Mat4* q, Vec3 p) {
    return q->m[0][0]*p.x*p.x + 2*q->m[0][1]*p.x*p.y + 2*q->m[0][2]*p.x*p.z + 2*q->m[0][3]*p.x +
           q->m[1][1]*p.y*p.y + 2*q->m[1][2]*p.y*p.z + 2*q->m[1][3]*p.y +
           q->m[2][2]*p.z*p.z + 2*q->m[2][3]*p.z +
           q->m[3][3];
}

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <initial_vertex_count> <target_vertex_count> <seed>\n", argv[0]);
        exit(1);
    }
    initial_vertex_count = atoi(argv[1]);
    target_vertex_count = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (initial_vertex_count <= target_vertex_count || target_vertex_count < 3) {
        fprintf(stderr, "Invalid vertex counts.\n");
        exit(1);
    }

    mt_seed(seed);

    vertices_pos = (Vec3*)malloc(initial_vertex_count * sizeof(Vec3));
    vertices_q = (Mat4*)malloc(initial_vertex_count * sizeof(Mat4));
    vertices_active = (int*)malloc(initial_vertex_count * sizeof(int));
    memset(vertices_q, 0, initial_vertex_count * sizeof(Mat4));

    // Generate vertices on a noisy sphere
    for (int i = 0; i < initial_vertex_count; ++i) {
        double u = mt_rand_double();
        double v = mt_rand_double();
        double theta = 2.0 * M_PI * u;
        double phi = acos(2.0 * v - 1.0);
        double r = 1.0 + (mt_rand_double() - 0.5) * 0.1; // Add some noise
        vertices_pos[i].x = r * sin(phi) * cos(theta);
        vertices_pos[i].y = r * sin(phi) * sin(theta);
        vertices_pos[i].z = r * cos(phi);
        vertices_active[i] = 1;
    }

    // Generate random triangles to calculate initial quadrics
    triangle_count = initial_vertex_count * 2;
    triangles = (int*)malloc(triangle_count * 3 * sizeof(int));
    for (int i = 0; i < triangle_count; ++i) {
        int v1 = mt_rand() % initial_vertex_count;
        int v2 = mt_rand() % initial_vertex_count;
        int v3 = mt_rand() % initial_vertex_count;
        triangles[i*3 + 0] = v1;
        triangles[i*3 + 1] = v2;
        triangles[i*3 + 2] = v3;

        Vec3 p1 = vertices_pos[v1];
        Vec3 p2 = vertices_pos[v2];
        Vec3 p3 = vertices_pos[v3];
        Vec3 norm = vec3_normalize(vec3_cross(vec3_sub(p2, p1), vec3_sub(p3, p1)));
        
        double a = norm.x, b = norm.y, c = norm.z;
        double d = -vec3_dot(norm, p1);

        Mat4 Kp = {{{a*a, a*b, a*c, a*d},{a*b, b*b, b*c, b*d},{a*c, b*c, c*c, c*d},{a*d, b*d, c*d, d*d}}};
        
        mat4_add(&vertices_q[v1], &vertices_q[v1], &Kp);
        mat4_add(&vertices_q[v2], &vertices_q[v2], &Kp);
        mat4_add(&vertices_q[v3], &vertices_q[v3], &Kp);
    }
    
    // Create a list of candidate edges from triangles
    edge_count = triangle_count * 3;
    edges = (Edge*)malloc(edge_count * sizeof(Edge));
    for (int i = 0; i < triangle_count; ++i) {
        int v0 = triangles[i*3 + 0], v1 = triangles[i*3 + 1], v2 = triangles[i*3 + 2];
        edges[i*3 + 0] = (Edge){v0, v1};
        edges[i*3 + 1] = (Edge){v1, v2};
        edges[i*3 + 2] = (Edge){v2, v0};
    }
}

void run_computation() {
    int current_vertex_count = initial_vertex_count;

    while (current_vertex_count > target_vertex_count) {
        double min_error = 1e30;
        int best_v1 = -1, best_v2 = -1;
        Vec3 optimal_pos = {0,0,0};

        for (int i = 0; i < edge_count; ++i) {
            int v1_idx = edges[i].v1;
            int v2_idx = edges[i].v2;

            if (!vertices_active[v1_idx] || !vertices_active[v2_idx] || v1_idx == v2_idx) {
                continue;
            }

            Mat4 q_sum;
            mat4_add(&q_sum, &vertices_q[v1_idx], &vertices_q[v2_idx]);
            
            Vec3 p;
            double m_inv[3][3];
            if (invert_3x3(m_inv, &q_sum)) {
                p.x = -m_inv[0][0] * q_sum.m[0][3] - m_inv[0][1] * q_sum.m[1][3] - m_inv[0][2] * q_sum.m[2][3];
                p.y = -m_inv[1][0] * q_sum.m[0][3] - m_inv[1][1] * q_sum.m[1][3] - m_inv[1][2] * q_sum.m[2][3];
                p.z = -m_inv[2][0] * q_sum.m[0][3] - m_inv[2][1] * q_sum.m[1][3] - m_inv[2][2] * q_sum.m[2][3];
            } else {
                Vec3 p1 = vertices_pos[v1_idx];
                Vec3 p2 = vertices_pos[v2_idx];
                p.x = (p1.x + p2.x) * 0.5;
                p.y = (p1.y + p2.y) * 0.5;
                p.z = (p1.z + p2.z) * 0.5;
            }

            double error = calculate_error(&q_sum, p);
            if (error < min_error) {
                min_error = error;
                best_v1 = v1_idx;
                best_v2 = v2_idx;
                optimal_pos = p;
            }
        }

        if (best_v1 == -1) break; // No more valid pairs to collapse

        // Perform the collapse: v2 -> v1
        mat4_add(&vertices_q[best_v1], &vertices_q[best_v1], &vertices_q[best_v2]);
        vertices_pos[best_v1] = optimal_pos;
        vertices_active[best_v2] = 0;

        current_vertex_count--;
    }

    final_checksum = 0.0;
    for (int i = 0; i < initial_vertex_count; ++i) {
        if (vertices_active[i]) {
            final_checksum += vertices_pos[i].x + vertices_pos[i].y + vertices_pos[i].z;
        }
    }
}

void cleanup() {
    free(vertices_pos);
    free(vertices_q);
    free(vertices_active);
    free(triangles);
    free(edges);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%f\n", final_checksum);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
