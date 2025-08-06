#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

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

// --- BENCHMARK SPECIFIC CODE --- 

#define EPSILON 1e-9

// --- Data Structures --- 
typedef struct {
    double x, y, z;
} Point3D;

typedef struct {
    int v[3]; // Vertex indices
    Point3D normal;
    double offset;
    int visible;
    int* outside_set; 
    int outside_set_count;
    int outside_set_capacity;
} Face;

typedef struct {
    int u, v; // vertex indices for an edge
} Edge;

// --- Global Data --- 
int g_num_points;
Point3D *g_points = NULL;
int g_final_result = 0;

// --- Geometric Helper Functions --- 
Point3D vec_sub(Point3D a, Point3D b) { return (Point3D){a.x - b.x, a.y - b.y, a.z - b.z}; }
Point3D cross_product(Point3D a, Point3D b) { return (Point3D){a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x}; }
double dot_product(Point3D a, Point3D b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
double vec_len_sq(Point3D a) { return dot_product(a, a); }
void normalize(Point3D* v) {
    double len = sqrt(dot_product(*v, *v));
    if (len > EPSILON) { v->x /= len; v->y /= len; v->z /= len; }
}
double dist_point_plane(Point3D p, Face f) { return dot_product(f.normal, p) + f.offset; }

// --- Quickhull Helper Functions --- 

// qsort comparison for edges
int compare_edges(const void* a, const void* b) {
    const Edge* e1 = (const Edge*)a;
    const Edge* e2 = (const Edge*)b;
    if (e1->u < e2->u) return -1;
    if (e1->u > e2->u) return 1;
    if (e1->v < e2->v) return -1;
    if (e1->v > e2->v) return 1;
    return 0;
}

void add_point_to_face(Face* f, int point_idx) {
    if (f->outside_set_count >= f->outside_set_capacity) {
        f->outside_set_capacity = f->outside_set_capacity == 0 ? 64 : f->outside_set_capacity * 2;
        f->outside_set = (int*)realloc(f->outside_set, f->outside_set_capacity * sizeof(int));
    }
    f->outside_set[f->outside_set_count++] = point_idx;
}

// --- Benchmark Functions --- 
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s num_points seed\n", argv[0]);
        exit(1);
    }
    g_num_points = atoi(argv[1]);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);

    if (g_num_points <= 4) {
        fprintf(stderr, "FATAL: num_points must be > 4\n");
        exit(1);
    }

    mt_seed(seed);

    g_points = (Point3D*)malloc(g_num_points * sizeof(Point3D));
    if (!g_points) {
        perror("Failed to allocate memory for points");
        exit(1);
    }

    // Generate points in a unit cube
    for (int i = 0; i < g_num_points; ++i) {
        g_points[i].x = (double)mt_rand() / UINT32_MAX;
        g_points[i].y = (double)mt_rand() / UINT32_MAX;
        g_points[i].z = (double)mt_rand() / UINT32_MAX;
    }
}

void run_computation() {
    // --- Data structures for the algorithm --- 
    Face* faces = NULL;
    int face_count = 0;
    int face_capacity = 128;
    faces = (Face*)calloc(face_capacity, sizeof(Face));

    // --- 1. Find initial tetrahedron --- 
    int initial_verts[4];
    // Find extreme points along axes
    int extremes[6] = {0, 0, 0, 0, 0, 0};
    for (int i = 1; i < g_num_points; ++i) {
        if (g_points[i].x < g_points[extremes[0]].x) extremes[0] = i;
        if (g_points[i].x > g_points[extremes[1]].x) extremes[1] = i;
        if (g_points[i].y < g_points[extremes[2]].y) extremes[2] = i;
        if (g_points[i].y > g_points[extremes[3]].y) extremes[3] = i;
        if (g_points[i].z < g_points[extremes[4]].z) extremes[4] = i;
        if (g_points[i].z > g_points[extremes[5]].z) extremes[5] = i;
    }

    // Find two most distant points among extremes
    double max_dist_sq = -1;
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            double d = vec_len_sq(vec_sub(g_points[extremes[i]], g_points[extremes[j]]));
            if (d > max_dist_sq) {
                max_dist_sq = d;
                initial_verts[0] = extremes[i];
                initial_verts[1] = extremes[j];
            }
        }
    }

    // Find third point maximizing distance to the line
    max_dist_sq = -1;
    for (int i = 0; i < g_num_points; ++i) {
        Point3D v_line = vec_sub(g_points[initial_verts[1]], g_points[initial_verts[0]]);
        Point3D v_p = vec_sub(g_points[i], g_points[initial_verts[0]]);
        double d = vec_len_sq(cross_product(v_line, v_p));
        if (d > max_dist_sq) {
            max_dist_sq = d;
            initial_verts[2] = i;
        }
    }

    // Find fourth point maximizing distance to the plane
    Point3D v1 = vec_sub(g_points[initial_verts[1]], g_points[initial_verts[0]]);
    Point3D v2 = vec_sub(g_points[initial_verts[2]], g_points[initial_verts[0]]);
    Point3D normal = cross_product(v1, v2);
    normalize(&normal);
    double offset = -dot_product(normal, g_points[initial_verts[0]]);
    max_dist_sq = -1;
    for(int i = 0; i < g_num_points; ++i) {
        double d = fabs(dot_product(normal, g_points[i]) + offset);
        if (d > max_dist_sq) {
            max_dist_sq = d;
            initial_verts[3] = i;
        }
    }

    // --- 2. Build initial hull --- 
    int temp_faces[4][3] = { {0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {1, 3, 2} };

    for (int i = 0; i < 4; ++i) {
        int p_idx[3] = { initial_verts[temp_faces[i][0]], initial_verts[temp_faces[i][1]], initial_verts[temp_faces[i][2]] };
        Point3D p0 = g_points[p_idx[0]];
        Point3D p1 = g_points[p_idx[1]];
        Point3D p2 = g_points[p_idx[2]];
        
        Face* f = &faces[face_count++];
        f->normal = cross_product(vec_sub(p1, p0), vec_sub(p2, p0));
        normalize(&f->normal);
        f->offset = -dot_product(f->normal, p0);
        
        // Orient face outwards. Check against the opposite vertex
        int other_v_idx = initial_verts[temp_faces[i ^ 1][2]]; // A simple way to get an opposite vertex
        if (dist_point_plane(g_points[other_v_idx], *f) > 0) {
            f->normal.x *= -1; f->normal.y *= -1; f->normal.z *= -1;
            f->offset *= -1;
            f->v[0] = p_idx[0]; f->v[1] = p_idx[2]; f->v[2] = p_idx[1];
        } else {
            f->v[0] = p_idx[0]; f->v[1] = p_idx[1]; f->v[2] = p_idx[2];
        }
    }

    // --- 3. Partition remaining points
    for (int i = 0; i < g_num_points; ++i) {
        int is_initial = 0;
        for(int j=0; j<4; ++j) if(i == initial_verts[j]) { is_initial=1; break; }
        if(is_initial) continue;

        int best_face_idx = -1;
        double max_dist = EPSILON;

        for (int j = 0; j < face_count; ++j) {
            double d = dist_point_plane(g_points[i], faces[j]);
            if (d > max_dist) {
                max_dist = d;
                best_face_idx = j;
            }
        }
        if (best_face_idx != -1) {
            add_point_to_face(&faces[best_face_idx], i);
        }
    }
    
    // --- 4. Main loop --- 
    int active_face_idx = -1;
    while(1) {
        active_face_idx = -1;
        for(int i = 0; i < face_count; ++i) {
            if(faces[i].outside_set_count > 0) {
                active_face_idx = i;
                break;
            }
        }
        if (active_face_idx == -1) break; // Done

        Face* active_face = &faces[active_face_idx];

        // Find eye point
        int eye_point_idx = -1;
        double max_dist = 0;
        for (int i = 0; i < active_face->outside_set_count; ++i) {
            double d = dist_point_plane(g_points[active_face->outside_set[i]], *active_face);
            if (d > max_dist) {
                max_dist = d;
                eye_point_idx = active_face->outside_set[i];
            }
        }

        // Find visible faces
        int* visible_face_indices = (int*)malloc(face_count * sizeof(int));
        int visible_count = 0;
        int* candidate_points = (int*)malloc(g_num_points * sizeof(int));
        int candidate_count = 0;

        for(int i = 0; i < face_count; ++i) {
            if (dist_point_plane(g_points[eye_point_idx], faces[i]) > EPSILON) {
                faces[i].visible = 1;
                visible_face_indices[visible_count++] = i;
                for(int j = 0; j < faces[i].outside_set_count; ++j) {
                    candidate_points[candidate_count++] = faces[i].outside_set[j];
                }
                free(faces[i].outside_set);
                faces[i].outside_set = NULL;
                faces[i].outside_set_count = 0;
            }
        }
        
        // Find horizon edges
        Edge* horizon_edges = (Edge*)malloc(visible_count * 3 * sizeof(Edge));
        int horizon_edge_count = 0;
        for(int i=0; i<visible_count; ++i) {
            Face* f = &faces[visible_face_indices[i]];
            int v[3] = {f->v[0], f->v[1], f->v[2]};
            Edge edges[3] = {{v[0], v[1]}, {v[1], v[2]}, {v[2], v[0]}};
            for(int j=0; j<3; ++j) {
              int u = edges[j].u < edges[j].v ? edges[j].u : edges[j].v;
              int v = edges[j].u > edges[j].v ? edges[j].u : edges[j].v;
              horizon_edges[horizon_edge_count++] = (Edge){u, v};
            }
        }
        qsort(horizon_edges, horizon_edge_count, sizeof(Edge), compare_edges);
        
        Edge* final_horizon = (Edge*)malloc(horizon_edge_count * sizeof(Edge));
        int final_horizon_count = 0;
        if (horizon_edge_count > 0) {
            for (int i = 0; i < horizon_edge_count; ) {
                int j = i + 1;
                while (j < horizon_edge_count && compare_edges(&horizon_edges[i], &horizon_edges[j]) == 0) j++;
                if ((j - i) == 1) {
                    final_horizon[final_horizon_count++] = horizon_edges[i];
                }
                i = j;
            }
        }
        free(horizon_edges);

        // Remove visible faces
        int new_face_count = 0;
        for(int i=0; i<face_count; ++i) {
            if(!faces[i].visible) {
                if (i != new_face_count) {
                    faces[new_face_count] = faces[i];
                }
                new_face_count++;
            } 
        }
        face_count = new_face_count;

        // Add new faces from horizon
        for (int i = 0; i < final_horizon_count; ++i) {
            if (face_count >= face_capacity) {
                face_capacity *= 2;
                faces = (Face*)realloc(faces, face_capacity * sizeof(Face));
            }
            Face* f = &faces[face_count++];
            memset(f, 0, sizeof(Face));
            f->v[0] = final_horizon[i].u; f->v[1] = final_horizon[i].v; f->v[2] = eye_point_idx;
            f->normal = cross_product(vec_sub(g_points[f->v[1]], g_points[f->v[0]]), vec_sub(g_points[f->v[2]], g_points[f->v[0]]));
            normalize(&f->normal);
            f->offset = -dot_product(f->normal, g_points[f->v[0]]);
        }
        free(final_horizon);

        // Re-assign candidate points
        for (int i = 0; i < candidate_count; ++i) {
            int p_idx = candidate_points[i];
            if (p_idx == eye_point_idx) continue;
            int best_face_idx = -1;
            double max_dist = EPSILON;
            for (int j = 0; j < face_count; ++j) {
                double d = dist_point_plane(g_points[p_idx], faces[j]);
                if (d > max_dist) {
                    max_dist = d;
                    best_face_idx = j;
                }
            }
            if (best_face_idx != -1) {
                add_point_to_face(&faces[best_face_idx], p_idx);
            }
        }

        free(candidate_points);
        free(visible_face_indices);
    }

    g_final_result = face_count;

    for (int i = 0; i < face_count; ++i) {
        free(faces[i].outside_set);
    }
    free(faces);
}

void cleanup() {
    if (g_points) {
        free(g_points);
        g_points = NULL;
    }
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", g_final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
