#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- Mersenne Twister (MT19937) Generator --- Do Not Modify ---
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
// --- End of Mersenne Twister ---

// --- Benchmark Data Structures ---
typedef struct {
    float x, y, z;
} Vertex;

typedef struct {
    int v[3];
} Face;

typedef struct {
    Vertex* vertices;
    int num_vertices;
    int capacity_vertices;

    Face* faces;
    int num_faces;
    int capacity_faces;
} Mesh;

// --- Simple Hash Map for Edge -> New Vertex Index Mapping ---
// Key: 64-bit integer combining two vertex indices. Value: index of the new midpoint vertex.
#define HM_EMPTY_KEY 0xFFFFFFFFFFFFFFFFULL

typedef struct {
    uint64_t key;
    int value;
} EdgeMapEntry;

typedef struct {
    EdgeMapEntry* entries;
    int capacity;
    int mask;
} EdgeMap;

static inline uint64_t create_edge_key(int v1, int v2) {
    if (v1 < v2) return ((uint64_t)v1 << 32) | (uint32_t)v2;
    return ((uint64_t)v2 << 32) | (uint32_t)v1;
}

static uint32_t hash_fnv1a(uint64_t key) {
    const uint64_t FNV_PRIME = 1099511628211ULL;
    const uint64_t FNV_OFFSET_BASIS = 14695981039346656037ULL;
    uint64_t hash = FNV_OFFSET_BASIS;
    for (int i = 0; i < 8; ++i) {
        hash ^= (key >> (i * 8)) & 0xFF;
        hash *= FNV_PRIME;
    }
    return (uint32_t)(hash ^ (hash >> 32));
}

EdgeMap* edge_map_create(int initial_capacity) {
    EdgeMap* map = (EdgeMap*)malloc(sizeof(EdgeMap));
    if (!map) return NULL;
    map->capacity = 1;
    while (map->capacity < initial_capacity) {
        map->capacity <<= 1;
    }
    map->mask = map->capacity - 1;
    map->entries = (EdgeMapEntry*)malloc(map->capacity * sizeof(EdgeMapEntry));
    if (!map->entries) {
        free(map);
        return NULL;
    }
    for (int i = 0; i < map->capacity; ++i) {
        map->entries[i].key = HM_EMPTY_KEY;
    }
    return map;
}

void edge_map_destroy(EdgeMap* map) {
    if (map) {
        free(map->entries);
        free(map);
    }
}

int edge_map_get(EdgeMap* map, uint64_t key) {
    uint32_t index = hash_fnv1a(key) & map->mask;
    while (map->entries[index].key != HM_EMPTY_KEY) {
        if (map->entries[index].key == key) {
            return map->entries[index].value;
        }
        index = (index + 1) & map->mask;
    }
    return -1; // Not found
}

void edge_map_put(EdgeMap* map, uint64_t key, int value) {
    uint32_t index = hash_fnv1a(key) & map->mask;
    while (map->entries[index].key != HM_EMPTY_KEY) {
        if (map->entries[index].key == key) {
            map->entries[index].value = value; // Update existing
            return;
        }
        index = (index + 1) & map->mask;
    }
    map->entries[index].key = key;
    map->entries[index].value = value;
}

// --- Global Benchmark State ---
static int g_num_subdivision_levels = 0;
static Mesh* g_mesh = NULL;
static float g_final_checksum = 0.0f;

// --- Benchmark Functions ---
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_vertices> <num_subdivision_levels> <seed>\n", argv[0]);
        exit(1);
    }

    int num_vertices = atoi(argv[1]);
    g_num_subdivision_levels = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    g_mesh = (Mesh*)malloc(sizeof(Mesh));
    if (!g_mesh) { exit(1); }

    g_mesh->num_vertices = num_vertices;
    g_mesh->capacity_vertices = num_vertices;
    g_mesh->vertices = (Vertex*)malloc(g_mesh->capacity_vertices * sizeof(Vertex));
    if (!g_mesh->vertices) { exit(1); }

    for (int i = 0; i < g_mesh->num_vertices; ++i) {
        g_mesh->vertices[i].x = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
        g_mesh->vertices[i].y = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
        g_mesh->vertices[i].z = ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
    }

    int num_faces = num_vertices * 2; 
    g_mesh->num_faces = num_faces;
    g_mesh->capacity_faces = num_faces;
    g_mesh->faces = (Face*)malloc(g_mesh->capacity_faces * sizeof(Face));
    if (!g_mesh->faces) { exit(1); }
    
    for (int i = 0; i < g_mesh->num_faces; ++i) {
        g_mesh->faces[i].v[0] = mt_rand() % g_mesh->num_vertices;
        g_mesh->faces[i].v[1] = mt_rand() % g_mesh->num_vertices;
        g_mesh->faces[i].v[2] = mt_rand() % g_mesh->num_vertices;
    }
}

void run_computation() {
    for (int level = 0; level < g_num_subdivision_levels; ++level) {
        Mesh* old_mesh = g_mesh;
        Mesh* new_mesh = (Mesh*)malloc(sizeof(Mesh));

        // Estimate new sizes. This is an over-estimate for vertices but safe.
        int num_edges_estimate = old_mesh->num_faces * 3;
        new_mesh->capacity_vertices = old_mesh->num_vertices + num_edges_estimate;
        new_mesh->capacity_faces = old_mesh->num_faces * 4;

        new_mesh->vertices = (Vertex*)malloc(new_mesh->capacity_vertices * sizeof(Vertex));
        new_mesh->faces = (Face*)malloc(new_mesh->capacity_faces * sizeof(Face));

        // Copy old vertices
        memcpy(new_mesh->vertices, old_mesh->vertices, old_mesh->num_vertices * sizeof(Vertex));
        new_mesh->num_vertices = old_mesh->num_vertices;
        new_mesh->num_faces = 0;

        // This simplified subdivision splits each edge at its midpoint.
        // A hash map is used to ensure each edge is split only once.
        EdgeMap* edge_map = edge_map_create(num_edges_estimate);
        
        for (int i = 0; i < old_mesh->num_faces; ++i) {
            int v_indices[3] = {old_mesh->faces[i].v[0], old_mesh->faces[i].v[1], old_mesh->faces[i].v[2]};
            int midpoint_indices[3];

            for (int j = 0; j < 3; ++j) {
                int v1_idx = v_indices[j];
                int v2_idx = v_indices[(j + 1) % 3];
                uint64_t edge_key = create_edge_key(v1_idx, v2_idx);

                int midpoint_idx = edge_map_get(edge_map, edge_key);
                if (midpoint_idx == -1) {
                    Vertex* v1 = &old_mesh->vertices[v1_idx];
                    Vertex* v2 = &old_mesh->vertices[v2_idx];
                    midpoint_idx = new_mesh->num_vertices++;
                    new_mesh->vertices[midpoint_idx].x = (v1->x + v2->x) * 0.5f;
                    new_mesh->vertices[midpoint_idx].y = (v1->y + v2->y) * 0.5f;
                    new_mesh->vertices[midpoint_idx].z = (v1->z + v2->z) * 0.5f;
                    edge_map_put(edge_map, edge_key, midpoint_idx);
                }
                midpoint_indices[j] = midpoint_idx;
            }
            
            // Create 4 new faces for the subdivided triangle
            // Face 1: v0, m01, m20
            new_mesh->faces[new_mesh->num_faces].v[0] = v_indices[0];
            new_mesh->faces[new_mesh->num_faces].v[1] = midpoint_indices[0];
            new_mesh->faces[new_mesh->num_faces].v[2] = midpoint_indices[2];
            new_mesh->num_faces++;

            // Face 2: v1, m12, m01
            new_mesh->faces[new_mesh->num_faces].v[0] = v_indices[1];
            new_mesh->faces[new_mesh->num_faces].v[1] = midpoint_indices[1];
            new_mesh->faces[new_mesh->num_faces].v[2] = midpoint_indices[0];
            new_mesh->num_faces++;
            
            // Face 3: v2, m20, m12
            new_mesh->faces[new_mesh->num_faces].v[0] = v_indices[2];
            new_mesh->faces[new_mesh->num_faces].v[1] = midpoint_indices[2];
            new_mesh->faces[new_mesh->num_faces].v[2] = midpoint_indices[1];
            new_mesh->num_faces++;

            // Face 4: m01, m12, m20 (center)
            new_mesh->faces[new_mesh->num_faces].v[0] = midpoint_indices[0];
            new_mesh->faces[new_mesh->num_faces].v[1] = midpoint_indices[1];
            new_mesh->faces[new_mesh->num_faces].v[2] = midpoint_indices[2];
            new_mesh->num_faces++;
        }

        // Cleanup old mesh and replace with the new one
        edge_map_destroy(edge_map);
        free(old_mesh->vertices);
        free(old_mesh->faces);
        free(old_mesh);
        g_mesh = new_mesh;
    }

    // Calculate a checksum to prevent dead code elimination
    g_final_checksum = 0.0f;
    for (int i = 0; i < g_mesh->num_vertices; ++i) {
        g_final_checksum += g_mesh->vertices[i].x + g_mesh->vertices[i].y + g_mesh->vertices[i].z;
    }
}

void cleanup() {
    if (g_mesh) {
        if (g_mesh->vertices) free(g_mesh->vertices);
        if (g_mesh->faces) free(g_mesh->faces);
        free(g_mesh);
    }
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", g_final_checksum);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
