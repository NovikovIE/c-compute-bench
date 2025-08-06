#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator ---
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
// --- End of MT19937 ---

// --- Benchmark Data Structures ---
typedef struct {
    float x, y;
} Vertex;

typedef struct {
    Vertex v[3];
} Triangle;

typedef struct {
    int num_triangles;
    int viewport_width;
    int viewport_height;
    Triangle* triangles;
    uint32_t* framebuffer;
    uint64_t final_result;
} BenchmarkData;

static BenchmarkData g_data;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_triangles viewport_width viewport_height seed\n", argv[0]);
        exit(1);
    }

    g_data.num_triangles = atoi(argv[1]);
    g_data.viewport_width = atoi(argv[2]);
    g_data.viewport_height = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    
    mt_seed(seed);

    g_data.triangles = (Triangle*)malloc(g_data.num_triangles * sizeof(Triangle));
    if (!g_data.triangles) {
        perror("Failed to allocate triangles");
        exit(1);
    }

    size_t framebuffer_size = (size_t)g_data.viewport_width * g_data.viewport_height;
    g_data.framebuffer = (uint32_t*)malloc(framebuffer_size * sizeof(uint32_t));
    if (!g_data.framebuffer) {
        perror("Failed to allocate framebuffer");
        free(g_data.triangles);
        exit(1);
    }

    for (int i = 0; i < g_data.num_triangles; ++i) {
        for (int j = 0; j < 3; ++j) {
            g_data.triangles[i].v[j].x = (float)(mt_rand() % g_data.viewport_width);
            g_data.triangles[i].v[j].y = (float)(mt_rand() % g_data.viewport_height);
        }
    }

    for (size_t i = 0; i < framebuffer_size; ++i) {
        g_data.framebuffer[i] = 0;
    }
    g_data.final_result = 0;
}

static inline float edge_function(const Vertex* a, const Vertex* b, const Vertex* p) {
    return (p->x - a->x) * (b->y - a->y) - (p->y - a->y) * (b->x - a->x);
}

void run_computation() {
    for (int i = 0; i < g_data.num_triangles; ++i) {
        const Vertex* v0 = &g_data.triangles[i].v[0];
        const Vertex* v1 = &g_data.triangles[i].v[1];
        const Vertex* v2 = &g_data.triangles[i].v[2];

        float min_x = fminf(v0->x, fminf(v1->x, v2->x));
        float max_x = fmaxf(v0->x, fmaxf(v1->x, v2->x));
        float min_y = fminf(v0->y, fminf(v1->y, v2->y));
        float max_y = fmaxf(v0->y, fmaxf(v1->y, v2->y));

        int start_x = (int)fmaxf(0.0f, floorf(min_x));
        int end_x = (int)fminf((float)g_data.viewport_width - 1.0f, ceilf(max_x));
        int start_y = (int)fmaxf(0.0f, floorf(min_y));
        int end_y = (int)fminf((float)g_data.viewport_height - 1.0f, ceilf(max_y));

        for (int y = start_y; y <= end_y; ++y) {
            for (int x = start_x; x <= end_x; ++x) {
                Vertex p = { (float)x + 0.5f, (float)y + 0.5f };

                float w0 = edge_function(v1, v2, &p);
                float w1 = edge_function(v2, v0, &p);
                float w2 = edge_function(v0, v1, &p);

                if ((w0 >= 0 && w1 >= 0 && w2 >= 0) || (w0 <= 0 && w1 <= 0 && w2 <= 0)) {
                    g_data.framebuffer[y * g_data.viewport_width + x] += (i + 1);
                }
            }
        }
    }

    uint64_t checksum = 0;
    size_t framebuffer_size = (size_t)g_data.viewport_width * g_data.viewport_height;
    for (size_t i = 0; i < framebuffer_size; ++i) {
        checksum += g_data.framebuffer[i];
    }
    g_data.final_result = checksum;
}

void cleanup() {
    free(g_data.triangles);
    free(g_data.framebuffer);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%llu\n", g_data.final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
