#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

/* --- MERSENNE TWISTER (Verbatim) --- */
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
/* --- END MERSENNE TWISTER --- */


/* --- BENCHMARK DATA AND GLOBALS --- */
typedef struct {
    double* coords;
} Point;

static int num_points;
static int dimensions;
static Point* points;

// The final result to be printed to stdout
static double final_min_dist;

/* --- HELPER FUNCTIONS --- */

// Calculate squared Euclidean distance between two points
static double dist_sq(const Point* p1, const Point* p2) {
    double sum_sq = 0.0;
    for (int i = 0; i < dimensions; i++) {
        double diff = p1->coords[i] - p2->coords[i];
        sum_sq += diff * diff;
    }
    return sum_sq;
}

// qsort comparison function for the first dimension
static int compare_x(const void* a, const void* b) {
    Point* p1 = (Point*)a;
    Point* p2 = (Point*)b;
    if (p1->coords[0] < p2->coords[0]) return -1;
    if (p1->coords[0] > p2->coords[0]) return 1;
    return 0;
}

// qsort comparison function for the second dimension
static int compare_y(const void* a, const void* b) {
    Point* p1 = (Point*)a;
    Point* p2 = (Point*)b;
    if (p1->coords[1] < p2->coords[1]) return -1;
    if (p1->coords[1] > p2->coords[1]) return 1;
    return 0;
}

// Brute-force method to find the smallest distance in a small set of points
static double brute_force(Point p[], int n) {
    double min_d_sq = DBL_MAX;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d_sq = dist_sq(&p[i], &p[j]);
            if (d_sq < min_d_sq) {
                min_d_sq = d_sq;
            }
        }
    }
    return min_d_sq;
}

// Find the minimum distance in a strip of points
static double strip_closest(Point strip[], int size, double d_sq) {
    double min_d_sq = d_sq;

    // Sort strip by the second dimension. For dimensions > 2, this is a heuristic.
    if (dimensions > 1) {
        qsort(strip, size, sizeof(Point), compare_y);
    }
    
    for (int i = 0; i < size; ++i) {
        for (int j = i + 1; j < size; ++j) {
            // Optimization: if the y-distance is already greater than the current min,
            // no need to check this point or any subsequent points in the sorted strip.
            if (dimensions > 1 && (strip[j].coords[1] - strip[i].coords[1]) * (strip[j].coords[1] - strip[i].coords[1]) >= min_d_sq) {
                break;
            }
            double current_dist_sq = dist_sq(&strip[i], &strip[j]);
            if (current_dist_sq < min_d_sq) {
                min_d_sq = current_dist_sq;
            }
        }
    }
    return min_d_sq;
}

// The main recursive divide-and-conquer function
static double closest_util(Point p_sorted_x[], int n) {
    if (n <= 3) {
        return brute_force(p_sorted_x, n);
    }

    int mid = n / 2;
    Point mid_point = p_sorted_x[mid];

    double dl_sq = closest_util(p_sorted_x, mid);
    double dr_sq = closest_util(p_sorted_x + mid, n - mid);
    double d_sq = fmin(dl_sq, dr_sq);

    Point* strip = (Point*)malloc(n * sizeof(Point));
    if (!strip) {
        perror("malloc strip failed");
        exit(1);
    }
    int j = 0;
    for (int i = 0; i < n; i++) {
        double x_diff = p_sorted_x[i].coords[0] - mid_point.coords[0];
        if (x_diff * x_diff < d_sq) {
            strip[j] = p_sorted_x[i];
            j++;
        }
    }

    double strip_min_sq = strip_closest(strip, j, d_sq);
    free(strip);

    return fmin(d_sq, strip_min_sq);
}


/* --- BENCHMARK FUNCTIONS --- */

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_points> <dimensions> <seed>\n", argv[0]);
        exit(1);
    }

    num_points = atoi(argv[1]);
    dimensions = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);
    
    if (num_points <= 1 || dimensions <= 0) {
        fprintf(stderr, "Error: num_points must be > 1 and dimensions must be > 0.\n");
        exit(1);
    }
    
    mt_seed(seed);

    points = (Point*) malloc(num_points * sizeof(Point));
    if (!points) {
        perror("Failed to allocate points array");
        exit(1);
    }

    for (int i = 0; i < num_points; i++) {
        points[i].coords = (double*) malloc(dimensions * sizeof(double));
        if (!points[i].coords) {
            perror("Failed to allocate coordinates array");
            for (int k = 0; k < i; k++) {
                free(points[k].coords);
            }
            free(points);
            exit(1);
        }
        for (int j = 0; j < dimensions; j++) {
            points[i].coords[j] = (double)mt_rand() / (double)UINT32_MAX * 10000.0;
        }
    }
}

void run_computation() {
    qsort(points, num_points, sizeof(Point), compare_x);

    double min_dist_squared = closest_util(points, num_points);

    final_min_dist = sqrt(min_dist_squared);
}

void cleanup() {
    for (int i = 0; i < num_points; i++) {
        free(points[i].coords);
    }
    free(points);
}


/* --- MAIN DRIVER --- */

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    
    printf("%f\n", final_min_dist);

    fprintf(stderr, "%.6f", time_taken);

    cleanup();

    return 0;
}
