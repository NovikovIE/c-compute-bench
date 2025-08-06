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

// --- Benchmark-specific code --- 

#define EPSILON 1e-9

// Data Structures
typedef struct {
    double x, y;
} Point;

typedef struct {
    Point p1, p2;
    int id;
} LineSegment;

typedef enum {
    LEFT_ENDPOINT,
    RIGHT_ENDPOINT,
    INTERSECTION
} EventType;

typedef struct {
    Point p;
    EventType type;
    LineSegment *s1;
    LineSegment *s2;
} Event;

// Global data structure
struct {
    int num_line_segments;
    LineSegment* segments;
    Event* events;
    size_t num_events;
    size_t events_capacity;
    LineSegment** sweep_line;
    size_t sweep_line_size;
    size_t sweep_line_capacity;
    double current_sweep_x;
    int total_intersections;
} g_data;

// --- Utility Functions ---

int point_compare(const Point* p1, const Point* p2) {
    if (fabs(p1->x - p2->x) > EPSILON) return p1->x < p2->x ? -1 : 1;
    if (fabs(p1->y - p2->y) > EPSILON) return p1->y < p2->y ? -1 : 1;
    return 0;
}

int event_compare(const void* a, const void* b) {
    Event* e1 = (Event*)a;
    Event* e2 = (Event*)b;
    return point_compare(&e1->p, &e2->p);
}

double get_y_at_x(LineSegment* s, double x) {
    if (fabs(s->p1.x - s->p2.x) < EPSILON) return fmin(s->p1.y, s->p2.y);
    return s->p1.y + (s->p2.y - s->p1.y) * (x - s->p1.x) / (s->p2.x - s->p1.x);
}

int find_intersection(LineSegment* s1, LineSegment* s2, Point* out_p) {
    double x1 = s1->p1.x, y1 = s1->p1.y;
    double x2 = s1->p2.x, y2 = s1->p2.y;
    double x3 = s2->p1.x, y3 = s2->p1.y;
    double x4 = s2->p2.x, y4 = s2->p2.y;

    double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if (fabs(den) < EPSILON) return 0; // Parallel or collinear

    double t_num = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4);
    double u_num = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3));

    double t = t_num / den;
    double u = u_num / den;

    if (t >= -EPSILON && t <= 1.0 + EPSILON && u >= -EPSILON && u <= 1.0 + EPSILON) {
        out_p->x = x1 + t * (x2 - x1);
        out_p->y = y1 + t * (y2 - y1);
        return 1;
    }
    return 0;
}

void add_event(Event new_event) {
    Point* new_p = &new_event.p;

    // Only add event if it's strictly to the right of the sweep line
    if (new_p->x < g_data.current_sweep_x + EPSILON) return;

    // Resize event queue if necessary
    if (g_data.num_events >= g_data.events_capacity) {
        g_data.events_capacity *= 2;
        g_data.events = (Event*)realloc(g_data.events, g_data.events_capacity * sizeof(Event));
        if (!g_data.events) { perror("realloc events failed"); exit(1); }
    }

    // Insert event while maintaining sorted order (binary search + memmove)
    size_t low = 0, high = g_data.num_events;
    size_t insert_idx = g_data.num_events;
    while(low < high) {
        size_t mid = low + (high-low)/2;
        if (event_compare(&new_event, &g_data.events[mid]) < 0) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    insert_idx = low;
    
    // Don't add duplicate intersection events
    if (insert_idx > 0 && event_compare(&new_event, &g_data.events[insert_idx - 1]) == 0) return;
    if (insert_idx < g_data.num_events && event_compare(&new_event, &g_data.events[insert_idx]) == 0) return;


    memmove(&g_data.events[insert_idx + 1], &g_data.events[insert_idx], (g_data.num_events - insert_idx) * sizeof(Event));
    g_data.events[insert_idx] = new_event;
    g_data.num_events++;
}

void check_and_add_intersection(LineSegment* s1, LineSegment* s2) {
    if (!s1 || !s2) return;
    Point p;
    if (find_intersection(s1, s2, &p)) {
        Event e = {p, INTERSECTION, s1, s2};
        add_event(e);
    }
}

size_t find_segment_in_sweep_line(LineSegment* s) {
    for (size_t i = 0; i < g_data.sweep_line_size; ++i) {
        if (g_data.sweep_line[i]->id == s->id) {
            return i;
        }
    }
    return g_data.sweep_line_size; // Not found
}

void handle_left_endpoint(size_t current_event_idx) {
    LineSegment* s = g_data.events[current_event_idx].s1;
    size_t insert_pos = 0;
    double y_at_sweep = get_y_at_x(s, g_data.current_sweep_x);

    while(insert_pos < g_data.sweep_line_size && get_y_at_x(g_data.sweep_line[insert_pos], g_data.current_sweep_x) < y_at_sweep - EPSILON) {
        insert_pos++;
    }

    if (g_data.sweep_line_size >= g_data.sweep_line_capacity) {
        g_data.sweep_line_capacity *= 2;
        g_data.sweep_line = (LineSegment**)realloc(g_data.sweep_line, g_data.sweep_line_capacity * sizeof(LineSegment*));
        if (!g_data.sweep_line) { perror("realloc sweep_line failed"); exit(1); }
    }

    memmove(&g_data.sweep_line[insert_pos + 1], &g_data.sweep_line[insert_pos], (g_data.sweep_line_size - insert_pos) * sizeof(LineSegment*));
    g_data.sweep_line[insert_pos] = s;
    g_data.sweep_line_size++;

    LineSegment* above = (insert_pos > 0) ? g_data.sweep_line[insert_pos - 1] : NULL;
    LineSegment* below = (insert_pos + 1 < g_data.sweep_line_size) ? g_data.sweep_line[insert_pos + 1] : NULL;

    check_and_add_intersection(s, above);
    check_and_add_intersection(s, below);
}

void handle_right_endpoint(size_t current_event_idx) {
    LineSegment* s = g_data.events[current_event_idx].s1;
    size_t idx = find_segment_in_sweep_line(s);

    if (idx < g_data.sweep_line_size) {
        LineSegment* above = (idx > 0) ? g_data.sweep_line[idx - 1] : NULL;
        LineSegment* below = (idx + 1 < g_data.sweep_line_size) ? g_data.sweep_line[idx + 1] : NULL;
        check_and_add_intersection(above, below);
        memmove(&g_data.sweep_line[idx], &g_data.sweep_line[idx + 1], (g_data.sweep_line_size - idx - 1) * sizeof(LineSegment*));
        g_data.sweep_line_size--;
    }
}

void handle_intersection(size_t current_event_idx) {
    g_data.total_intersections++;
    Event* e = &g_data.events[current_event_idx];
    LineSegment *s1 = e->s1; LineSegment *s2 = e->s2;
    
    size_t idx1 = find_segment_in_sweep_line(s1);
    size_t idx2 = find_segment_in_sweep_line(s2);

    if (idx1 >= g_data.sweep_line_size || idx2 >= g_data.sweep_line_size) return; 

    // Ensure idx1 is the upper segment
    if (get_y_at_x(g_data.sweep_line[idx1], g_data.current_sweep_x) < get_y_at_x(g_data.sweep_line[idx2], g_data.current_sweep_x)) {
        size_t temp = idx1; idx1 = idx2; idx2 = temp;
    }
    
    LineSegment* upper_neighbor = (idx1 > 0) ? g_data.sweep_line[idx1 - 1] : NULL;
    LineSegment* lower_neighbor = (idx2 + 1 < g_data.sweep_line_size) ? g_data.sweep_line[idx2 + 1] : NULL;

    check_and_add_intersection(g_data.sweep_line[idx2], upper_neighbor);
    check_and_add_intersection(g_data.sweep_line[idx1], lower_neighbor);

    g_data.sweep_line[idx1] = s2;
    g_data.sweep_line[idx2] = s1;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_line_segments> <seed>\n", argv[0]);
        exit(1);
    }

    g_data.num_line_segments = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);
    mt_seed(seed);

    // Allocate segments
    g_data.segments = (LineSegment*)malloc(g_data.num_line_segments * sizeof(LineSegment));
    if (!g_data.segments) { perror("malloc segments failed"); exit(1); }

    // Allocate events queue (with extra space for intersections)
    g_data.num_events = g_data.num_line_segments * 2;
    g_data.events_capacity = g_data.num_events * 2; // Initial guess for capacity
    g_data.events = (Event*)malloc(g_data.events_capacity * sizeof(Event));
    if (!g_data.events) { perror("malloc events failed"); exit(1); }

    // Generate segments and initial events
    for (int i = 0; i < g_data.num_line_segments; ++i) {
        Point p1 = { (double)mt_rand() / UINT32_MAX, (double)mt_rand() / UINT32_MAX };
        Point p2 = { (double)mt_rand() / UINT32_MAX, (double)mt_rand() / UINT32_MAX };

        if (point_compare(&p1, &p2) > 0) {
            Point temp = p1; p1 = p2; p2 = temp;
        }

        g_data.segments[i] = (LineSegment){p1, p2, i};
        g_data.events[2 * i] = (Event){p1, LEFT_ENDPOINT, &g_data.segments[i], NULL};
        g_data.events[2 * i + 1] = (Event){p2, RIGHT_ENDPOINT, &g_data.segments[i], NULL};
    }

    // Sort initial events
    qsort(g_data.events, g_data.num_events, sizeof(Event), event_compare);

    g_data.total_intersections = 0;
    g_data.sweep_line = NULL;
}

void run_computation() {
    // Allocate sweep line
    g_data.sweep_line_capacity = g_data.num_line_segments > 16 ? 16 : g_data.num_line_segments;
    g_data.sweep_line = (LineSegment**)malloc(g_data.sweep_line_capacity * sizeof(LineSegment*));
    if (!g_data.sweep_line && g_data.sweep_line_capacity > 0) { perror("malloc sweep line failed"); exit(1); }
    g_data.sweep_line_size = 0;

    size_t event_idx = 0;
    while(event_idx < g_data.num_events) {
        Event current_event = g_data.events[event_idx];
        g_data.current_sweep_x = current_event.p.x;

        // To handle vertical clusters of events, process all events at this same point together
        size_t end_of_cluster = event_idx;
        while(end_of_cluster + 1 < g_data.num_events && point_compare(&g_data.events[event_idx].p, &g_data.events[end_of_cluster+1].p) == 0) {
            end_of_cluster++;
        }

        for (size_t i = event_idx; i <= end_of_cluster; ++i) {
             if (g_data.events[i].type == RIGHT_ENDPOINT) handle_right_endpoint(i);
        }
        for (size_t i = event_idx; i <= end_of_cluster; ++i) {
            if (g_data.events[i].type == LEFT_ENDPOINT) handle_left_endpoint(i);
            else if (g_data.events[i].type == INTERSECTION) handle_intersection(i);
        }
        
        event_idx = end_of_cluster + 1;
    }
}

void cleanup() {
    free(g_data.segments);
    free(g_data.events);
    free(g_data.sweep_line);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%d\n", g_data.total_intersections);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}