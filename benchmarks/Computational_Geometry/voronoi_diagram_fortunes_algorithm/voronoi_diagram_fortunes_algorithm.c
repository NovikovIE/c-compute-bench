/*
 * Benchmark: voronoi_diagram_fortunes_algorithm
 * 
 * Description: This program implements Fortune's sweep-line algorithm to compute the
 * Voronoi diagram for a given set of 2D points (sites). The Voronoi diagram partitions
 * the plane into regions, where each region consists of points closer to one particular
 * site than to any other. This algorithm is fundamental in computational geometry and has
 * applications in areas like computer graphics, geographic information systems, and robotics.
 *
 * Algorithm Overview:
 * 1.  A vertical sweep-line moves across the plane (top to bottom in this implementation).
 * 2.  An "event queue" (a priority queue) stores two types of events:
 *     - Site events: When the sweep-line reaches a new site point.
 *     - Circle events: When the sweep-line becomes tangent to a circle defined by three
 *       sites whose corresponding arcs are adjacent on the "beach line". This signals the
 *       creation of a Voronoi vertex.
 * 3.  The "beach line" is a sequence of parabolic arcs. It is the locus of points
 *     equidistant from a site and the sweep-line. This implementation uses a simple
 *     doubly-linked list for the beach line, making it an O(n^2) algorithm.
 * 4.  The output is a set of Voronoi vertices. The final result for the benchmark is the
 *     sum of the coordinates of these vertices.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>

#define BOX_SIZE 10000.0
#define EPSILON 1e-9

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

// --- Data Structures ---
typedef struct {
    double x, y;
} Point;

typedef struct Site {
    Point p;
} Site;

struct Arc;
struct Edge;

typedef struct Event {
    double y;           // Event y-coordinate (sweep-line position)
    Point p;            // Point associated with event (site or Voronoi vertex)
    struct Arc *arc;    // Arc to be removed for a circle event
    int is_circle_event; // Flag: 1 for circle event, 0 for site event
    int is_valid;       // Flag to invalidate circle events
    Site* site;         // Site that causes the event
} Event;

typedef struct Arc {
    Point site_p;
    struct Arc *prev, *next;
    Event *circle_event;
} Arc;

// --- Global Benchmark Data ---
static int num_sites;
static Site* sites = NULL;
static double final_result = 0.0;

// Priority queue (min-heap) for events
static Event** event_queue = NULL;
static int event_queue_size = 0;
static int event_queue_capacity = 0;

// Doubly linked list for the beach line
static Arc* beach_line_root = NULL;

// Result data
static Point* voronoi_vertices = NULL;
static int num_vertices = 0;
static int vertices_capacity = 0;

// --- Helper Functions Prototypes ---
static void heap_push(Event* e);
static Event* heap_pop();
static double get_parabola_intersection_x(Point p1, Point p2, double sweep_y);
static void check_for_circle_event(Arc* arc);
static void handle_site_event(Event* e);
static void handle_circle_event(Event* e);
static void add_vertex(Point p);

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_sites> <seed>\n", argv[0]);
        exit(1);
    }
    num_sites = atoi(argv[1]);
    uint32_t seed = atoi(argv[2]);
    mt_seed(seed);

    sites = (Site*)malloc(num_sites * sizeof(Site));
    if (!sites) { fprintf(stderr, "Failed to allocate sites\n"); exit(1); }

    // Initial capacity for event queue should be at least num_sites
    event_queue_capacity = num_sites * 2;
    event_queue = (Event**)malloc(event_queue_capacity * sizeof(Event*));
    if (!event_queue) { fprintf(stderr, "Failed to allocate event queue\n"); exit(1); }

    // Voronoi vertices can be up to 2n-5
    vertices_capacity = num_sites * 2;
    voronoi_vertices = (Point*)malloc(vertices_capacity * sizeof(Point));
    if (!voronoi_vertices) { fprintf(stderr, "Failed to allocate vertices\n"); exit(1); }

    for (int i = 0; i < num_sites; i++) {
        sites[i].p.x = (double)mt_rand() / UINT32_MAX * BOX_SIZE;
        sites[i].p.y = (double)mt_rand() / UINT32_MAX * BOX_SIZE;
        
        Event* e = (Event*)malloc(sizeof(Event));
        e->y = sites[i].p.y;
        e->p = sites[i].p;
        e->site = &sites[i];
        e->is_circle_event = 0;
        e->is_valid = 1;
        e->arc = NULL;
        heap_push(e);
    }
}

void run_computation() {
    while (event_queue_size > 0) {
        Event* e = heap_pop();
        if (!e->is_valid) {
            free(e);
            continue;
        }

        if (e->is_circle_event) {
            handle_circle_event(e);
        } else {
            handle_site_event(e);
        }
        free(e);
    }

    // Accumulate result
    final_result = 0.0;
    for(int i = 0; i < num_vertices; i++) {
        final_result += voronoi_vertices[i].x + voronoi_vertices[i].y;
    }
}

void cleanup() {
    free(sites);
    for (int i = 0; i < event_queue_size; i++) {
        free(event_queue[i]);
    }
    free(event_queue);

    Arc* arc = beach_line_root;
    while(arc != NULL) {
        Arc* next = arc->next;
        free(arc);
        arc = next;
    }

    free(voronoi_vertices);
}

// --- Main and Timing ---
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

// --- Helper Function Implementations ---

static void add_vertex(Point p) {
    if (num_vertices >= vertices_capacity) {
        vertices_capacity *= 2;
        voronoi_vertices = (Point*)realloc(voronoi_vertices, vertices_capacity * sizeof(Point));
    }
    voronoi_vertices[num_vertices++] = p;
}

static void heap_push(Event* e) {
    if (event_queue_size >= event_queue_capacity) {
        event_queue_capacity *= 2;
        event_queue = (Event**)realloc(event_queue, event_queue_capacity * sizeof(Event*));
    }
    event_queue[event_queue_size] = e;
    int i = event_queue_size++;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (event_queue[p]->y < event_queue[i]->y) break;
        Event *tmp = event_queue[p];
        event_queue[p] = event_queue[i];
        event_queue[i] = tmp;
        i = p;
    }
}

static Event* heap_pop() {
    if (event_queue_size == 0) return NULL;

    Event* top = event_queue[0];
    event_queue[0] = event_queue[--event_queue_size];
    int i = 0;
    while (1) {
        int l = 2 * i + 1;
        int r = 2 * i + 2;
        int smallest = i;
        if (l < event_queue_size && event_queue[l]->y < event_queue[smallest]->y) smallest = l;
        if (r < event_queue_size && event_queue[r]->y < event_queue[smallest]->y) smallest = r;
        if (smallest == i) break;
        Event *tmp = event_queue[i];
        event_queue[i] = event_queue[smallest];
        event_queue[smallest] = tmp;
        i = smallest;
    }
    return top;
}

static double get_parabola_intersection_x(Point p1, Point p2, double sweep_y) {
    double x1 = p1.x, y1 = p1.y;
    double x2 = p2.x, y2 = p2.y;
    double dx = x1 - x2;
    double dy = y1 - y2;

    if (fabs(dy) < EPSILON) {
        return (x1 + x2) / 2.0;
    }

    double a = 1.0 / (2.0 * (y1 - sweep_y)) - 1.0 / (2.0 * (y2 - sweep_y));
    double b = -2.0 * (x1 / (2.0 * (y1 - sweep_y)) - x2 / (2.0 * (y2 - sweep_y)));
    double c = (x1 * x1 + y1 * y1 - sweep_y * sweep_y) / (2.0 * (y1 - sweep_y)) - (x2 * x2 + y2 * y2 - sweep_y * sweep_y) / (2.0 * (y2 - sweep_y));

    // simplified since y^2 term of parabola cancels
    double p1_denom = 2 * (p1.y - sweep_y);
    double p2_denom = 2 * (p2.y - sweep_y);

    double a_simple = 1/p1_denom - 1/p2_denom;
    double b_simple = -2*(p1.x/p1_denom - p2.x/p2_denom);
    double c_simple = (p1.x*p1.x + p1.y*p1.y - sweep_y*sweep_y)/p1_denom 
                    - (p2.x*p2.x + p2.y*p2.y - sweep_y*sweep_y)/p2_denom;
    
    if (fabs(a_simple) < EPSILON) { // Parabolas don't intersect, happens for collinear sites horizontal to sweep line 
        return (p1.x + p2.x)/2.0; 
    } 
    double delta = b_simple*b_simple - 4*a_simple*c_simple;
    return (-b_simple - sqrt(delta))/(2*a_simple);
}


static void check_for_circle_event(Arc *arc) {
    if (!arc || !arc->prev || !arc->next) return;

    Point p1 = arc->prev->site_p;
    Point p2 = arc->site_p;
    Point p3 = arc->next->site_p;

    // check for collinear points, which don't form a circle
    if ((p2.y - p1.y) * (p3.x - p2.x) - (p3.y - p2.y) * (p2.x - p1.x) > -EPSILON) {
        return;
    }

    double D = 2 * (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y));
    if (fabs(D) < EPSILON) return; // Collinear

    Point center;
    center.x = ((p1.x * p1.x + p1.y * p1.y) * (p2.y - p3.y) + (p2.x * p2.x + p2.y * p2.y) * (p3.y - p1.y) + (p3.x * p3.x + p3.y * p3.y) * (p1.y - p2.y)) / D;
    center.y = ((p1.x * p1.x + p1.y * p1.y) * (p3.x - p2.x) + (p2.x * p2.x + p2.y * p2.y) * (p1.x - p3.x) + (p3.x * p3.x + p3.y * p3.y) * (p2.x - p1.x)) / D;
    
    double radius = sqrt((p1.x - center.x) * (p1.x - center.x) + (p1.y - center.y) * (p1.y - center.y));
    double event_y = center.y - radius;

    Event* e = (Event*)malloc(sizeof(Event));
    e->y = event_y;
    e->p = center;
    e->arc = arc;
    e->is_circle_event = 1;
    e->is_valid = 1;

    arc->circle_event = e;
    heap_push(e);
}

static void handle_site_event(Event* e_site) {
    Point site_p = e_site->site->p;
    if (!beach_line_root) {
        beach_line_root = (Arc*)malloc(sizeof(Arc));
        beach_line_root->site_p = site_p;
        beach_line_root->prev = beach_line_root->next = NULL;
        beach_line_root->circle_event = NULL;
        return;
    }

    // Find arc on beach line directly above the new site
    Arc* arc_above = beach_line_root;
    while (arc_above->next) {
        double intersection_x = get_parabola_intersection_x(arc_above->site_p, arc_above->next->site_p, site_p.y);
        if (site_p.x < intersection_x) break;
        arc_above = arc_above->next;
    }

    // Invalidate the circle event for the arc we are breaking
    if (arc_above->circle_event) {
        arc_above->circle_event->is_valid = 0;
    }

    // Break arc_above and insert new arc
    Arc* new_arc = (Arc*)malloc(sizeof(Arc));
    new_arc->site_p = site_p;
    new_arc->circle_event = NULL;

    Arc* arc_right = (Arc*)malloc(sizeof(Arc));
    arc_right->site_p = arc_above->site_p;
    arc_right->circle_event = NULL;

    new_arc->prev = arc_above;
    new_arc->next = arc_right;
    
    arc_right->prev = new_arc;
    arc_right->next = arc_above->next;
    if(arc_above->next) arc_above->next->prev = arc_right;
    arc_above->next = new_arc;

    // Check for new circle events
    check_for_circle_event(arc_above);
    check_for_circle_event(arc_right);
}

static void handle_circle_event(Event* e) {
    Arc* arc = e->arc;
    add_vertex(e->p);
    
    if (arc->prev) arc->prev->next = arc->next;
    if (arc->next) arc->next->prev = arc->prev;
    
    if (arc == beach_line_root) {
        beach_line_root = arc->next ? arc->next : arc->prev;
    }
    
    if (arc->prev && arc->prev->circle_event) arc->prev->circle_event->is_valid = 0;
    if (arc->next && arc->next->circle_event) arc->next->circle_event->is_valid = 0;

    check_for_circle_event(arc->prev);
    check_for_circle_event(arc->next);
    
    free(arc);
}
