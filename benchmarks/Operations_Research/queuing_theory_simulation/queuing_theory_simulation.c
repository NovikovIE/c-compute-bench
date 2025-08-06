#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- BEGIN MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK SPECIFIC CODE ---

#define ARRIVAL 0
#define DEPARTURE 1

// Simulation parameters
static int NUM_SERVERS;
static double ARRIVAL_RATE;
static double SERVICE_RATE;
static long NUM_EVENTS;

// Data structures for the simulation
typedef struct {
    double time;
    int type;
} Event;

// Global state
static Event* event_heap; // Min-heap for event list
static long heap_size;
static long heap_capacity;

static double* waiting_queue; // FIFO queue for customer arrival times
static long queue_head, queue_tail, queue_count, queue_capacity;

static double total_wait_time = 0.0;

// --- Helper Functions ---

double exponential_rand(double rate) {
    double u;
    do {
        u = (double)mt_rand() / (double)UINT32_MAX;
    } while (u == 1.0); // Avoid log(0)
    return -log(1.0 - u) / rate;
}

void _swap_events(long i, long j) {
    Event temp = event_heap[i];
    event_heap[i] = event_heap[j];
    event_heap[j] = temp;
}

void _sift_up(long index) {
    if (index == 0) return;
    long parent_index = (index - 1) / 2;
    if (event_heap[index].time < event_heap[parent_index].time) {
        _swap_events(index, parent_index);
        _sift_up(parent_index);
    }
}

void heap_push(Event event) {
    if (heap_size >= heap_capacity) return; // Should not happen with proper capacity setting
    event_heap[heap_size] = event;
    heap_size++;
    _sift_up(heap_size - 1);
}

void _sift_down(long index) {
    long left_child_idx = 2 * index + 1;
    long right_child_idx = 2 * index + 2;
    long smallest = index;

    if (left_child_idx < heap_size && event_heap[left_child_idx].time < event_heap[smallest].time) {
        smallest = left_child_idx;
    }
    if (right_child_idx < heap_size && event_heap[right_child_idx].time < event_heap[smallest].time) {
        smallest = right_child_idx;
    }
    if (smallest != index) {
        _swap_events(index, smallest);
        _sift_down(smallest);
    }
}

Event heap_pop() {
    Event popped_event = event_heap[0];
    heap_size--;
    event_heap[0] = event_heap[heap_size];
    _sift_down(0);
    return popped_event;
}

void queue_push(double arrival_time) {
    if (queue_count >= queue_capacity) return; // Should not happen
    waiting_queue[queue_tail] = arrival_time;
    queue_tail = (queue_tail + 1) % queue_capacity;
    queue_count++;
}

double queue_pop() {
    if (queue_count == 0) return -1.0; // Error
    double arrival_time = waiting_queue[queue_head];
    queue_head = (queue_head + 1) % queue_capacity;
    queue_count--;
    return arrival_time;
}

// --- Core Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s num_servers arrival_rate service_rate num_events seed\n", argv[0]);
        exit(1);
    }

    NUM_SERVERS = atoi(argv[1]);
    ARRIVAL_RATE = atof(argv[2]);
    SERVICE_RATE = atof(argv[3]);
    NUM_EVENTS = atol(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);

    mt_seed(seed);

    heap_capacity = NUM_EVENTS + 2; // Should be enough for pending events
    event_heap = (Event*)malloc(heap_capacity * sizeof(Event));
    if (!event_heap) {
        perror("Failed to allocate event heap");
        exit(1);
    }
    heap_size = 0;

    queue_capacity = NUM_EVENTS + 1;
    waiting_queue = (double*)malloc(queue_capacity * sizeof(double));
    if (!waiting_queue) {
        perror("Failed to allocate waiting queue");
        exit(1);
    }
    queue_head = 0;
    queue_tail = 0;
    queue_count = 0;

    // Start with a single arrival event
    Event first_event = {exponential_rand(ARRIVAL_RATE), ARRIVAL};
    heap_push(first_event);
}

void run_computation() {
    double current_time = 0.0;
    long events_processed = 0;
    int busy_servers = 0;

    while (events_processed < NUM_EVENTS) {
        if (heap_size == 0) break;

        Event current_event = heap_pop();
        current_time = current_event.time;

        if (current_event.type == ARRIVAL) {
            // Schedule the next arrival
            double next_arrival_time = current_time + exponential_rand(ARRIVAL_RATE);
            Event next_arrival_event = {next_arrival_time, ARRIVAL};
            heap_push(next_arrival_event);

            // Process the current arrival
            if (busy_servers < NUM_SERVERS) {
                busy_servers++;
                double service_time = exponential_rand(SERVICE_RATE);
                Event departure_event = {current_time + service_time, DEPARTURE};
                heap_push(departure_event);
            } else {
                queue_push(current_time);
            }
        } else { // DEPARTURE
            if (queue_count > 0) {
                double arrival_time = queue_pop();
                total_wait_time += current_time - arrival_time;

                // This server immediately starts serving the next person in the queue
                double service_time = exponential_rand(SERVICE_RATE);
                Event departure_event = {current_time + service_time, DEPARTURE};
                heap_push(departure_event);
            } else {
                busy_servers--;
            }
        }
        events_processed++;
    }
}

void cleanup() {
    free(event_heap);
    free(waiting_queue);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", total_wait_time);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
