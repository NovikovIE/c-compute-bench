#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

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

// Data structure for each subject in the survival analysis
typedef struct {
    int time_to_event;
    int event_observed; // 1 if event occurred, 0 if censored
} Subject;

// --- Global Benchmark Data ---
int NUM_SUBJECTS;
Subject *subjects = NULL;
double final_result = 0.0;

// Comparison function for qsort to sort subjects by time.
// Primary sort key: time_to_event (ascending)
// Secondary sort key: event_observed (descending, to process events before censorings at the same time)
int compare_subjects(const void *a, const void *b) {
    const Subject *subj_a = (const Subject*)a;
    const Subject *subj_b = (const Subject*)b;
    if (subj_a->time_to_event != subj_b->time_to_event) {
        return (subj_a->time_to_event > subj_b->time_to_event) - (subj_a->time_to_event < subj_b->time_to_event);
    }
    return subj_b->event_observed - subj_a->event_observed;
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_subjects> <num_events> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_SUBJECTS = atoi(argv[1]);
    int target_events = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    if (NUM_SUBJECTS <= 0 || target_events < 0 || target_events > NUM_SUBJECTS) {
        fprintf(stderr, "FATAL: Invalid arguments. num_subjects > 0 and 0 <= num_events <= num_subjects\n");
        exit(1);
    }

    mt_seed(seed);

    subjects = (Subject*)malloc(NUM_SUBJECTS * sizeof(Subject));
    if (subjects == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for subjects.\n");
        exit(1);
    }

    // Create an array specifying event status for each subject, then shuffle it.
    // This ensures we have exactly `target_events`.
    int *event_statuses = (int*)malloc(NUM_SUBJECTS * sizeof(int));
    if (event_statuses == NULL) {
        fprintf(stderr, "FATAL: Memory allocation failed for event_statuses.\n");
        free(subjects);
        exit(1);
    }
    for (int i = 0; i < NUM_SUBJECTS; ++i) {
        event_statuses[i] = (i < target_events) ? 1 : 0;
    }

    // Shuffle event statuses using Fisher-Yates algorithm
    for (int i = NUM_SUBJECTS - 1; i > 0; --i) {
        uint32_t j = mt_rand() % (i + 1);
        int temp = event_statuses[i];
        event_statuses[i] = event_statuses[j];
        event_statuses[j] = temp;
    }

    // Assign data to subjects. Max time is arbitrary, chosen for variety.
    const int max_time = 5000;
    for (int i = 0; i < NUM_SUBJECTS; ++i) {
        subjects[i].time_to_event = mt_rand() % max_time + 1;
        subjects[i].event_observed = event_statuses[i];
    }

    free(event_statuses);
}

void run_computation() {
    // The primary computational work is sorting the large dataset.
    qsort(subjects, NUM_SUBJECTS, sizeof(Subject), compare_subjects);

    double survival_prob = 1.0;
    double result_accumulator = 0.0;
    int subjects_at_risk = NUM_SUBJECTS;

    // Iterate through the sorted subjects to calculate the Kaplan-Meier estimator.
    int i = 0;
    while (i < NUM_SUBJECTS) {
        int current_time = subjects[i].time_to_event;
        int events_at_this_time = 0;
        int processed_at_this_time = 0;

        // Group subjects at the current time point
        int j = i;
        while (j < NUM_SUBJECTS && subjects[j].time_to_event == current_time) {
            if (subjects[j].event_observed) {
                events_at_this_time++;
            }
            processed_at_this_time++;
            j++;
        }

        // Apply Kaplan-Meier formula if there were events and subjects at risk.
        if (events_at_this_time > 0 && subjects_at_risk > 0) {
             survival_prob *= (1.0 - ((double)events_at_this_time / (double)subjects_at_risk));
        }

        // Accumulate result to prevent dead-code elimination.
        // Summing the survival probability at each distinct time step is a reasonable proxy.
        result_accumulator += survival_prob;

        // Update the number of subjects at risk for the next time point.
        subjects_at_risk -= processed_at_this_time;

        // Advance index to the next distinct time point.
        i = j;
    }

    final_result = result_accumulator;
}

void cleanup() {
    if (subjects != NULL) {
        free(subjects);
        subjects = NULL;
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

    // Print the final accumulated result to stdout
    printf("%.6f\n", final_result);

    // Print the execution time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
