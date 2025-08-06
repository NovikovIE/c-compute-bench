#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

// --- Benchmark Globals ---
#define ALPHABET_SIZE 26

typedef struct {
    int len;
    int link;
    int next[ALPHABET_SIZE];
} State;

char* text;
int text_length;
State* sa;
int max_states;
long long final_result;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();

// --- Benchmark Setup ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <text_length> <seed>\n", argv[0]);
        exit(1);
    }
    text_length = atoi(argv[1]);
    uint32_t seed = (uint32_t)atoi(argv[2]);

    if(text_length <= 0) {
        fprintf(stderr, "FATAL: text_length must be positive.\n");
        exit(1);
    }

    mt_seed(seed);

    text = (char*)malloc((text_length + 1) * sizeof(char));
    if (!text) {
        fprintf(stderr, "FATAL: Failed to allocate memory for text.\n");
        exit(1);
    }

    for (int i = 0; i < text_length; ++i) {
        text[i] = (char)('a' + (mt_rand() % ALPHABET_SIZE));
    }
    text[text_length] = '\0';

    // Max states for a string of length N is 2*N-1 (for N>1)
    max_states = 2 * text_length;
    sa = (State*)malloc(max_states * sizeof(State));
    if (!sa) {
        fprintf(stderr, "FATAL: Failed to allocate memory for suffix automaton.\n");
        free(text);
        exit(1);
    }
}

// --- Core Computation ---
void run_computation() {
    int sz = 1;
    int last = 0;
    
    // Initialize the root state (state 0)
    sa[0].len = 0;
    sa[0].link = -1;
    memset(sa[0].next, -1, sizeof(sa[0].next));

    for (int i = 0; i < text_length; ++i) {
        int c = text[i] - 'a';
        int cur = sz++;
        sa[cur].len = sa[last].len + 1;
        memset(sa[cur].next, -1, sizeof(sa[cur].next));

        int p = last;
        while (p != -1 && sa[p].next[c] == -1) {
            sa[p].next[c] = cur;
            p = sa[p].link;
        }

        if (p == -1) {
            sa[cur].link = 0;
        } else {
            int q = sa[p].next[c];
            if (sa[q].len == sa[p].len + 1) {
                sa[cur].link = q;
            } else {
                int clone = sz++;
                sa[clone].len = sa[p].len + 1;
                memcpy(sa[clone].next, sa[q].next, sizeof(sa[q].next));
                sa[clone].link = sa[q].link;

                while (p != -1 && sa[p].next[c] == q) {
                    sa[p].next[c] = clone;
                    p = sa[p].link;
                }
                sa[q].link = clone;
                sa[cur].link = clone;
            }
        }
        last = cur;
    }

    long long checksum = 0;
    for (int i = 0; i < sz; i++) {
        checksum += sa[i].len + sa[i].link;
        for (int j = 0; j < ALPHABET_SIZE; ++j) {
            if (sa[i].next[j] != -1) {
                checksum += sa[i].next[j];
            }
        }
    }
    final_result = checksum;
}

// --- Cleanup ---
void cleanup() {
    free(text);
    free(sa);
}

// --- Main ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%lld\n", final_result);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
