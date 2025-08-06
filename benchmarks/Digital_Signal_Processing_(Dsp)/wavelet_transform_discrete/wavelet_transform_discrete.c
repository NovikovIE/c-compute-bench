/* 
 * Copyright 2024 The Computer Benchmarks Organization.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

/************************************************************************
 *                 MERSENNE TWISTER (DO NOT MODIFY)                     *
 ************************************************************************/
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

/************************************************************************
 *                        BENCHMARK SPECIFIC CODE                       *
 ************************************************************************/

// Benchmark parameters
static int signal_length;
static int decomposition_levels;

// Data structures
static double* signal_data;
static double* workspace;

// Result accumulator
static double final_result = 0.0;

// --- BENCHMARK SETUP ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <signal_length> <decomposition_levels> <seed>\n", argv[0]);
        exit(1);
    }

    signal_length = atoi(argv[1]);
    decomposition_levels = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    // Input validation
    if (signal_length <= 0 || (signal_length & (signal_length - 1)) != 0) {
        fprintf(stderr, "ERROR: signal_length must be a positive power of 2.\n");
        exit(1);
    }
    if (decomposition_levels <= 0 || signal_length < (1L << decomposition_levels)) {
        fprintf(stderr, "ERROR: decomposition_levels is invalid or too large for the given signal_length.\n");
        exit(1);
    }

    // Seed the random number generator
    mt_seed(seed);

    // Allocate memory for the signal and a workspace buffer
    signal_data = (double*)malloc(signal_length * sizeof(double));
    workspace = (double*)malloc(signal_length * sizeof(double));
    if (!signal_data || !workspace) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize the signal with random values between 0.0 and 1.0
    for (int i = 0; i < signal_length; i++) {
        signal_data[i] = (double)mt_rand() / (double)UINT32_MAX;
    }
}

// --- BENCHMARK COMPUTATION ---
void run_computation() {
    // This function performs a 1D Discrete Wavelet Transform (DWT)
    // using the Haar wavelet for a specified number of decomposition levels.
    // The transform separates a signal into approximation (low-frequency)
    // and detail (high-frequency) coefficients.
    
    int current_length = signal_length;
    double norm = 1.0 / sqrt(2.0);
    
    for (int level = 0; level < decomposition_levels; ++level) {
        int half_length = current_length / 2;
        for (int i = 0; i < half_length; ++i) {
            // Approximation coefficients (low-pass filter)
            workspace[i] = (signal_data[2 * i] + signal_data[2 * i + 1]) * norm;
            // Detail coefficients (high-pass filter)
            workspace[i + half_length] = (signal_data[2 * i] - signal_data[2 * i + 1]) * norm;
        }

        // Copy the transformed data (approximation and detail coefficients) 
        // back into the primary signal array for the next level.
        for (int i = 0; i < current_length; ++i) {
            signal_data[i] = workspace[i];
        }

        // The next decomposition operates on the approximation coefficients,
        // which are now in the first half of the array.
        current_length = half_length;
    }

    // Accumulate the final transformed signal to prevent dead code elimination.
    // The final array contains coefficients from all levels of decomposition.
    double sum = 0.0;
    for (int i = 0; i < signal_length; ++i) {
        sum += signal_data[i];
    }
    final_result = sum;
}

// --- BENCHMARK CLEANUP ---
void cleanup() {
    free(signal_data);
    free(workspace);
}

// --- MAIN FUNCTION ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", final_result);

    // Print total time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
