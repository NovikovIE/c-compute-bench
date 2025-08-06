#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>

// --- Mersenne Twister (MT19937) Generator (Verbatim) ---
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

// Parameters
int num_bonds;
int spline_complexity; // Represents the number of iterations for the yield solver

// Data structures
typedef struct {
    double maturity;     // Time to maturity in years
    double coupon_rate;  // Annual coupon rate
    double market_price; // Observed market price
    double face_value;   // Principal amount repaid at maturity
} Bond;

Bond* bonds;       // Array of bonds to be priced
double* yields;    // Array to store the computed yield-to-maturity for each bond

// Final result
double result_checksum = 0.0;

// --- Benchmark Functions ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_bonds> <spline_complexity> <seed>\n", argv[0]);
        exit(1);
    }

    num_bonds = atoi(argv[1]);
    spline_complexity = atoi(argv[2]);
    uint32_t seed = (uint32_t)atoi(argv[3]);

    mt_seed(seed);

    bonds = (Bond*)malloc(num_bonds * sizeof(Bond));
    yields = (double*)malloc(num_bonds * sizeof(double));

    if (!bonds || !yields) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < num_bonds; i++) {
        // Generate somewhat realistic bond data
        bonds[i].maturity = 1.0 + (double)i * (40.0 / num_bonds) + (mt_rand() / (double)UINT32_MAX - 0.5) * 0.2;
        if (bonds[i].maturity < 0.5) bonds[i].maturity = 0.5; // Ensure positive maturity
        bonds[i].coupon_rate = 0.01 + (mt_rand() / (double)UINT32_MAX) * 0.07; // 1% to 8% coupon
        bonds[i].face_value = 1000.0;
        // Price is a function of a simplified underlying yield curve plus noise
        double true_yield = 0.02 + 0.001 * bonds[i].maturity; 
        bonds[i].market_price = (bonds[i].coupon_rate / true_yield) * bonds[i].face_value * (1.0 - pow(1.0 + true_yield, -bonds[i].maturity)) + bonds[i].face_value * pow(1.0 + true_yield, -bonds[i].maturity);
        bonds[i].market_price *= (0.99 + (mt_rand() / (double)UINT32_MAX) * 0.02); // Add market noise +/- 1%

        // Initialize yield guess
        yields[i] = 0.05; // Initial guess of 5% for all bonds
    }
}

void run_computation() {
    const int MAX_GLOBAL_ITERATIONS = 10; // Number of times to refine the entire curve
    
    for (int iter = 0; iter < MAX_GLOBAL_ITERATIONS; ++iter) {
        for (int i = 0; i < num_bonds; i++) {
            double current_yield = yields[i];

            // Use Newton-Raphson method to solve for Yield-to-Maturity (YTM)
            // 'spline_complexity' controls the number of iterations (depth of search)
            for (int j = 0; j < spline_complexity; j++) {
                double price_model = 0.0;
                double price_derivative = 0.0; // aka Macaulay Duration * Price / (1+y)

                // Sum of present values of coupon payments (annual coupons for simplicity)
                for (int t = 1; t <= (int)floor(bonds[i].maturity); t++) {
                    double coupon_payment = bonds[i].coupon_rate * bonds[i].face_value;
                    price_model += coupon_payment * pow(1.0 + current_yield, -(double)t);
                    price_derivative += -((double)t) * coupon_payment * pow(1.0 + current_yield, -((double)t + 1.0));
                }

                // Present value of face value payment
                price_model += bonds[i].face_value * pow(1.0 + current_yield, -bonds[i].maturity);
                price_derivative += -bonds[i].maturity * bonds[i].face_value * pow(1.0 + current_yield, -(bonds[i].maturity + 1.0));

                double price_error = price_model - bonds[i].market_price;
                
                // Avoid division by zero
                if (fabs(price_derivative) < 1e-8) {
                    break;
                }

                // Newton-Raphson update step
                current_yield -= price_error / price_derivative;
            }
            yields[i] = current_yield;
        }
    }

    // Calculate a checksum to prevent dead code elimination and provide a result
    double checksum = 0.0;
    for (int i = 0; i < num_bonds; i++) {
        checksum += yields[i];
    }
    result_checksum = checksum;
}

void cleanup() {
    free(bonds);
    free(yields);
}

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print final result to stdout
    printf("%f\n", result_checksum);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
