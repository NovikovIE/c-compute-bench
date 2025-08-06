#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- START MERSENNE TWISTER (DO NOT MODIFY) ---
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

// --- BENCHMARK DATA & PARAMETERS ---
typedef struct {
    int num_samples;      // n
    int num_predictors;   // p
    int mcmc_iterations;

    // Input data
    double *X; // n x p matrix (flattened)
    double *y; // n-element vector

    // Precomputed values from input data
    double *XTX; // p x p matrix
    double *XTy; // p-element vector
    double yTy; // y' * y scalar

    // MCMC chain storage
    double *beta_samples;   // iters x p matrix
    double *sigma2_samples; // iters-element vector

    // Final single result to prevent dead-code elimination
    double final_result;
} BenchmarkData;

static BenchmarkData g_data;

// --- UTILITY FUNCTIONS ---

// Generate a random double in [0, 1)
double rand_double() {
    return (double)mt_rand() / (double)UINT32_MAX;
}

// Box-Muller transform to generate two standard normal random variables
void generate_normal_pair(double *z0, double *z1) {
    double u1, u2;
    do {
        u1 = rand_double();
        u2 = rand_double();
    } while (u1 <= 1e-9); // Avoid log(0)
    *z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    *z1 = sqrt(-2.0 * log(u1)) * sin(2.0 * M_PI * u2);
}

// Marsaglia's squeeze method for Gamma(shape, 1) where shape > 1
double generate_gamma(double shape) {
    if (shape <= 0.0) return NAN;
    if (shape == 1.0) return -log(1.0 - rand_double()); // Exponentialverteilung
    if (shape < 1.0) {
        // Uses algorithm from Marsaglia & Tsang (2000) for shape < 1
        double u, v, x;
        u = rand_double();
        v = generate_gamma(shape + 1.0);
        return pow(u, 1.0/shape) * v;
    }

    // Marsaglia's Squeeze method for shape > 1
    double d = shape - 1.0 / 3.0;
    double c = 1.0 / sqrt(9.0 * d);
    double x, v, u;
    while (1) {
        do {
            generate_normal_pair(&x, &v); // Reuse v as the second normal
            v = 1.0 + c * x;
        } while (v <= 0.0);
        v = v * v * v;
        u = rand_double();
        if (u < 1.0 - 0.0331 * x * x * x * x) return (d * v);
        if (log(u) < 0.5 * x * x + d * (1.0 - v + log(v))) return (d * v);
    }
}

// Cholesky decomposition of a symmetric positive-definite matrix A (p x p)
// Produces lower-triangular matrix L such that A = L * L^T
int cholesky(const double *A, double *L, int p) {
    memset(L, 0, p * p * sizeof(double));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < j; k++) {
                sum += L[i * p + k] * L[j * p + k];
            }
            if (i == j) {
                double d = A[i * p + i] - sum;
                if (d <= 0.0) return -1; // Not positive-definite
                L[i * p + j] = sqrt(d);
            } else {
                L[i * p + j] = (A[i * p + j] - sum) / L[j * p + j];
            }
        }
    }
    return 0;
}

// Solves L*x = b for x, where L is a lower-triangular matrix
void forward_subst(const double *L, const double *b, double *x, int p) {
    for (int i = 0; i < p; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i * p + j] * x[j];
        }
        x[i] = (b[i] - sum) / L[i * p + i];
    }
}

// Solves L^T*x = b for x, where L is a lower-triangular matrix
void backward_subst(const double *L, const double *b, double *x, int p) {
    for (int i = p - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < p; j++) {
            sum += L[j * p + i] * x[j];
        }
        x[i] = (b[i] - sum) / L[i * p + i];
    }
}


// --- BENCHMARK FUNCTIONS ---

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s num_samples num_predictors mcmc_iterations seed\n", argv[0]);
        exit(1);
    }

    g_data.num_samples = atoi(argv[1]);
    g_data.num_predictors = atoi(argv[2]);
    g_data.mcmc_iterations = atoi(argv[3]);
    uint32_t seed = (uint32_t)atoi(argv[4]);
    mt_seed(seed);

    int n = g_data.num_samples;
    int p = g_data.num_predictors;

    // Allocate memory
    g_data.X = (double*)malloc(n * p * sizeof(double));
    g_data.y = (double*)malloc(n * sizeof(double));
    g_data.beta_samples = (double*)malloc(g_data.mcmc_iterations * p * sizeof(double));
    g_data.sigma2_samples = (double*)malloc(g_data.mcmc_iterations * sizeof(double));
    g_data.XTX = (double*)malloc(p * p * sizeof(double));
    g_data.XTy = (double*)malloc(p * sizeof(double));
    double *beta_true = (double*)malloc(p * sizeof(double));

    // Generate synthetic data
    for (int i = 0; i < p; i++) {
        beta_true[i] = (rand_double() - 0.5) * 2.0; // Uniform(-1, 1)
    }
    double true_sigma = 2.0;

    for (int i = 0; i < n; i++) {
        double y_hat = 0.0;
        for (int j = 0; j < p; j++) {
            g_data.X[i * p + j] = rand_double() * 2.0; // Uniform(0, 2)
            y_hat += g_data.X[i * p + j] * beta_true[j];
        }
        double z0, z1;
        generate_normal_pair(&z0, &z1);
        g_data.y[i] = y_hat + true_sigma * z0;
    }

    free(beta_true);

    // Precompute X'X, X'y, and y'y
    for (int i = 0; i < p; i++) {
        for (int j = i; j < p; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += g_data.X[k * p + i] * g_data.X[k * p + j];
            }
            g_data.XTX[i * p + j] = sum;
            g_data.XTX[j * p + i] = sum; // Symmetric
        }
        double sum_y = 0.0;
        for (int k = 0; k < n; k++) {
            sum_y += g_data.X[k * p + i] * g_data.y[k];
        }
        g_data.XTy[i] = sum_y;
    }

    g_data.yTy = 0.0;
    for (int i = 0; i < n; i++) {
        g_data.yTy += g_data.y[i] * g_data.y[i];
    }
}

void run_computation() {
    int n = g_data.num_samples;
    int p = g_data.num_predictors;
    int iters = g_data.mcmc_iterations;

    // Temporary arrays for computation
    double *beta_current = (double*)calloc(p, sizeof(double));
    double sigma2_current = 1.0;
    double *precision_matrix = (double*)malloc(p * p * sizeof(double));
    double *cholesky_L = (double*)malloc(p * p * sizeof(double));
    double *mean_beta = (double*)malloc(p * sizeof(double));
    double *b_vec = (double*)malloc(p * sizeof(double));
    double *z_vec = (double*)malloc(p * sizeof(double));
    double *tmp_vec = (double*)malloc(p * sizeof(double));

    // Gibbs sampling loop
    for (int i = 0; i < iters; i++) {
        // 1. Sample beta | sigma^2, y, X
        // Precision = XTX/sigma^2 + I
        double inv_sigma2 = 1.0 / sigma2_current;
        for (int r = 0; r < p; r++) {
            for (int c = 0; c < p; c++) {
                precision_matrix[r * p + c] = g_data.XTX[r * p + c] * inv_sigma2;
            }
            precision_matrix[r * p + r] += 1.0; // Add identity
        }

        cholesky(precision_matrix, cholesky_L, p);
        
        // Solve for mean_beta: Precision * mean_beta = XTy/sigma^2
        for(int j=0; j<p; ++j) b_vec[j] = g_data.XTy[j] * inv_sigma2;
        forward_subst(cholesky_L, b_vec, tmp_vec, p); // tmp_vec = L^-1 * b_vec
        backward_subst(cholesky_L, tmp_vec, mean_beta, p); // mean_beta = (L^T)^-1 * tmp_vec

        // Sample from N(mean_beta, Precision^-1)
        for(int j=0; j<p; j+=2) {
            if (j + 1 < p) {
                generate_normal_pair(&z_vec[j], &z_vec[j+1]);
            } else {
                double z0, z1; // z1 is discarded
                generate_normal_pair(&z0, &z1);
                z_vec[j] = z0; 
            }
        }
        backward_subst(cholesky_L, z_vec, tmp_vec, p); // tmp_vec = (L^T)^-1 * z
        for(int j=0; j<p; ++j) {
            beta_current[j] = mean_beta[j] + tmp_vec[j];
            g_data.beta_samples[i * p + j] = beta_current[j];
        }

        // 2. Sample sigma^2 | beta, y, X
        // Sum of squared errors: sse = y'y - 2y'Xb + b'X'Xb
        double yTx_b = 0.0;
        for(int j=0; j<p; ++j) yTx_b += g_data.XTy[j] * beta_current[j];

        double b_XTX_b = 0.0;
        for(int j=0; j<p; ++j) {
            tmp_vec[j] = 0.0;
            for(int k=0; k<p; ++k) {
                 tmp_vec[j] += g_data.XTX[j * p + k] * beta_current[k];
            }
        }
        for(int j=0; j<p; ++j) b_XTX_b += beta_current[j] * tmp_vec[j];

        double sse = g_data.yTy - 2.0 * yTx_b + b_XTX_b;
        
        // Draw from IG(n/2, sse/2)
        double shape_a = (double)n / 2.0;
        double scale_b = sse / 2.0;
        sigma2_current = scale_b / generate_gamma(shape_a);
        g_data.sigma2_samples[i] = sigma2_current;
    }

    // Accumulate results to prevent dead code elimination
    g_data.final_result = 0.0;
    for(int i = 0; i < iters; i++){
        for(int j = 0; j < p; j++){
            g_data.final_result += g_data.beta_samples[i * p + j];
        }
        g_data.final_result += g_data.sigma2_samples[i];
    }

    // Free temporary memory
    free(beta_current); free(precision_matrix); free(cholesky_L); free(mean_beta); 
    free(b_vec); free(z_vec); free(tmp_vec);
}

void cleanup() {
    free(g_data.X);
    free(g_data.y);
    free(g_data.beta_samples);
    free(g_data.sigma2_samples);
    free(g_data.XTX);
    free(g_data.XTy);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    printf("%f\n", g_data.final_result);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
