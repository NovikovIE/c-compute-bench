#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- MERSENNE TWISTER (DO NOT MODIFY) ---
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

// Benchmark-specific globals
static int M, N;          // Matrix dimensions M >= N
static double *A;         // Input matrix M x N
static double *A_copy;    // Copy of A for in-place computation
static double *S;         // Singular values (vector of size N)
static double *V;         // Right singular vectors (N x N matrix)
static double final_result;

// Helper macros and functions for SVD
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double max(double a, double b) {
    return a > b ? a : b;
}

// Numerically stable computation of sqrt(a^2 + b^2)
static double pythag(double a, double b) {
    double absa = fabs(a);
    double absb = fabs(b);
    if (absa > absb) {
        return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
    } else {
        return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb)));
    }
}

/**
 * @brief Singular Value Decomposition using Golub-Kahan-Reinsch algorithm.
 * 
 * This implementation is adapted from public domain versions inspired by the
 * algorithm described in "Numerical Recipes in C". It computes the SVD of a 
 * matrix A = U*S*V^T. It assumes M >= N.
 * @param a M-by-N matrix. On output, `a` is replaced by the matrix U.
 * @param m Number of rows of matrix `a`.
 * @param n Number of columns of matrix `a`.
 * @param s Output vector of N singular values.
 * @param v Output N-by-N orthogonal matrix V.
 */
void perform_svd(double *a, int m, int n, double *s, double *v) {
    int i, j, k, l = 0, its;
    double anorm, c, f, g, h, scale, x, y, z;
    
    double *rv1 = (double *)malloc(n * sizeof(double));
    if (!rv1) { fprintf(stderr, "FATAL: SVD temporary allocation failed\n"); exit(1); }
    
    // Householder reduction to bidiagonal form
    g = scale = anorm = 0.0;
    for (i = 0; i < n; i++) {
        l = i + 1;
        rv1[i] = scale * g;
        g = s[i] = scale = 0.0;
        if (i < m) {
            for (k = i; k < m; k++) scale += fabs(a[k * n + i]);
            if (scale != 0.0) {
                double sum_sq = 0.0;
                for (k = i; k < m; k++) { a[k * n + i] /= scale; sum_sq += a[k * n + i] * a[k * n + i]; }
                f = a[i * n + i];
                g = -SIGN(sqrt(sum_sq), f);
                h = f * g - sum_sq;
                a[i * n + i] = f - g;
                for (j = l; j < n; j++) {
                    double t_sum = 0.0;
                    for (k = i; k < m; k++) t_sum += a[k * n + i] * a[k * n + j];
                    f = t_sum / h;
                    for (k = i; k < m; k++) a[k * n + j] += f * a[k * n + i];
                }
                for (k = i; k < m; k++) a[k * n + i] *= scale;
            }
        }
        s[i] = scale * g;
        g = scale = 0.0;
        if (i < m && i != n - 1) {
            for (k = l; k < n; k++) scale += fabs(a[i * n + k]);
            if (scale != 0.0) {
                double sum_sq = 0.0;
                for (k = l; k < n; k++) { a[i * n + k] /= scale; sum_sq += a[i * n + k] * a[i * n + k]; }
                f = a[i * n + l];
                g = -SIGN(sqrt(sum_sq), f);
                h = f * g - sum_sq;
                a[i * n + l] = f - g;
                for (k = l; k < n; k++) rv1[k] = a[i * n + k] / h;
                for (j = l; j < m; j++) {
                    double t_sum = 0.0;
                    for (k = l; k < n; k++) t_sum += a[j * n + k] * a[i * n + k];
                    for (k = l; k < n; k++) a[j * n + k] += t_sum * rv1[k];
                }
                for (k = l; k < n; k++) a[i * n + k] *= scale;
            }
        }
        anorm = max(anorm, (fabs(s[i]) + fabs(rv1[i])));
    }

    // Accumulation of right-hand transformations
    for (i = n - 1; i >= 0; i--) {
        if (i < n - 1) {
            if (g != 0.0) {
                for (j = l; j < n; j++) v[j * n + i] = (a[i * n + j] / a[i * n + l]) / g;
                for (j = l; j < n; j++) {
                    double t_sum = 0.0;
                    for (k = l; k < n; k++) t_sum += a[i * n + k] * v[k * n + j];
                    for (k = l; k < n; k++) v[k * n + j] += t_sum * v[k * n + i];
                }
            }
            for (j = l; j < n; j++) v[i * n + j] = v[j * n + i] = 0.0;
        }
        v[i * n + i] = 1.0;
        g = rv1[i];
        l = i;
    }

    // Diagonalization by QR iteration
    for (k = n - 1; k >= 0; k--) {
        for (its = 0; its < 30; its++) {
            int flag = 1;
            for (l = k; l >= 0; l--) {
                if (l==0 || fabs(rv1[l]) < 1e-12 * anorm) { flag = 0; break; }
                if (fabs(s[l - 1]) < 1e-12 * anorm) break;
            }
            if (flag) {
                c = 0.0; f = 1.0;
                for (i = l; i <= k; i++) {
                    g = rv1[i]; y = s[i]; h = f * g;
                    if (fabs(h) > 1e-12 * anorm) {
                         g = y; z = pythag(f, h); s[i] = z;
                         c = f / z; f = y / z; rv1[i] = h / z;
                         for (j = 0; j < n; j++) { 
                             x = v[j * n + l - 1]; y = v[j * n + i];
                             v[j * n + l - 1] = x * c + y * f; v[j * n + i] = -x * f + y * c; 
                         } 
                    } else { rv1[i] = 0.0; }
                }
            }
            z = s[k];
            if (l == k) {
                if (z < 0.0) { s[k] = -z; for (j = 0; j < n; j++) v[j * n + k] = -v[j * n + k]; }
                break;
            }
            if (its >= 29) { /* Not a fatal error */ }

            x = s[l]; y = s[k - 1]; g = rv1[k - 1]; h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * (y / (f + SIGN(g, f)) - h)) / x;
            c = f = 1.0;
            for (j = l; j < k; j++) {
                i = j + 1; g = rv1[i]; y = s[i];
                h = f * g; g = c * g; z = pythag(f, h);
                rv1[j] = z; c = f / z; f = h / z;
                x = x * c + g * f; g = -x * f + g * c; h = y * f; y = y * c;
                for (int jj = 0; jj < n; jj++) {
                     z = v[jj * n + j]; v[jj * n + j] = z * c + v[jj * n + i] * f;
                     v[jj * n + i] = -z * f + v[jj * n + i] * c;
                }
                z = pythag(x, h); s[j] = z;
                if(z) { c = x / z; f = h / z; }
                x = c * g + f * y; y = -f * g + c * y;
            }
            rv1[l] = 0.0; rv1[k] = x; s[k] = y;
        }
    }
    free(rv1);
}

void setup_benchmark(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <matrix_height> <matrix_width> <seed>\n", argv[0]);
        exit(1);
    }

    M = atoi(argv[1]);
    N = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    if (M <= 0 || N <= 0) {
        fprintf(stderr, "FATAL: Matrix dimensions must be positive.\n");
        exit(1);
    }
    
    if (M < N) {
        fprintf(stderr, "FATAL: This implementation requires matrix_height >= matrix_width.\n");
        exit(1);
    }

    mt_seed(seed);

    A = (double *)malloc((size_t)M * N * sizeof(double));
    A_copy = (double *)malloc((size_t)M * N * sizeof(double));
    S = (double *)malloc((size_t)N * sizeof(double));
    V = (double *)malloc((size_t)N * N * sizeof(double));

    if (!A || !A_copy || !S || !V) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = ((double)mt_rand() / (double)(UINT32_MAX / 2.0)) - 1.0;
        }
    }
}

void run_computation() {
    memcpy(A_copy, A, (size_t)M * N * sizeof(double));
    
    // The matrix 'A_copy' is destroyed by the SVD routine.
    // It returns singular values in 'S' and the V matrix in 'V'.
    // The modified 'A_copy' matrix becomes U.
    perform_svd(A_copy, M, N, S, V);

    final_result = 0.0;
    for (int i = 0; i < N; i++) {
        final_result += S[i];
    }
}

void cleanup() {
    free(A);
    free(A_copy);
    free(S);
    free(V);
}

int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout to prevent dead code elimination
    printf("%f\n", final_result);

    // Print timing to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
