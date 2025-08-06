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
// --- END MERSENNE TWISTER ---

// --- Global Benchmark Parameters ---
int batch_size;
int sequence_length;
int embedding_dim;
int num_heads;
int d_head;
int ffn_hidden_dim;

// --- Global Data Structures ---
// Inputs
float *input_tensor;

// Model Weights
float *w_q, *w_k, *w_v, *w_o;
float *b_q, *b_k, *b_v, *b_o;
float *w_ffn1, *b_ffn1;
float *w_ffn2, *b_ffn2;
float *ln1_gamma, *ln1_beta;
float *ln2_gamma, *ln2_beta;

// Intermediate Buffers
float *q_proj, *k_proj, *v_proj;
float *attn_scores;
float *attn_output;
float *attn_heads_concat;
float *add_norm_1_output;
float *ffn1_output;
float *add_norm_2_input;

// Final Output
float *output_tensor;
double final_result = 0.0;

// --- Helper Functions ---
float random_float() {
    return ((float)mt_rand() / (float)UINT32_MAX) * 2.0f - 1.0f;
}

void matmul_add_bias(float *out, const float *inp, const float *weight, const float *bias, int M, int N, int K) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            float sum = 0.0f;
            for (int k = 0; k < K; k++) {
                sum += inp[i * K + k] * weight[k * N + j];
            }
            out[i * N + j] = sum + bias[j];
        }
    }
}

void elementwise_add(float *out, const float *a, const float *b, long long n) {
    for (long long i = 0; i < n; i++) {
        out[i] = a[i] + b[i];
    }
}

void layer_norm(float *out, const float *inp, const float *gamma, const float *beta, int C) {
    const float epsilon = 1e-5f;
    float mean = 0.0f;
    for (int i = 0; i < C; i++) {
        mean += inp[i];
    }
    mean /= C;

    float variance = 0.0f;
    for (int i = 0; i < C; i++) {
        float diff = inp[i] - mean;
        variance += diff * diff;
    }
    variance /= C;

    float inv_std = 1.0f / sqrtf(variance + epsilon);
    for (int i = 0; i < C; i++) {
        out[i] = gamma[i] * (inp[i] - mean) * inv_std + beta[i];
    }
}

void softmax(float *x, int size) {
    float max_val = x[0];
    for (int i = 1; i < size; i++) {
        if (x[i] > max_val) max_val = x[i];
    }
    float sum_exp = 0.0f;
    for (int i = 0; i < size; i++) {
        x[i] = expf(x[i] - max_val);
        sum_exp += x[i];
    }
    for (int i = 0; i < size; i++) {
        x[i] /= sum_exp;
    }
}

void relu(float *x, long long n) {
    for (long long i = 0; i < n; i++) {
        x[i] = fmaxf(0.0f, x[i]);
    }
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 6) {
        fprintf(stderr, "Usage: %s batch_size sequence_length embedding_dim num_heads seed\n", argv[0]);
        exit(1);
    }

    batch_size = atoi(argv[1]);
    sequence_length = atoi(argv[2]);
    embedding_dim = atoi(argv[3]);
    num_heads = atoi(argv[4]);
    uint32_t seed = (uint32_t)atoi(argv[5]);
    mt_seed(seed);

    if (embedding_dim % num_heads != 0) {
        fprintf(stderr, "embedding_dim must be divisible by num_heads\n");
        exit(1);
    }
    d_head = embedding_dim / num_heads;
    ffn_hidden_dim = embedding_dim * 4;

    // Allocate memory
    long long batch_seq_embed = (long long)batch_size * sequence_length * embedding_dim;
    long long batch_seq_ffn = (long long)batch_size * sequence_length * ffn_hidden_dim;
    long long embed_embed = (long long)embedding_dim * embedding_dim;
    long long embed_ffn = (long long)embedding_dim * ffn_hidden_dim;
    long long ffn_embed = (long long)ffn_hidden_dim * embedding_dim;
    long long batch_head_seq_seq = (long long)batch_size * num_heads * sequence_length * sequence_length;
    
    input_tensor = (float*)malloc(batch_seq_embed * sizeof(float));
    w_q = (float*)malloc(embed_embed * sizeof(float));
    w_k = (float*)malloc(embed_embed * sizeof(float));
    w_v = (float*)malloc(embed_embed * sizeof(float));
    w_o = (float*)malloc(embed_embed * sizeof(float));
    b_q = (float*)malloc(embedding_dim * sizeof(float));
    b_k = (float*)malloc(embedding_dim * sizeof(float));
    b_v = (float*)malloc(embedding_dim * sizeof(float));
    b_o = (float*)malloc(embedding_dim * sizeof(float));
    w_ffn1 = (float*)malloc(embed_ffn * sizeof(float));
    b_ffn1 = (float*)malloc(ffn_hidden_dim * sizeof(float));
    w_ffn2 = (float*)malloc(ffn_embed * sizeof(float));
    b_ffn2 = (float*)malloc(embedding_dim * sizeof(float));
    ln1_gamma = (float*)malloc(embedding_dim * sizeof(float));
    ln1_beta = (float*)malloc(embedding_dim * sizeof(float));
    ln2_gamma = (float*)malloc(embedding_dim * sizeof(float));
    ln2_beta = (float*)malloc(embedding_dim * sizeof(float));

    q_proj = (float*)malloc(batch_seq_embed * sizeof(float));
    k_proj = (float*)malloc(batch_seq_embed * sizeof(float));
    v_proj = (float*)malloc(batch_seq_embed * sizeof(float));
    attn_scores = (float*)malloc(batch_head_seq_seq * sizeof(float));
    attn_output = (float*)malloc(batch_seq_embed * sizeof(float));
    attn_heads_concat = (float*)malloc(batch_seq_embed * sizeof(float));
    add_norm_1_output = (float*)malloc(batch_seq_embed * sizeof(float));
    ffn1_output = (float*)malloc(batch_seq_ffn * sizeof(float));
    output_tensor = (float*)malloc(batch_seq_embed * sizeof(float));

    // Initialize with random data
    for(long long i = 0; i < batch_seq_embed; i++) input_tensor[i] = random_float();
    for(long long i = 0; i < embed_embed; i++) w_q[i] = random_float();
    for(long long i = 0; i < embed_embed; i++) w_k[i] = random_float();
    for(long long i = 0; i < embed_embed; i++) w_v[i] = random_float();
    for(long long i = 0; i < embed_embed; i++) w_o[i] = random_float();
    for(int i = 0; i < embedding_dim; i++) { b_q[i] = random_float(); b_k[i] = random_float(); b_v[i] = random_float(); b_o[i] = random_float(); }
    for(long long i = 0; i < embed_ffn; i++) w_ffn1[i] = random_float();
    for(int i = 0; i < ffn_hidden_dim; i++) b_ffn1[i] = random_float();
    for(long long i = 0; i < ffn_embed; i++) w_ffn2[i] = random_float();
    for(int i = 0; i < embedding_dim; i++) b_ffn2[i] = random_float();
    for(int i = 0; i < embedding_dim; i++) { ln1_gamma[i] = 1.0f; ln1_beta[i] = 0.0f; ln2_gamma[i] = 1.0f; ln2_beta[i] = 0.0f; }
}

void run_computation() {
    float scale = 1.0f / sqrtf((float)d_head);
    long long seq_embed_stride = sequence_length * embedding_dim;

    // Project Q, K, V for all batches at once
    matmul_add_bias(q_proj, input_tensor, w_q, b_q, batch_size * sequence_length, embedding_dim, embedding_dim);
    matmul_add_bias(k_proj, input_tensor, w_k, b_k, batch_size * sequence_length, embedding_dim, embedding_dim);
    matmul_add_bias(v_proj, input_tensor, w_v, b_v, batch_size * sequence_length, embedding_dim, embedding_dim);

    for (int b = 0; b < batch_size; b++) {
        for (int h = 0; h < num_heads; h++) {
            // Pointers to the current batch's data for Q, K, V
            const float *q_b = q_proj + b * seq_embed_stride;
            const float *k_b = k_proj + b * seq_embed_stride;
            const float *v_b = v_proj + b * seq_embed_stride;
            
            float *current_head_scores = attn_scores + (b * num_heads + h) * sequence_length * sequence_length;

            // 1. Scaled Dot-Product Attention: (Q*K^T)/sqrt(d_k)
            for (int i = 0; i < sequence_length; i++) {
                for (int j = 0; j < sequence_length; j++) {
                    float score = 0.0f;
                    for (int k = 0; k < d_head; k++) {
                        score += q_b[i * embedding_dim + h * d_head + k] * k_b[j * embedding_dim + h * d_head + k];
                    }
                    current_head_scores[i * sequence_length + j] = score * scale;
                }
            }
            
            // 2. Softmax
            for (int i = 0; i < sequence_length; i++) {
                softmax(current_head_scores + i * sequence_length, sequence_length);
            }

            // 3. Matmul with V
            float *current_head_concat = attn_heads_concat + b * seq_embed_stride;
            for (int i = 0; i < sequence_length; i++) {
                for (int k = 0; k < d_head; k++) {
                    float value = 0.0f;
                    for (int j = 0; j < sequence_length; j++) {
                        value += current_head_scores[i * sequence_length + j] * v_b[j * embedding_dim + h * d_head + k];
                    }
                    current_head_concat[i * embedding_dim + h * d_head + k] = value;
                }
            }
        }
    }

    // 4. Final attention projection
    matmul_add_bias(attn_output, attn_heads_concat, w_o, b_o, batch_size * sequence_length, embedding_dim, embedding_dim);

    // 5. Residual Connection + LayerNorm_1
    float *residual_1_input = (float*) malloc(seq_embed_stride * sizeof(float));
    for (int b = 0; b < batch_size; b++) {
        const float* current_input = input_tensor + b * seq_embed_stride;
        const float* current_attn_out = attn_output + b * seq_embed_stride;
        float* current_addnorm1_out = add_norm_1_output + b * seq_embed_stride;
        elementwise_add(residual_1_input, current_input, current_attn_out, seq_embed_stride);
        for (int i=0; i < sequence_length; i++) {
            layer_norm(current_addnorm1_out + i * embedding_dim, residual_1_input + i * embedding_dim, ln1_gamma, ln1_beta, embedding_dim);
        }
    }
    free(residual_1_input);

    // 6. Feed-Forward Network: Layer 1 + ReLU
    matmul_add_bias(ffn1_output, add_norm_1_output, w_ffn1, b_ffn1, batch_size * sequence_length, ffn_hidden_dim, embedding_dim);
    relu(ffn1_output, (long long)batch_size * sequence_length * ffn_hidden_dim);

    // 7. Feed-Forward Network: Layer 2
    // This will go into `output_tensor` temporarily
    matmul_add_bias(output_tensor, ffn1_output, w_ffn2, b_ffn2, batch_size * sequence_length, embedding_dim, ffn_hidden_dim);

    // 8. Residual Connection + LayerNorm_2
    float *residual_2_input = (float*) malloc(seq_embed_stride * sizeof(float));
    for (int b = 0; b < batch_size; b++) {
        const float* current_addnorm1_out = add_norm_1_output + b * seq_embed_stride;
        const float* current_ffn_out = output_tensor + b * seq_embed_stride;
        float* final_current_out = output_tensor + b * seq_embed_stride;
        elementwise_add(residual_2_input, current_addnorm1_out, current_ffn_out, seq_embed_stride);
        for (int i=0; i < sequence_length; i++) {
             layer_norm(final_current_out + i * embedding_dim, residual_2_input + i * embedding_dim, ln2_gamma, ln2_beta, embedding_dim);
        }
    }
    free(residual_2_input);

    // Accumulate final result to prevent dead code elimination
    long long total_elements = (long long)batch_size * sequence_length * embedding_dim;
    for (long long i = 0; i < total_elements; i++) {
        final_result += output_tensor[i];
    }
}

void cleanup() {
    free(input_tensor);
    free(w_q); free(w_k); free(w_v); free(w_o);
    free(b_q); free(b_k); free(b_v); free(b_o);
    free(w_ffn1); free(b_ffn1);
    free(w_ffn2); free(b_ffn2);
    free(ln1_gamma); free(ln1_beta);
    free(ln2_gamma); free(ln2_beta);
    free(q_proj); free(k_proj); free(v_proj);
    free(attn_scores);
    free(attn_output);
    free(attn_heads_concat);
    free(add_norm_1_output);
    free(ffn1_output);
    free(output_tensor);
}

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