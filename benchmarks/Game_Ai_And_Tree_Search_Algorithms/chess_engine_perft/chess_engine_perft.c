#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>

// --- Mersenne Twister (Do Not Modify - Include This Verbatim) ---
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
// --- End of Mersenne Twister ---

// --- Chess Engine Definitions ---
#define EMPTY 0
#define PAWN 1
#define KNIGHT 2
#define BISHOP 3
#define ROOK 4
#define QUEEN 5
#define KING 6

#define WHITE 8
#define BLACK 16

#define MAX_MOVES 256

typedef struct {
    int from;
    int to;
    int capture;
    int promotion;
} Move;

// --- Global Benchmark Data ---
static int g_depth;
static unsigned long long g_nodes;
static int g_board[64];
static int g_side_to_move;

// --- Function Prototypes ---
void setup_benchmark(int argc, char *argv[]);
void run_computation();
void cleanup();
unsigned long long perft(int depth);
void setup_initial_board();
void generate_moves(Move move_list[], int* move_count, int side);
void make_move(const Move* move);
void undo_move(const Move* move);

// --- Main Function ---
int main(int argc, char *argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("%llu\n", g_nodes);
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}

// --- Benchmark Functions ---
void setup_benchmark(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <depth> <seed>\n", argv[0]);
        exit(1);
    }

    g_depth = atoi(argv[1]);
    uint32_t seed = (uint32_t)strtoul(argv[2], NULL, 10);
    mt_seed(seed);

    setup_initial_board();
    g_nodes = 0;
}

void run_computation() {
    g_nodes = perft(g_depth);
}

void cleanup() {
    // No heap memory was allocated, so nothing to free.
}

// --- Chess Logic Implementations ---

void setup_initial_board() {
    const int initial_board[64] = {
        BLACK | ROOK,   BLACK | KNIGHT, BLACK | BISHOP, BLACK | QUEEN,  BLACK | KING,   BLACK | BISHOP, BLACK | KNIGHT, BLACK | ROOK,
        BLACK | PAWN,   BLACK | PAWN,   BLACK | PAWN,   BLACK | PAWN,   BLACK | PAWN,   BLACK | PAWN,   BLACK | PAWN,   BLACK | PAWN,
        EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,
        EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,
        EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,
        EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,          EMPTY,
        WHITE | PAWN,   WHITE | PAWN,   WHITE | PAWN,   WHITE | PAWN,   WHITE | PAWN,   WHITE | PAWN,   WHITE | PAWN,   WHITE | PAWN,
        WHITE | ROOK,   WHITE | KNIGHT, WHITE | BISHOP, WHITE | QUEEN,  WHITE | KING,   WHITE | BISHOP, WHITE | KNIGHT, WHITE | ROOK,
    };
    for (int i = 0; i < 64; ++i) {
        g_board[i] = initial_board[i];
    }
    g_side_to_move = WHITE;
}

void make_move(const Move* move) {
    g_board[move->to] = g_board[move->from];
    g_board[move->from] = EMPTY;
    if (move->promotion) {
        g_board[move->to] = (g_board[move->to] & (WHITE | BLACK)) | move->promotion;
    }
    g_side_to_move = (g_side_to_move == WHITE) ? BLACK : WHITE;
}

void undo_move(const Move* move) {
    g_side_to_move = (g_side_to_move == WHITE) ? BLACK : WHITE;
    int piece_type = (g_side_to_move | PAWN);
    if (move->promotion) {
        g_board[move->from] = piece_type;
    } else {
        g_board[move->from] = g_board[move->to];
    }
    g_board[move->to] = move->capture;
}

static void add_move(Move list[], int* count, int from, int to, int promotion) {
    list[*count].from = from;
    list[*count].to = to;
    list[*count].capture = g_board[to];
    list[*count].promotion = promotion;
    (*count)++;
}

/*
 * This is a pseudo-legal move generator. It generates all possible moves
 * for a piece without checking if the move leaves the king in check.
 * It also simplifies promotions to always be a queen and omits castling/en-passant.
 * This is sufficient for a computationally intensive perft benchmark.
 */
void generate_moves(Move move_list[], int* move_count, int side) {
    int opponent = (side == WHITE) ? BLACK : WHITE;

    for (int from = 0; from < 64; from++) {
        if ((g_board[from] & side) == 0) continue;

        int piece = g_board[from] & 7;
        int rank = from / 8;
        int file = from % 8;

        switch (piece) {
            case PAWN: {
                int dir = (side == WHITE) ? -8 : 8;
                int start_rank = (side == WHITE) ? 6 : 1;
                int promotion_rank = (side == WHITE) ? 0 : 7;

                // Single push
                if (g_board[from + dir] == EMPTY) {
                    if (rank == promotion_rank) add_move(move_list, move_count, from, from + dir, QUEEN);
                    else add_move(move_list, move_count, from, from + dir, 0);
                    // Double push
                    if (rank == start_rank && g_board[from + 2 * dir] == EMPTY) {
                        add_move(move_list, move_count, from, from + 2 * dir, 0);
                    }
                }
                // Captures
                if (file > 0 && (g_board[from + dir - 1] & opponent)) {
                     if (rank == promotion_rank) add_move(move_list, move_count, from, from + dir - 1, QUEEN);
                     else add_move(move_list, move_count, from, from + dir - 1, 0);
                }
                if (file < 7 && (g_board[from + dir + 1] & opponent)) {
                     if (rank == promotion_rank) add_move(move_list, move_count, from, from + dir + 1, QUEEN);
                     else add_move(move_list, move_count, from, from + dir + 1, 0);
                }
                break;
            }
            case KNIGHT: {
                const int offsets[] = {-17, -15, -10, -6, 6, 10, 15, 17};
                for (int i = 0; i < 8; ++i) {
                    int to = from + offsets[i];
                    if (to < 0 || to >= 64) continue;
                    int to_file = to % 8;
                    int from_file_dist = abs(to_file - file);
                    if (from_file_dist > 2) continue;
                    if ((g_board[to] & side) == 0) {
                        add_move(move_list, move_count, from, to, 0);
                    }
                }
                break;
            }
            case BISHOP: case ROOK: case QUEEN: {
                int start_dir = (piece == BISHOP) ? 4 : 0;
                int end_dir = (piece == ROOK) ? 4 : 8;
                const int directions[] = {-8, -1, 1, 8, -9, -7, 7, 9}; // N, W, E, S, NW, NE, SE, SW

                for (int i = start_dir; i < end_dir; ++i) {
                    int dir = directions[i];
                    int to = from;
                    while (1) {
                        to += dir;
                        if (to < 0 || to >= 64) break;
                        int to_file = to % 8;
                        int from_file = (to-dir) % 8;
                        if(abs(to_file - from_file) > 1 && (dir == -1 || dir==1 || dir == -9 || dir == -7 || dir == 7 || dir == 9)) break; // wrapped around

                        if (g_board[to] == EMPTY) {
                            add_move(move_list, move_count, from, to, 0);
                        } else {
                            if (g_board[to] & opponent) {
                                add_move(move_list, move_count, from, to, 0);
                            }
                            break;
                        }
                    }
                }
                break;
            }
            case KING: {
                const int offsets[] = {-9, -8, -7, -1, 1, 7, 8, 9};
                for (int i = 0; i < 8; ++i) {
                    int to = from + offsets[i];
                    if (to < 0 || to >= 64) continue;
                    int to_file = to % 8;
                    int from_file_dist = abs(to_file - file);
                    if (from_file_dist > 1) continue;

                    if ((g_board[to] & side) == 0) {
                        add_move(move_list, move_count, from, to, 0);
                    }
                }
                break;
            }
        }
    }
}

unsigned long long perft(int depth) {
    if (depth == 0) {
        return 1ULL;
    }

    Move move_list[MAX_MOVES];
    int move_count = 0;
    generate_moves(move_list, &move_count, g_side_to_move);

    if (depth == 1) {
        return (unsigned long long)move_count;
    }
    
    unsigned long long nodes = 0ULL;
    for (int i = 0; i < move_count; i++) {
        make_move(&move_list[i]);
        nodes += perft(depth - 1);
        undo_move(&move_list[i]);
    }

    return nodes;
}
