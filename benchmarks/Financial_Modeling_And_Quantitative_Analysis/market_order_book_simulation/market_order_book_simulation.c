#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

// --- BEGIN: Mersenne Twister (DO NOT MODIFY) ---
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
// --- END: Mersenne Twister ---

// --- Benchmark Data Structures ---
typedef enum {
    LIMIT_BUY, LIMIT_SELL, MARKET_BUY, MARKET_SELL
} OrderType;

typedef struct {
    double price;
    uint32_t quantity;
} PriceLevel;

typedef struct {
    OrderType type;
    double price;     // Used for LIMIT orders
    uint32_t quantity;
} MarketEvent;

// --- Global Variables ---
long long NUM_MARKET_EVENTS;
int ORDER_BOOK_DEPTH;

MarketEvent* market_events; // Array of events to process
PriceLevel* bid_book;       // Sorted array (descending price)
PriceLevel* ask_book;       // Sorted array (ascending price)

double total_volume_traded = 0.0; // The final result

// --- Function Declarations ---
void setup_benchmark(int argc, char* argv[]);
void run_computation();
void cleanup();
static void insert_order(PriceLevel* book, double price, uint32_t quantity, int is_bid);

// --- Benchmark Implementation ---
void setup_benchmark(int argc, char* argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <num_market_events> <order_book_depth> <seed>\n", argv[0]);
        exit(1);
    }

    NUM_MARKET_EVENTS = atoll(argv[1]);
    ORDER_BOOK_DEPTH = atoi(argv[2]);
    uint32_t seed = atoi(argv[3]);

    mt_seed(seed);

    market_events = (MarketEvent*) malloc(NUM_MARKET_EVENTS * sizeof(MarketEvent));
    bid_book = (PriceLevel*) calloc(ORDER_BOOK_DEPTH, sizeof(PriceLevel));
    ask_book = (PriceLevel*) calloc(ORDER_BOOK_DEPTH, sizeof(PriceLevel));

    if (!market_events || !bid_book || !ask_book) {
        fprintf(stderr, "FATAL: Memory allocation failed.\n");
        exit(1);
    }

    // Initialize the order book with a plausible spread
    double center_price = 100.0;
    double tick_size = 0.01;
    for (int i = 0; i < ORDER_BOOK_DEPTH; ++i) {
        bid_book[i].price = center_price - (i * tick_size);
        bid_book[i].quantity = (mt_rand() % 1000) + 100; // 100 to 1099 shares
        ask_book[i].price = center_price + tick_size + (i * tick_size);
        ask_book[i].quantity = (mt_rand() % 1000) + 100;
    }

    // Generate market events
    for (long long i = 0; i < NUM_MARKET_EVENTS; ++i) {
        market_events[i].type = (OrderType)(mt_rand() % 4);
        market_events[i].quantity = (mt_rand() % 100) + 1; // 1 to 100 shares

        if (market_events[i].type == LIMIT_BUY || market_events[i].type == LIMIT_SELL) {
            double price_offset = ((double)(mt_rand() % 200) - 100.0) * tick_size; // +/- 1.00 from center
             market_events[i].price = center_price + price_offset;
             market_events[i].price = (int)(market_events[i].price * 100.0 + 0.5) / 100.0; // round to 2 decimal places
        } else {
            market_events[i].price = 0; // Not used for market orders
        }
    }
}

void cleanup() {
    free(market_events);
    free(bid_book);
    free(ask_book);
}

// Inserts an order into a sorted book. is_bid=1 for bid book, 0 for ask book.
static void insert_order(PriceLevel* book, double price, uint32_t quantity, int is_bid) {
    int insert_pos = -1;
    for (int i = 0; i < ORDER_BOOK_DEPTH; ++i) {
        if (book[i].price == price) { // Found matching price level
            book[i].quantity += quantity;
            return;
        }
        if (is_bid ? (price > book[i].price) : (price < book[i].price)) {
            if (book[i].price == 0 && insert_pos == -1) { // Empty slot at end
                insert_pos = i;
            } else if (insert_pos == -1) { 
                insert_pos = i;
            }
        } else if (book[i].price == 0 && insert_pos == -1) { // Empty slot
             insert_pos = i;
        }
    }

    if (insert_pos != -1) {
        // Shift elements to make space for the new order
        memmove(&book[insert_pos + 1], &book[insert_pos], (ORDER_BOOK_DEPTH - 1 - insert_pos) * sizeof(PriceLevel));
        book[insert_pos].price = price;
        book[insert_pos].quantity = quantity;
    }
    // If no position found, order is not good enough to make the book (dropped)
}

void run_computation() {
    for (long long i = 0; i < NUM_MARKET_EVENTS; ++i) {
        MarketEvent* event = &market_events[i];
        uint32_t remaining_qty = event->quantity;

        switch (event->type) {
            case LIMIT_BUY:
            case MARKET_BUY: {
                // Match against ask book
                while (remaining_qty > 0 && ask_book[0].price > 0 && 
                       (event->type == MARKET_BUY || event->price >= ask_book[0].price)) {

                    uint32_t trade_qty = (remaining_qty < ask_book[0].quantity) ? remaining_qty : ask_book[0].quantity;
                    total_volume_traded += trade_qty * ask_book[0].price;
                    remaining_qty -= trade_qty;
                    ask_book[0].quantity -= trade_qty;

                    if (ask_book[0].quantity == 0) {
                        // Level cleared, shift book up
                        memmove(&ask_book[0], &ask_book[1], (ORDER_BOOK_DEPTH - 1) * sizeof(PriceLevel));
                        ask_book[ORDER_BOOK_DEPTH - 1].price = 0;
                        ask_book[ORDER_BOOK_DEPTH - 1].quantity = 0;
                    }
                }
                // If it was a limit order with quantity left, add to bid book
                if (event->type == LIMIT_BUY && remaining_qty > 0) {
                    insert_order(bid_book, event->price, remaining_qty, 1);
                }
                break;
            }

            case LIMIT_SELL:
            case MARKET_SELL: {
                 // Match against bid book
                while (remaining_qty > 0 && bid_book[0].price > 0 &&
                       (event->type == MARKET_SELL || event->price <= bid_book[0].price)) {

                    uint32_t trade_qty = (remaining_qty < bid_book[0].quantity) ? remaining_qty : bid_book[0].quantity;
                    total_volume_traded += trade_qty * bid_book[0].price;
                    remaining_qty -= trade_qty;
                    bid_book[0].quantity -= trade_qty;

                    if (bid_book[0].quantity == 0) {
                        // Level cleared, shift book up
                        memmove(&bid_book[0], &bid_book[1], (ORDER_BOOK_DEPTH - 1) * sizeof(PriceLevel));
                        bid_book[ORDER_BOOK_DEPTH - 1].price = 0;
                        bid_book[ORDER_BOOK_DEPTH - 1].quantity = 0;
                    }
                }
                // If it was a limit order with quantity left, add to ask book
                if (event->type == LIMIT_SELL && remaining_qty > 0) {
                    insert_order(ask_book, event->price, remaining_qty, 0);
                }
                break;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    struct timespec start, end;

    setup_benchmark(argc, argv);

    clock_gettime(CLOCK_MONOTONIC, &start);
    run_computation();
    clock_gettime(CLOCK_MONOTONIC, &end);

    cleanup();

    double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    // Print result to stdout
    fprintf(stdout, "%f\n", total_volume_traded);

    // Print time to stderr
    fprintf(stderr, "%.6f", time_taken);

    return 0;
}
