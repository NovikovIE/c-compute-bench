# c-compute-bench

**c-compute-bench** is a curated suite of diverse, CPU-bound computational benchmarks written in C. It is designed for researchers, compiler engineers, and hardware architects to reliably measure and compare the performance of different CPU architectures and compiler optimization flags.

What makes this suite unique is its **generation and refinement methodology**. The entire codebase was created, debugged, and performance-tuned through an advanced, automated pipeline powered by Large Language Models (LLMs), ensuring both code quality and measurement accuracy.

Every benchmark in this collection adheres to a strict structural requirement:
1.  **`setup_benchmark()`**: A dedicated function that allocates memory and generates all necessary input data.
2.  **`run_computation()`**: A function that executes **only** the core computational algorithm on the pre-generated data.
3.  **`cleanup()`**: A function that frees all allocated memory.

A native C timing mechanism, built into every program, measures **only the pure execution time of the `run_computation()` phase.** This guarantees that you are benchmarking the algorithm itself, not the overhead of its setup.

## Key Features of the Suite

*   **Diverse Workloads:** Includes dozens of benchmarks across various domains like Numerical Methods, Computer Graphics, Graph Algorithms, Cryptography, and more.
*   **High-Quality Measurements:** By isolating the computation phase, the results are more accurate and reliable for comparing subtle differences between compiler flags or CPU architectures. Random number generation is strictly controlled, using only the Mersenne Twister algorithm to ensure reproducibility on different platforms.
*   **Time-Targeted:** Every benchmark has been iteratively tuned so that its pure computation time is approximately **1 second** on a reference modern CPU. This provides a substantial and consistent workload for stable measurements.
*   **Rigorously Validated:** The entire suite has been tested with AddressSanitizer (ASan) and UndefinedBehaviorSanitizer (UBSan) to eliminate common memory errors and undefined behavior, ensuring code correctness.
*   **Self-Contained and Portable:** Each benchmark is a single, self-contained C file with no external data dependencies, making it easy to compile and run anywhere with GCC or Clang.
*   **Transparent Generation:** The LLM-powered scripts used to generate, refine, and test this suite are included in the repository for full transparency and potential extension.

## Getting Started: Using the Benchmarks

### 1. Prerequisites

*   A C compiler: `gcc` (v11+) or `clang`.
*   **For Linux (Debian/Ubuntu):** `sudo apt update && sudo apt install build-essential`

### 2. How to Run

The suite includes a powerful execution script (`run.py`) that handles compilation and testing for you.

#### Example 1: Basic Performance Measure

```bash
python run.py
```

This will compile and run each benchmark and save a detailed report to `benchmark_resulsts.csv`. The key column to analyze is `time_seconds`.

#### Example 2: Running a Single Benchmark

To test only the `matrix_multiplication` program:

```bash
python run.py -p matrix_multiplication
```

#### Example 3: Verifying Code Quality with Sanitizers

This project was already validated, but you can re-verify it on your own system. The script uses `clang` on macOS for best results.

```bash
# Check for memory errors (e.g., out-of-bounds access)
python run.py --sanitizer address

# Check for undefined behavior (e.g., integer overflow)
python run.py --sanitizer undefined
```
This mode will report any errors found by the sanitizer.

### 3. Understanding the Results

After running `run.py`, you will get a `.csv` file with the following important columns:

*   `program`: The name of the benchmark.
*   `cflags`: The compiler flags used for that run.
*   `status`: `SUCCESS`, `CRASHED`, or `FAILED (SANITIZER)`.
*   `time_seconds`: **The pure computation time in seconds.** This is the primary metric for performance analysis.
*   `reason`: A detailed error message if the status is not `SUCCESS`.

## For Developers: The Generation Pipeline

While this repository is presented as a ready-to-use benchmark suite, the tools that built it are included. This allows you to extend the suite with new benchmarks or modify the existing ones.

*   `generate_programs.py`: Generates the initial C code from high-level descriptions.
*   `tune_time.py`: The iterative time-targeting and auto-tuning script.
*   `autofix.py`: Script to auto-fix buggy code.

scripts are written for Yandex inner API, so to use this scripts you need to change logic to use public gemini API or other any LLM API.

## Directory Structure

```
.
├── benchmarks/         # The benchmark suite source code.
├── build/              # Directory for compiled executables (created on run).
├── run.py          # The primary script for running and testing the benchmarks.
├── generate_programs.py
├── tune_time.py
├── autofix.py
├── all_benchmarks.py
└── README.md           # This file
```

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.