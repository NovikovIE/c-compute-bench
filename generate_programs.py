import os
import sys
import json
import requests
import time
import concurrent.futures
import textwrap

from all_benchmarks import all_benchmarks

from requests.packages.urllib3.exceptions import InsecureRequestWarning

requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

# --- Configuration ---
API_KEY = os.getenv("SOY_TOKEN")
API_URL = "https://api.eliza.yandex.net/openrouter/v1/chat/completions"
MODEL_NAME = "google/gemini-2.5-pro-preview-03-25"
OUTPUT_DIR = "benchmarks"
HEADERS = {
    "Authorization": f"OAuth {API_KEY}",
    "Content-Type": "application/json"
}

def create_prompt(theme_name, theme_desc, prog_name, params):
    """
    Generates a new, high-quality prompt that enforces separation of concerns
    and built-in timing for the computation phase.
    """
    c_params = params + ["seed"] if "seed" not in params else params
    param_list = "\n".join([f"    - `argv[{i+1}]`: {p}" for i, p in enumerate(c_params)])

    # ---> ГЛАВНОЕ ИЗМЕНЕНИЕ: ПОЛНОСТЬЮ НОВЫЙ ПРОМПТ <---
    prompt_text = textwrap.dedent(f"""
        You are an expert C programmer creating high-quality, self-contained CPU benchmarks.
        Your task is to generate a C program that strictly separates data setup from computation and includes precise internal timing.
        Your entire response must be a single, raw JSON object.

        **Benchmark Details:**
        - **Theme:** {theme_name}
        - **Description:** {theme_desc}
        - **Program Name:** {prog_name}
        - **Parameters:** {', '.join(params)}

        **Your Response MUST be a single JSON object:**
        ```json
        {{
          "c_code": "...",
          "args": "...",
          "args_values": "..."
        }}
        ```

        ---
        **C Program - General Requirements:**

        1.  **Structure:** The C code MUST be organized into the following three functions:
            - `void setup_benchmark(int argc, char *argv[])`: Parses `argv`, allocates memory, and generates all input data.
            - `void run_computation()`: Executes the core algorithm on the data prepared by `setup_benchmark`.
            - `void cleanup()`: Frees all memory allocated in `setup_benchmark`.
            These functions should use global pointers or a global struct to share data.

        2.  **Timing (CRITICAL):**
            - The `main` function must measure the execution time of **ONLY** the `run_computation()` call.
            - Use `clock_gettime(CLOCK_MONOTONIC, ...)` for timing. Do not time `setup_benchmark` or `cleanup`.
            - The full program MUST include `<time.h>`.

        3.  **Output (CRITICAL):**
            - The program must produce **two** distinct outputs:
              1. **`stdout`**: The final computational result (an accumulated value to prevent dead code elimination), followed by a newline. Example: `printf("%f\\n", result);`
              2. **`stderr`**: The measured time of `run_computation()` in seconds, as a floating-point number, with no newline. Example: `fprintf(stderr, "%.6f", time_in_seconds);`

        ---
        **C Program - Detailed Function Requirements:**

        **1. Data and Randomness:**
            - Use the provided Mersenne Twister (MT19937) generator for any random data.
            - The `seed` for the generator MUST be passed via `argv` and used in `setup_benchmark`.
            - Allocate large data structures on the heap (`malloc`).

        **2. `main` function:**
            - The `main` function should be very simple: call `setup_benchmark`, time the call to `run_computation`, call `cleanup`, and then print the outputs.
            - Example `main` structure:
              ```c
              int main(int argc, char *argv[]) {{
                  struct timespec start, end;
                  
                  setup_benchmark(argc, argv);
                  
                  clock_gettime(CLOCK_MONOTONIC, &start);
                  run_computation();
                  clock_gettime(CLOCK_MONOTONIC, &end);
                  
                  cleanup();
                  
                  double time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
                  
                  // Print result to stdout
                  printf("%d\\n", final_result); // Or %f, etc.
                  
                  // Print time to stderr
                  fprintf(stderr, "%.6f", time_taken);
                  
                  return 0;
              }}
              ```

        **3. Arguments:**
            - The program must accept command-line arguments as specified below:
        {param_list}
            - The last argument is ALWAYS the `seed`.

        <hr>
        
        **Mersenne Twister Generator (Do Not Modify - Include This Verbatim):**
        ```c
        #include <stdint.h>
        // stdio.h, stdlib.h, time.h should be included at the top of the file
        
        #define MT_N 624
        #define MT_M 397
        #define MT_MATRIX_A 0x9908b0dfUL
        #define MT_UPPER_MASK 0x80000000UL
        #define MT_LOWER_MASK 0x7fffffffUL
        
        static uint32_t mt[MT_N];
        static int mt_index = MT_N + 1;
        
        void mt_seed(uint32_t seed) {{
            mt[0] = seed;
            for (mt_index = 1; mt_index < MT_N; mt_index++) {{
                mt[mt_index] = (1812433253UL * (mt[mt_index - 1] ^ (mt[mt_index - 1] >> 30)) + mt_index);
            }}
        }}
        
        uint32_t mt_rand(void) {{
            uint32_t y;
            static const uint32_t mag01[2] = {{0x0UL, MT_MATRIX_A}};
            if (mt_index >= MT_N) {{
                if (mt_index > MT_N) {{
                     fprintf(stderr, "FATAL: Mersenne Twister not seeded.");
                     exit(1);
                }}
                for (int i = 0; i < MT_N - MT_M; i++) {{
                    y = (mt[i] & MT_UPPER_MASK) | (mt[i + 1] & MT_LOWER_MASK);
                    mt[i] = mt[i + MT_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
                }}
                for (int i = MT_N - MT_M; i < MT_N - 1; i++) {{
                    y = (mt[i] & MT_UPPER_MASK) | (mt[i + 1] & MT_LOWER_MASK);
                    mt[i] = mt[i + (MT_M - MT_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
                }}
                y = (mt[MT_N - 1] & MT_UPPER_MASK) | (mt[0] & MT_LOWER_MASK);
                mt[MT_N - 1] = mt[MT_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
                mt_index = 0;
            }}
            y = mt[mt_index++];
            y ^= (y >> 11); y ^= (y << 7) & 0x9d2c5680UL; y ^= (y << 15) & 0xefc60000UL; y ^= (y >> 18);
            return y;
        }}
        ```

        ---
        **Argument Requirements (`args` and `args_values`):**

        1.  **`args` (string):**
            - Comma-separated string of ALL parameter NAMES, including `seed`.
            - Example: `{",".join(params)},seed`

        2.  **`args_values` (string):**
            - Comma-separated string of VALUES.
            - **CRITICAL:** Choose values so that the **`run_computation()` function** takes **approximately 1 second** on a modern CPU.
            - Use a fixed integer like `12345` for the `seed`.
    """)
    return prompt_text


def make_api_request(prompt_content, prog_name):
    payload = {"model": MODEL_NAME, "messages": [{"role": "user", "content": prompt_content}]}
    generated_content_str = ""
    try:
        response = requests.post(API_URL, headers=HEADERS, json=payload, timeout=600, verify=False)
        response.raise_for_status()
        response_json = response.json()
        generated_content_str = response_json.get("response", {}).get("choices", [{}])[0].get("message", {}).get("content", "")
        if not generated_content_str:
            print(f"  [{prog_name}] [Error] Received empty content from API.")
            return None
        if generated_content_str.strip().startswith("```json"):
            cleaned_str = '\n'.join(generated_content_str.strip().split('\n')[1:-1])
        else:
            cleaned_str = generated_content_str
        return json.loads(cleaned_str)
    except requests.exceptions.RequestException as e:
        print(f"  [{prog_name}] [Error] Network or API error: {e}")
    except json.JSONDecodeError as e:
        print(f"  [{prog_name}] [Error] Failed to parse JSON: {e}\n[Raw Content]: {generated_content_str[:500]}...")
    except (KeyError, IndexError) as e:
        print(f"  [{prog_name}] [Error] Unexpected API response: {e}\n[Raw Response]: {response.text}")
    return None

def process_benchmark(theme, theme_desc, prog_name, params):
    prog_dir = os.path.join(OUTPUT_DIR, theme.replace(" ", "_").title(), prog_name)
    c_file_path = os.path.join(prog_dir, f"{prog_name}.c")
    args_file_path = os.path.join(prog_dir, "args.txt")
    args_values_file_path = os.path.join(prog_dir, "args_values.txt")

    if os.path.exists(c_file_path) and os.path.exists(args_file_path) and os.path.exists(args_values_file_path):
        print(f"  Skipping '{prog_name}' (already exists).")
        return (prog_name, "skipped")

    print(f"  Generating '{prog_name}'...")
    os.makedirs(prog_dir, exist_ok=True)
    prompt = create_prompt(theme, theme_desc, prog_name, params)
    result_json = make_api_request(prompt, prog_name)
    
    if result_json and "c_code" in result_json and "args" in result_json and "args_values" in result_json:
        if 'seed' not in result_json["args"].split(','):
            print(f"    -> Failed for '{prog_name}'. 'seed' is missing in returned 'args'.")
            return (prog_name, "failed")
        with open(c_file_path, "w", encoding="utf-8") as f:
            f.write(result_json["c_code"])
        with open(args_file_path, "w", encoding="utf-8") as f:
            f.write(result_json["args"])
        with open(args_values_file_path, "w", encoding="utf-8") as f:
            f.write(result_json["args_values"])
        print(f"    -> Successfully saved '{prog_name}'")
        return (prog_name, "success")
    else:
        print(f"    -> Failed to generate code for '{prog_name}'.")
        if result_json:
             print(f"       Missing keys in JSON: {list(result_json.keys())}")
        return (prog_name, "failed")

def main():
    if not API_KEY or API_KEY == '...':
        print("Error: API_KEY is not set in the script."); sys.exit(1)
    MAX_WORKERS = 10
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    tasks = []
    print("Collecting tasks to generate...")
    for theme, (theme_desc, programs) in all_benchmarks.items():
        theme_dir = os.path.join(OUTPUT_DIR, theme.replace(" ", "_").title())
        os.makedirs(theme_dir, exist_ok=True)
        for prog_name, params in programs.items():
            tasks.append((theme_desc, prog_name, params))
    print(f"Found {len(tasks)} total tasks to process.\n")
    futures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        print(f"Submitting tasks to a pool of {MAX_WORKERS} workers...")
        for theme_desc, prog_name, params in tasks:
            found_theme = next(t for t, (td, p) in all_benchmarks.items() if prog_name in p)
            future = executor.submit(process_benchmark, found_theme, theme_desc, prog_name, params)
            futures.append(future)
        print("\n--- Waiting for results ---")
        success_count, skipped_count, failed_count = 0, 0, 0
        for future in concurrent.futures.as_completed(futures):
            try:
                _, status = future.result()
                if status == "success": success_count += 1
                elif status == "skipped": skipped_count += 1
                else: failed_count += 1
            except Exception as e:
                print(f"A task generated an unexpected exception: {e}"); failed_count += 1
    print("\n--- Generation Complete ---")
    print(f"Summary: {success_count} succeeded, {failed_count} failed, {skipped_count} skipped.")

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")
