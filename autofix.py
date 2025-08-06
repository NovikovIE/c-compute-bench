import os
import sys
import csv
import json
import requests
import shutil
from datetime import datetime
import textwrap
import concurrent.futures

# --- КОНФИГУРАЦИЯ ---
API_KEY = os.getenv("SOY_TOKEN")
API_URL = "https://api.eliza.yandex.net/openrouter/v1/chat/completions"
MODEL_NAME = "google/gemini-2.5-pro-preview-03-25"
MAX_WORKERS = 8
HEADERS = {
    "Authorization": f"OAuth {API_KEY}",
    "Content-Type": "application/json"
}
try:
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
except ImportError:
    pass

BENCHMARKS_ROOT = "benchmarks"
RESULTS_CSV_PATH = "benchmark_results.csv"

def create_fix_prompt(theme, prog_name, c_code, args_str, args_values_str, error_message):
    """
    Создает специализированный промпт для исправления, понимающий НОВУЮ структуру бенчмарков.
    """
    # Эти требования - точная копия из нового скрипта генерации.
    original_requirements = textwrap.dedent("""\
        1.  **Structure:** The C code MUST be organized into `setup_benchmark()`, `run_computation()`, and `cleanup()` functions.
        2.  **Timing:** The `main` function must time ONLY the `run_computation()` call using `clock_gettime(CLOCK_MONOTONIC, ...)`.
        3.  **Output:** The program must print the final computational result to `stdout` and the measured time of `run_computation()` to `stderr`.
        4.  **Data:** All data must be generated internally using the provided Mersenne Twister generator, seeded from an `argv` parameter. Large data must be allocated on the heap.
        5.  **Goal:** The program must compile with `gcc` and finish with a `0` exit code. Crashing with signals like -11 (Segfault) is a critical failure.
    """)

    prompt_template = textwrap.dedent("""\
        You are an expert C language debugger and a JSON creating machine. Your task is to analyze a buggy C benchmark program that failed during testing, fix it, and return ONLY the corrected C code in a specific JSON format.

        --- PROBLEM TO SOLVE ---
        - **Program Name:** `{prog_name}`
        - **Theme:** {theme}
        - **Failure Reason:**
        ```
        {error_message}
        ```
        - **Context - Arguments used during crash:**
          - Argument Names (`args.txt`): `{args_str}`
          - Argument Values (`args_values.txt`): `{args_values_str}`

        - **ORIGINAL (BUGGY) C CODE (`{prog_name}.c`):**
        ```c
        {c_code}
        ```
        
        - **ORIGINAL BENCHMARK REQUIREMENTS (MUST BE PRESERVED):**
        {original_requirements}

        **YOUR TASK:**
        1.  Analyze the bug in the C code, considering the error message and the arguments. Common bugs are null pointer dereferences, out-of-bounds array access, incorrect `malloc` sizes, or logical errors.
        2.  Fix the C code. The fix must be robust and strictly adhere to all the original benchmark requirements listed above. Do not alter the program's structure (`setup`/`compute`/`cleanup`) or its timing/output mechanism.
        3.  Return your response as a single, raw JSON object containing ONLY the full, corrected C code. DO NOT ADD ANY ANALYZE OR DESCRIPTION BEFORE JSON, ADD ONLY IT

        **EXAMPLE OF THE REQUIRED OUTPUT FORMAT:**
        ```json
        {{
          "fixed_c_code": "#include <stdio.h>\\n#include <stdlib.h>\\n// ... full, corrected C code as a single escaped string ...\\nint main(int argc, char *argv[]) {{\\n    // ... fixed logic ...\\n    return 0;\\n}}"
        }}
        ```
    """)
    
    return prompt_template.format(
        theme=theme, prog_name=prog_name, error_message=error_message,
        args_str=args_str, args_values_str=args_values_str, c_code=c_code,
        original_requirements=original_requirements
    )

def make_api_request(prompt_content, prog_name):
    # Эта функция без изменений
    payload = {"model": MODEL_NAME, "messages": [{"role": "user", "content": prompt_content}], "temperature": 0.1}
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
        print(f"\n  [{prog_name}] [Error] Network/API error: {e}")
    except json.JSONDecodeError as e:
        print(f"\n  [{prog_name}] [Error] Failed to parse JSON: {e}\n[Raw Content]: {generated_content_str[:500]}...")
    except (KeyError, IndexError) as e:
        print(f"\n  [{prog_name}] [Error] Unexpected API response: {e}\n[Raw Response]: {response.text}")
    return None

def read_failed_benchmarks(csv_path):
    # Эта функция без изменений, так как статусы ошибок те же
    failed_tasks = []
    try:
        with open(csv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['status'] in ["CRASHED", "FAILED", 'FAILED (SANITIZER)']:
                    failed_tasks.append({
                        "theme": row['theme'], "program": row['program'],
                        "cflags": row['cflags'], "reason": row['reason']
                    })
    except FileNotFoundError:
        print(f"FATAL: Results file not found at '{csv_path}'."); sys.exit(1)
    return failed_tasks

def fix_benchmark_task(task):
    # Эта функция без изменений, она универсальна
    theme, prog_name = task['theme'], task['program']
    print(f"  -> Starting fix for [{prog_name}] (failed with '{task['cflags']}')")
    
    prog_dir = os.path.join(BENCHMARKS_ROOT, theme.replace(" ", "_").title(), prog_name)
    c_file_path = os.path.join(prog_dir, f"{prog_name}.c")
    args_file_path = os.path.join(prog_dir, "args.txt")
    args_values_file_path = os.path.join(prog_dir, "args_values.txt")
    
    try:
        with open(c_file_path, 'r', encoding='utf-8') as f: c_code = f.read()
        with open(args_file_path, 'r', encoding='utf-8') as f: args_str = f.read().strip()
        with open(args_values_file_path, 'r', encoding='utf-8') as f: args_values_str = f.read().strip()
    except FileNotFoundError as e:
        print(f"  [{prog_name}] [Error] Could not find source file: {e}."); return False
        
    error_message = task['reason'][:2000]
    prompt = create_fix_prompt(theme, prog_name, c_code, args_str, args_values_str, error_message)
    api_result = make_api_request(prompt, prog_name)

    if api_result and "fixed_c_code" in api_result and api_result["fixed_c_code"]:
        backup_path = f"{c_file_path}.{datetime.now().strftime('%Y%m%d_%H%M%S')}.bak"
        shutil.copy(c_file_path, backup_path)
        with open(c_file_path, 'w', encoding='utf-8') as f:
            f.write(api_result["fixed_c_code"])
        print(f"  <- SUCCESS: Applied fix for '{prog_name}'. Backup saved.")
        return True
    else:
        print(f"  <- FAILED: No valid fix for '{prog_name}'. File not changed.")
        return False

def main():
    # Эта функция без изменений, она универсальна
    if not API_KEY or API_KEY == '...':
        print("FATAL: API_KEY is not set."); sys.exit(1)
    failed_tasks_raw = read_failed_benchmarks(RESULTS_CSV_PATH)
    if not failed_tasks_raw:
        print("✅ Congratulations! No 'CRASHED' or 'FAILED' benchmarks found."); return
    tasks_to_fix = list({task['program']: task for task in failed_tasks_raw}.values())
    print(f"Found {len(tasks_to_fix)} unique programs that need fixing:")
    for i, task in enumerate(tasks_to_fix):
        print(f"  {i+1}. {task['program']} (example failure with cflags '{task['cflags']}')")
    try:
        choice = input(f"\nAttempt to fix these {len(tasks_to_fix)} programs in parallel? (yes/no): ").lower()
        if choice != 'yes': print("Aborting."); return
    except KeyboardInterrupt: print("\nAborting."); return

    success_fixes, failed_fixes = 0, 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_task = {executor.submit(fix_benchmark_task, task): task for task in tasks_to_fix}
        print(f"\n--- Submitting {len(tasks_to_fix)} jobs to the pool. Waiting for results... ---")
        for future in concurrent.futures.as_completed(future_to_task):
            task_name = future_to_task[future]['program']
            try:
                if future.result(): success_fixes += 1
                else: failed_fixes += 1
            except Exception as exc:
                print(f"  [{task_name}] Generated an unexpected exception: {exc}"); failed_fixes += 1
    print("\n--- Auto-Fixing Complete ---")
    print(f"Summary: {success_fixes} fixes applied, {failed_fixes} fixes failed.")
    print("\nRECOMMENDATION: Re-run './run_all.py' to verify the fixes.")

if __name__ == "__main__":
    main()

