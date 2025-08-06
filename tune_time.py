import os
import sys
import subprocess
import time
import argparse
import csv
import re
import signal
import json
import shutil
from datetime import datetime
import textwrap
import concurrent.futures
import requests
import re 
# --- ОБЩАЯ КОНФИГУРАЦИЯ ---
API_KEY = os.getenv("SOY_TOKEN")
API_URL = "https://api.eliza.yandex.net/openrouter/v1/chat/completions"
MODEL_NAME = "google/gemini-2.5-flash"
BASE_DIR = os.getcwd()
BENCHMARKS_ROOT = os.path.join(BASE_DIR, "benchmarks")
BUILD_DIR = os.path.join(BASE_DIR, "build")
COMPILER = "gcc-13"

# --- КОНФИГУРАЦИЯ ИТЕРАТИВНОГО ПРОЦЕССА ---
MAX_REFINEMENT_ATTEMPTS = 10
MAX_API_RETRIES = 5
MIN_TIME, MAX_TIME = 0.5, 1.5
MAX_WORKERS = 8 # Параллельное уточнение для N бенчмарков
CFLAGS_FOR_REFINEMENT = "-O3" # Используем -O3 как эталон для замера времени

HEADERS = {"Authorization": f"OAuth {API_KEY}", "Content-Type": "application/json"}
try:
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
except ImportError:
    pass

# --- ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ (из предыдущих скриптов) ---

def run_command(command, timeout_sec=None):
    try:
        return subprocess.run(command, capture_output=True, text=True, timeout=timeout_sec), None
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT"
    except Exception as e:
        return None, str(e)

def sanitize_cflags_for_filename(cflags):
    sanitized = re.sub(r'^-', '', cflags)
    return re.sub(r'[\s/=-]+', '_', sanitized)

def discover_benchmarks(root_dir, specific_program=None):
    # Эта функция без изменений, но нужна здесь
    tasks = []
    if not os.path.isdir(root_dir): return []
    for theme_dir in os.listdir(root_dir):
        theme_path = os.path.join(root_dir, theme_dir)
        if not os.path.isdir(theme_path): continue
        for prog_name in os.listdir(theme_path):
            if specific_program and prog_name != specific_program: continue
            prog_path = os.path.join(theme_path, prog_name)
            c_file = os.path.join(prog_path, f"{prog_name}.c")
            if os.path.exists(c_file):
                tasks.append({"theme": theme_dir.replace("_", " ").title(), "prog_name": prog_name})
    return tasks

def run_single_test(task, cflags, timeout):
    # Слегка упрощенная версия из run_all.py
    prog_name = task['prog_name']
    prog_dir = os.path.join(BENCHMARKS_ROOT, task['theme'].replace(" ", "_").title(), prog_name)
    c_file_path = os.path.join(prog_dir, f"{prog_name}.c")
    args_values_file = os.path.join(prog_dir, "args_values.txt")

    try:
        with open(args_values_file, 'r') as f:
            args_values = [arg.strip() for arg in f.read().strip().split(',')]
    except FileNotFoundError:
        return {"status": "FAILED", "stage": "Setup", "reason": f"{args_values_file} not found."}

    executable_file = os.path.join(BUILD_DIR, f"{prog_name}_{sanitize_cflags_for_filename(cflags)}")
    
    # Компиляция
    compile_cmd = [COMPILER] + cflags.split() + [c_file_path, "-o", executable_file, "-lm"]
    process, err_str = run_command(compile_cmd)
    if err_str or process.returncode != 0:
        return {"status": "FAILED", "stage": "Compilation", "reason": err_str or f"GCC Error:\n{process.stderr}"}
    
    # Выполнение
    exec_cmd = [executable_file] + args_values
    process, err_str = run_command(exec_cmd, timeout_sec=timeout)

    # Анализ результата
    if err_str == "TIMEOUT": return {"status": "TIMED OUT"}
    if not process: return {"status": "FAILED", "stage": "Execution", "reason": f"Subprocess failed: {err_str}"}
    if process.returncode != 0:
        if process.returncode < 0:
            try: signal_name = signal.Signals(-process.returncode).name
            except ValueError: signal_name = "Unknown Signal"
            return {"status": "CRASHED", "reason": f"Signal {-process.returncode} ({signal_name})"}
        else:
            return {"status": "FAILED", "reason": f"Exit code {process.returncode}. Stderr: {process.stderr.strip()}"}
    
    try:
        computation_time = float(process.stderr.strip())
        return {"status": "SUCCESS", "time_seconds": computation_time}
    except (ValueError, TypeError):
        return {"status": "FAILED", "reason": f"Could not parse time from stderr: '{process.stderr.strip()}'"}

# --- НОВЫЕ ФУНКЦИИ ДЛЯ ИТЕРАТИВНОЙ ПОДСТРОЙКИ ---

def create_refinement_prompt(prog_name, c_code, args_str, run_history):
    """Создает промпт, который передает модели историю попыток."""
    history_str = ""
    for i, run in enumerate(run_history):
        result_desc = ""
        if run['status'] == 'SUCCESS':
            time = run['time_seconds']
            if time < MIN_TIME: result_desc = f"{time:.4f}s (Too Fast)"
            elif time > MAX_TIME: result_desc = f"{time:.4f}s (Too Slow)"
            else: result_desc = f"{time:.4f}s (SUCCESS)" # Should not happen in loop
        else:
            result_desc = f"{run['status']}: {run.get('reason', 'N/A')}"
        history_str += f"  - Attempt {i+1}: Args = {run['args_values']}, Result: {result_desc}\n"
    
    prompt = textwrap.dedent(f"""
        You are an intelligent performance tuning engineer. Your task is to iteratively find the right arguments for a C benchmark to make its `run_computation` phase take about 1.0 second.

        **Target:**
        - **Goal:** Execution time between {MIN_TIME:.2f}s and {MAX_TIME:.2f}s.
        - **Program:** `{prog_name}`
        - **Argument Names:** `{args_str}`

        **History of Previous Attempts:**
        {history_str}

        **C Source Code (for complexity analysis):**
        ```c
        {c_code}
        ```

        **YOUR TASK:**
        Based on the history, analyze the code's complexity and propose the **next** set of argument values to try.
        - If the last attempt was too fast or slow, extrapolate/interpolate a better value.
        - If the last attempt CRASHED or FAILED, the arguments were likely invalid (e.g., too large, causing memory errors, or an edge case). Propose a more conservative/safer set of values.
        - Do NOT change the `seed` value.
        - Respond ONLY with a JSON object containing the new comma-separated values.

        **Required JSON Output Format:**
        ```json
        {{
          "new_args_values": "..."
        }}
        ```
    """)
    return prompt

def make_api_request(prompt_content, prog_name):
    """
    Отправляет запрос к API с экспоненциальной задержкой и УМНЫМ парсингом JSON.
    Находит JSON блок ```json ... ``` в любом месте ответа.
    """
    payload = {"model": MODEL_NAME, "messages": [{"role": "user", "content": prompt_content}], "temperature": 0.4}
    raw_content_str = ""

    for attempt in range(MAX_API_RETRIES):
        try:
            response = requests.post(API_URL, headers=HEADERS, json=payload, timeout=60, verify=False)
            response.raise_for_status()

            response_json = response.json()
            raw_content_str = response_json.get("response", {}).get("choices", [{}])[0].get("message", {}).get("content", "")
            
            if not raw_content_str:
                print(f"\n  [{prog_name}] API returned empty content on attempt {attempt + 1}. Retrying...")
                continue
                
            # ---> НОВАЯ, БОЛЕЕ НАДЕЖНАЯ ЛОГИКА ПАРСИНГА <---
            # Ищем блок, заключенный в ```json ... ```
            match = re.search(r"```json\s*(\{.*?\})\s*```", raw_content_str, re.DOTALL)
            
            json_text = ""
            if match:
                # Если нашли блок, берем его содержимое
                json_text = match.group(1)
            else:
                # Если блока нет, может быть, модель вернула чистый JSON. Пробуем парсить весь ответ.
                json_text = raw_content_str

            return json.loads(json_text)

        except requests.exceptions.RequestException as e:
            print(f"\n  [{prog_name}] API request failed: {e}. ({attempt + 1}/{MAX_API_RETRIES})")
        
        except json.JSONDecodeError as e:
            print(f"\n  [{prog_name}] [FATAL API ERROR] Failed to parse JSON. Error: {e}")
            print(f"  RAW MODEL OUTPUT FOR {prog_name}:\n---\n{raw_content_str}\n---")
            return None
        
        except (KeyError, IndexError) as e:
            print(f"\n  [{prog_name}] [FATAL API ERROR] Unexpected response structure: {e}")
            print(f"  RAW API RESPONSE FOR {prog_name}:\n---\n{response.text}\n---")
            return None

    print(f"\n  [{prog_name}] API request failed after {MAX_API_RETRIES} retries.")
    return None

def refine_benchmark(task, timeout):
    """
    Выполняет итеративный процесс подстройки для одного бенчмарка.
    """
    prog_name = task['prog_name']
    prog_dir = os.path.join(BENCHMARKS_ROOT, task['theme'].replace(" ", "_").title(), prog_name)
    print(f"--- Refining [{prog_name}] ---")
    
    # Загружаем "статичные" данные один раз
    try:
        with open(os.path.join(prog_dir, f"{prog_name}.c"), 'r') as f: c_code = f.read()
        with open(os.path.join(prog_dir, "args.txt"), 'r') as f: args_names = f.read().strip()
    except FileNotFoundError as e:
        print(f"  [ERROR] Missing source file for {prog_name}: {e}"); return "SETUP_FAILED"

    run_history = []
    for attempt in range(1, MAX_REFINEMENT_ATTEMPTS + 1):
        print(f"  [Attempt {attempt}/{MAX_REFINEMENT_ATTEMPTS}] ", end='')
        
        # Запускаем тест с текущими параметрами
        result = run_single_test(task, CFLAGS_FOR_REFINEMENT, timeout)
        
        # Сохраняем результат в историю
        with open(os.path.join(prog_dir, "args_values.txt"), 'r') as f:
            current_args = f.read().strip()
        run_history.append({"args_values": current_args, **result})

        if result['status'] == 'SUCCESS' and MIN_TIME <= result['time_seconds'] <= MAX_TIME:
            print(f"Result: SUCCESS, Time: {result['time_seconds']:.4f}s (IN-RANGE)")
            print(f"✅ Refined [{prog_name}] successfully in {attempt} attempts.")
            return "SUCCESS"
        else:
            status_desc = f"{result['status']}, Time: {result.get('time_seconds', 'N/A'):.4f}s" if result.get('time_seconds') else result['status']
            print(f"Result: {status_desc} (OUT OF RANGE OR FAILED)")
        
        if attempt == MAX_REFINEMENT_ATTEMPTS: break # Последняя попытка, выходим

        # Готовим и отправляем запрос в API
        print("    -> Adjusting parameters via LLM...")
        prompt = create_refinement_prompt(prog_name, c_code, args_names, run_history)
        api_result = make_api_request(prompt, prog_name)

        if not api_result or "new_args_values" not in api_result:
            print(f"  [ERROR] Failed to get valid adjustment from LLM for {prog_name}")
            return "API_FAILED"
        
        new_values = api_result['new_args_values']
        print(f"    -> LLM proposed new args: {new_values}")
        
        # Обновляем файл с аргументами для следующей итерации
        args_values_file = os.path.join(prog_dir, "args_values.txt")
        shutil.copy(args_values_file, f"{args_values_file}.{attempt}.bak") # Делаем бэкап
        with open(args_values_file, 'w') as f:
            f.write(new_values)
            
    print(f"❌ Failed to refine [{prog_name}] within {MAX_REFINEMENT_ATTEMPTS} attempts.")
    return "FAILED_TO_REFINE"

def main():
    parser = argparse.ArgumentParser(description="Iteratively refine benchmark execution times.")
    parser.add_argument("-p", "--program", help="Refine only a single benchmark.")
    parser.add_argument("--timeout", type=float, default=30.0, help="Timeout per run.")
    args = parser.parse_args()

    if not API_KEY or API_KEY == '...':
        print("FATAL: API_KEY is not set."); sys.exit(1)
        
    tasks = discover_benchmarks(BENCHMARKS_ROOT, args.program)
    if not tasks: print("No benchmarks found to refine."); return
    
    os.makedirs(BUILD_DIR, exist_ok=True)
    print(f"Found {len(tasks)} benchmarks. Starting refinement process with CFLAGS='{CFLAGS_FOR_REFINEMENT}'...")
    
    final_statuses = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_task = {executor.submit(refine_benchmark, task, args.timeout): task for task in tasks}

        for future in concurrent.futures.as_completed(future_to_task):
            task_name = future_to_task[future]['prog_name']
            try:
                status = future.result()
                final_statuses[task_name] = status
            except Exception as e:
                print(f"An error occurred while refining {task_name}: {e}")
                final_statuses[task_name] = "EXCEPTION"

    print("\n" + "="*50)
    print(" " * 15 + "REFINEMENT SUMMARY")
    print("="*50)
    # Считаем итоги
    summary = {}
    for status in final_statuses.values():
        summary[status] = summary.get(status, 0) + 1
    
    for status, count in summary.items():
        print(f"- {status}: {count}")

    print("="*50)
    print("\nRECOMMENDATION: Now run 'run_all.py' to get the final performance results across all compiler flags.")

if __name__ == "__main__":
    main()
