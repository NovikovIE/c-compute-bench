import os
import sys
import subprocess
import time
import argparse
import csv
import re
import signal
import platform

# --- КОНФИГУРАЦИЯ ---
BASE_DIR = os.getcwd()
BENCHMARKS_ROOT = os.path.join(BASE_DIR, "benchmarks")
BUILD_DIR = os.path.join(BASE_DIR, "build")
COMPILER = "gcc-13" 

def run_command(command, timeout_sec=None):
    # Эта функция без изменений
    try:
        process = subprocess.run(command, capture_output=True, text=True, timeout=timeout_sec)
        return process, None
    except subprocess.TimeoutExpired:
        return None, "TIMEOUT"
    except Exception as e:
        return None, str(e)

def sanitize_cflags_for_filename(cflags):
    # Эта функция без изменений
    sanitized = re.sub(r'^-', '', cflags)
    sanitized = re.sub(r'[\s/=-]+', '_', sanitized)
    return sanitized

def discover_benchmarks(root_dir, specific_program=None):
    """Автоматически находит все валидные бенчмарки."""
    tasks = []
    print(f"🔍 Discovering benchmarks in '{root_dir}'...")
    if not os.path.isdir(root_dir):
        print(f"Error: Benchmark root directory not found at '{root_dir}'"); return []
    for theme_dir in os.listdir(root_dir):
        theme_path = os.path.join(root_dir, theme_dir)
        if not os.path.isdir(theme_path): continue
        for prog_name in os.listdir(theme_path):
            if specific_program and prog_name != specific_program: continue
            prog_path = os.path.join(theme_path, prog_name)
            if not os.path.isdir(prog_path): continue
            c_file = os.path.join(prog_path, f"{prog_name}.c")
            args_values_file = os.path.join(prog_path, "args_values.txt")
            if os.path.exists(c_file) and os.path.exists(args_values_file):
                try:
                    with open(args_values_file, 'r') as f:
                        args_values = [arg.strip() for arg in f.read().strip().split(',')]
                    tasks.append({
                        "theme": theme_dir.replace("_", " ").title(), "prog_name": prog_name,
                        "c_file_path": c_file, "args_values": args_values
                    })
                except Exception as e:
                    print(f"  - Warning: Could not process {prog_name}. Reason: {e}")
    print(f"✅ Found {len(tasks)} valid benchmarks.")
    return tasks

# ---> ГЛАВНЫЕ ИЗМЕНЕНИЯ В ЭТОЙ ФУНКЦИИ <---
def run_single_test(task, cflags, force_compile, timeout, sanitizer=None):
    """
    Компилирует и выполняет один тест. Для санитайзеров на macOS
    автоматически переключается на нативный компилятор 'clang'.
    """
    
    # --- НОВАЯ ЛОГИКА ВЫБОРА КОМПИЛЯТОРА И ФЛАГОВ ---
    local_compiler = COMPILER
    final_cflags_list = cflags.split()
    
    if sanitizer:
        # Для санитайзеров всегда используем -O1 -g для лучшей отладки
        final_cflags_list = ["-O1", "-g", f"-fsanitize={sanitizer}"]

        # Если мы на macOS, переключаемся на clang, который нативно поддерживает санитайзеры
        if platform.system() == "Darwin":
            print(f"\n[Info] On macOS, switching to 'clang' for sanitizer run '{sanitizer}'.")
            local_compiler = "clang"
            # Для clang на macOS НЕ НУЖНЫ никакие дополнительные флаги типа -L или -static-lib...
            
    final_cflags_str = " ".join(final_cflags_list)

    prog_name, theme = task['prog_name'], task['theme']
    # Уменьшим количество выводимой информации, чтобы не загромождать лог
    # print(f"--- Testing: [{theme}] / [{prog_name}] with CFLAGS='{final_cflags_str}' ---")
    
    sanitized_flags = sanitize_cflags_for_filename(final_cflags_str)
    executable_file = os.path.join(BUILD_DIR, f"{prog_name}_{sanitized_flags}")

    # --- 1. Компиляция ---
    if not force_compile and os.path.exists(executable_file) and os.path.getmtime(executable_file) >= os.path.getmtime(task['c_file_path']):
        # print("  [1/2] Skipping compilation (up-to-date)")
        pass
    else:
        print(f"--- Compiling: [{prog_name}] with '{local_compiler}' and CFLAGS='{final_cflags_str}'")
        # Используем `local_compiler`, который может быть 'gcc-13' или 'clang'
        compile_cmd = [local_compiler] + final_cflags_list + [task['c_file_path'], "-o", executable_file, "-lm"]
        process, err_str = run_command(compile_cmd)
        if err_str or process.returncode != 0:
            print("  -> FAILED")
            reason = err_str or f"Compiler failed with code {process.returncode}:\n{process.stderr}"
            return {"status": "FAILED", "stage": "Compilation", "reason": reason}
        print("  -> OK")
    
    # --- 2. Выполнение и Анализ (остальная часть функции без изменений) ---
    print(f"--- Executing: [{prog_name}]")
    exec_cmd = [executable_file] + task['args_values']
    
    run_timeout = timeout * 10 if sanitizer else timeout
    process, err_str = run_command(exec_cmd, timeout_sec=run_timeout)

    if err_str == "TIMEOUT":
        print("  -> TIMED OUT")
        return {"status": "TIMED OUT", "stage": "Execution", "reason": f"Exceeded {run_timeout}s limit"}
    if not process:
        print("  -> FAILED (Command Error)"); return {"status": "FAILED", "stage": "Execution", "reason": f"Subprocess failed: {err_str}"}

    if sanitizer and process.returncode != 0 and process.stderr:
        print("  -> SANITIZER ERROR")
        return {"status": "FAILED (SANITIZER)", "stage": "Execution", "reason": process.stderr.strip(), "exit_code": process.returncode}
        
    if process.returncode != 0:
        if process.returncode < 0:
            print("  -> CRASHED")
            signal_name = "Unknown Signal"
            try: signal_name = signal.Signals(-process.returncode).name
            except ValueError: pass
            return {"status": "CRASHED", "stage": "Execution", "reason": f"Crashed with signal {-process.returncode} ({signal_name})", "exit_code": process.returncode}
        else:
            print("  -> FAILED")
            return {"status": "FAILED", "stage": "Execution", "reason": f"Exited with code {process.returncode}. Stderr: {process.stderr.strip()}", "exit_code": process.returncode}

    if sanitizer:
        print("  -> OK (No errors found)")
        return {"status": "SUCCESS", "result": "Sanitizer run completed without errors."}

    try:
        computation_time = float(process.stderr.strip())
        print(f"  -> OK (Time: {computation_time:.4f}s)")
        return {"status": "SUCCESS", "time_seconds": computation_time, "result": process.stdout.strip().replace('\n', ' ')}
    except (ValueError, TypeError):
        print("  -> FAILED (Bad Output)")
        return {"status": "FAILED", "stage": "Execution", "reason": f"Could not parse time from stderr: '{process.stderr.strip()}'"}


def main():
    parser = argparse.ArgumentParser(description="Compile and run C benchmarks.", formatter_class=argparse.RawTextHelpFormatter)
    
    # ---> НОВЫЙ АРГУМЕНТ <---
    parser.add_argument(
        "--sanitizer",
        choices=["address", "undefined", "thread", "leak"],
        help="Run tests with a specific sanitizer instead of performance measurement."
    )
    # Остальные аргументы без изменений
    parser.add_argument("--cflags", nargs='+', default=["-O0", "-O2", "-O3"], help="Compiler flags sets to use for performance runs.")
    parser.add_argument("-p", "--program", help="Run a single benchmark.")
    parser.add_argument("-f", "--force-compile", action="store_true", help="Force recompilation.")
    parser.add_argument("-o", "--output", default="benchmark_results.csv", help="Output CSV file.")
    parser.add_argument("--timeout", type=float, default=20.0, help="Execution timeout for performance runs.")
    args = parser.parse_args()

    # ---> НОВАЯ ЛОГИКА ВЫБОРА РЕЖИМА <---
    if args.sanitizer:
        print(f"Sanitizer mode enabled: address={args.sanitizer}")
        print("Performance measurement will be disabled. CFLAGS will be overridden.")
        # Для санитайзера достаточно одного "запуска"
        cflags_to_run = [f"-fsanitize={args.sanitizer}"]
    else:
        cflags_to_run = args.cflags

    tasks = discover_benchmarks(BENCHMARKS_ROOT, args.program)
    if not tasks: sys.exit(1)

    os.makedirs(BUILD_DIR, exist_ok=True)
    
    try:
        csv_file = open(args.output, 'w', newline='', encoding='utf-8')
        header = ["theme", "program", "cflags", "status", "time_seconds", "result", "stage", "exit_code", "reason"]
        csv_writer = csv.DictWriter(csv_file, fieldnames=header)
        csv_writer.writeheader()
        print(f"\nLogging results to '{args.output}'")
    except IOError as e:
        print(f"FATAL: Could not open output file '{args.output}'. Reason: {e}"); sys.exit(1)

    all_results = []
    for task in tasks:
        # Теперь цикл идет по `cflags_to_run`, который зависит от режима
        for cflag_set in cflags_to_run:
            result = run_single_test(task, cflag_set, args.force_compile, args.timeout, args.sanitizer)
            full_result = {"theme": task['theme'], "program": task['prog_name'], "cflags": cflag_set, **result}
            for key in header:
                if key not in full_result: full_result[key] = ''
            all_results.append(full_result)
            csv_writer.writerow(full_result)
    csv_file.close()

    # Итоговая таблица без изменений
    print("\n" + "="*100)
    print(" " * 40 + "BENCHMARK SUMMARY")
    print("="*100)
    print(f"{'STATUS':<20} {'THEME':<20} {'PROGRAM':<25} {'CFLAGS/SANITIZER':<20} {'DETAILS'}")
    print("-" * 100)
    for res in sorted(all_results, key=lambda x: (x['theme'], x['program'], x['cflags'])):
        details = ""
        if res['status'] == "SUCCESS":
            details = f"Time: {res.get('time_seconds', 'N/A')}s | Result: {str(res.get('result', ''))[:30]}"
        else:
            reason_short = str(res.get('reason', 'Unknown')).splitlines()[0]
            details = f"[{res.get('stage', 'N/A')}] {reason_short[:50]}"
        print(f"[{res['status']:<18}] {res['theme']:<20} {res['program']:<25} {res['cflags']:<20} {details}")
    print("="*100); print(f"📈 Detailed report saved to '{args.output}'")


if __name__ == "__main__":
    main()
