import time
import collections
import dataclasses
import struct
import numpy as np
import random
import sys

# =============================================================================
# --- Benchmark Configuration ---
# =============================================================================
# Data size for collections used in tests
DATA_SIZE = 1_000

# Number of repetitions for each operation.
# Adjusted per category because some operations are much slower than others.
REPEATS_FAST_OPS = 100_000   # For very quick, in-place operations
REPEATS_MEDIUM_OPS = 10_000  # For operations creating lists/objects
REPEATS_SLOW_OPS = 100       # For very slow operations like a bad sort

print(f"--- CPython Performance Benchmark (Direct Timing Method) ---")
print(f"Python Version: {sys.version}")
print(f"Numpy Version: {np.__version__}")
print(f"Configuration: Data Size = {DATA_SIZE}\n")

# =============================================================================
# --- Helper Function to Run and Print Benchmarks ---
# =============================================================================
def run_and_print_results(title, benchmarks, *args):
    """
    A helper to run benchmark functions, time them, and print formatted results.
    'benchmarks' is a list of tuples: (name, function_to_run).
    """
    print(f"--- {title} ---")
    results = {}
    for name, func in benchmarks:
        # Simple warm-up run to account for any initial overhead
        try:
            func(*args)
            
            # Actual timed run
            start_time = time.perf_counter()
            func(*args)
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            results[name] = elapsed_time
        except Exception as e:
            results[name] = -1 # Indicate failure
            print(f"  Warning: Benchmark '{name}' failed with error: {e}")


    # Find the fastest result to use as a baseline for comparison
    valid_results = {k: v for k, v in results.items() if v >= 0}
    if not valid_results:
        print("  All benchmarks in this category failed.")
        print("\n")
        return

    baseline_name = min(valid_results, key=valid_results.get)
    baseline_time = valid_results[baseline_name]

    # Print sorted results
    for name, timing in sorted(results.items(), key=lambda item: item[1]):
        if timing < 0:
            print(f"  {name:<40}: FAILED")
            continue
        relative_speed = f"({timing / baseline_time:.2f}x slower)" if name != baseline_name else "(baseline)"
        print(f"  {name:<40}: {timing:.6f} seconds {relative_speed}")
    print("\n")


# =============================================================================
# 1. Loops vs. Comprehensions
# =============================================================================

# 1a. Single list creation
def bench_for_loop_append(data, n_repeats):
    for _ in range(n_repeats):
        new_list = []
        for i in data:
            new_list.append(i * 2)

def bench_list_comprehension(data, n_repeats):
    for _ in range(n_repeats):
        new_list = [i * 2 for i in data]

def bench_set_comprehension(data, n_repeats):
    for _ in range(n_repeats):
        new_set = {i * 2 for i in data}

# 1b. Multiple list creation
def bench_single_loop_multiple_appends(data, n_repeats):
    for _ in range(n_repeats):
        list_a, list_b, list_c = [], [], []
        for i in data:
            list_a.append(i * 2)
            list_b.append(i * 3)
            list_c.append(i * 4)

def bench_multiple_comprehensions(data, n_repeats):
    for _ in range(n_repeats):
        list_a = [i * 2 for i in data]
        list_b = [i * 3 for i in data]
        list_c = [i * 4 for i in data]

# =============================================================================
# 2. Heterogeneous Data Structures (Struct-like objects)
# =============================================================================

# --- Define the structures first ---
class MyClass:
    def __init__(self, id, name, value):
        self.id, self.name, self.value = id, name, value

@dataclasses.dataclass
class MyDataClass:
    id: int
    name: str
    value: float

NamedTuple = collections.namedtuple('NamedTuple', ['id', 'name', 'value'])
numpy_dtype = np.dtype([('id', 'i4'), ('name', 'U10'), ('value', 'f8')])
struct_format = 'i10sd'
packer = struct.Struct(struct_format)

# --- Creation Benchmarks ---
def bench_create_dict(n_repeats):
    for _ in range(n_repeats):
        structs = [{'id': i, 'name': 'test', 'value': i * 3.14} for i in range(DATA_SIZE)]

def bench_create_tuple(n_repeats):
    for _ in range(n_repeats):
        structs = [(i, 'test', i * 3.14) for i in range(DATA_SIZE)]

def bench_create_namedtuple(n_repeats):
    for _ in range(n_repeats):
        structs = [NamedTuple(i, 'test', i * 3.14) for i in range(DATA_SIZE)]

def bench_create_class(n_repeats):
    for _ in range(n_repeats):
        structs = [MyClass(i, 'test', i * 3.14) for i in range(DATA_SIZE)]

def bench_create_dataclass(n_repeats):
    for _ in range(n_repeats):
        structs = [MyDataClass(i, 'test', i * 3.14) for i in range(DATA_SIZE)]

# --- Access Benchmarks ---
def bench_access_dict(data, n_repeats):
    for _ in range(n_repeats):
        total = 0
        for s in data: total += s['value']

def bench_access_tuple(data, n_repeats):
    for _ in range(n_repeats):
        total = 0
        for s in data: total += s[2]

def bench_access_namedtuple(data, n_repeats):
    for _ in range(n_repeats):
        total = 0
        for s in data: total += s.value

def bench_access_class(data, n_repeats):
    for _ in range(n_repeats):
        total = 0
        for s in data: total += s.value

def bench_access_dataclass(data, n_repeats):
    for _ in range(n_repeats):
        total = 0
        for s in data: total += s.value

def bench_access_numpy(data, n_repeats):
    for _ in range(n_repeats):
        total = np.sum(data['value'])

# =============================================================================
# 3. Piecewise Function Evaluation
# =============================================================================
def piecewise_py(x):
    if x < -1: return -1.0
    elif -1 <= x <= 1: return x * x
    else: return 1.0

def bench_piecewise_python_loop(data, n_repeats):
    for _ in range(n_repeats):
        results = [piecewise_py(x) for x in data]

def bench_piecewise_numpy_masking(data, n_repeats):
    for _ in range(n_repeats):
        conds = [data < -1, (data >= -1) & (data <= 1), data > 1]
        funcs = [-1.0, lambda x: x*x, 1.0]
        results = np.piecewise(data, conds, funcs)

# =============================================================================
# 4. Built-in Functions vs. Manual Implementation
# =============================================================================
def is_even(n): return n % 2 == 0
def double(n): return n * 2

def bench_map_function(data, n_repeats):
    for _ in range(n_repeats):
        list(map(double, data))

def bench_map_comprehension(data, n_repeats):
    for _ in range(n_repeats):
        [double(x) for x in data]

def bench_filter_function(data, n_repeats):
    for _ in range(n_repeats):
        list(filter(is_even, data))

def bench_filter_comprehension(data, n_repeats):
    for _ in range(n_repeats):
        [x for x in data if is_even(x)]

def bench_builtin_sum(data, n_repeats):
    for _ in range(n_repeats):
        sum(data)

def bench_manual_sum(data, n_repeats):
    for _ in range(n_repeats):
        total = 0
        for x in data: total += x

def bench_builtin_sorted(data, n_repeats):
    for _ in range(n_repeats):
        sorted(data) # In-place sort would modify, so use sorted()

# =============================================================================
# 5. Function Call Overhead
# =============================================================================
def operation_in_func(x, y): return x * y + 5

def bench_inline_operation(data, n_repeats):
    for _ in range(n_repeats):
        result = 0
        for i in data:
            result += i * 2 + 5

def bench_function_call_operation(data, n_repeats):
    for _ in range(n_repeats):
        result = 0
        for i in data:
            result += operation_in_func(i, 2)


# =============================================================================
# --- MAIN EXECUTION BLOCK ---
# =============================================================================
if __name__ == "__main__":
    # --- Setup common data ---
    source_data = list(range(DATA_SIZE))
    numpy_data = np.linspace(-3, 3, DATA_SIZE)

    # --- Run Suite 1: Loops vs Comprehensions ---
    run_and_print_results(
        "1a. Loop vs. Comprehension (Single Collection)",
        [
            ("For Loop with .append()", bench_for_loop_append),
            ("List Comprehension", bench_list_comprehension),
            ("Set Comprehension", bench_set_comprehension),
        ],
        source_data, REPEATS_MEDIUM_OPS
    )

    run_and_print_results(
        "1b. Single Loop vs. Multiple Comprehensions",
        [
            ("Single Loop, Multiple Appends", bench_single_loop_multiple_appends),
            ("Multiple Comprehensions", bench_multiple_comprehensions),
        ],
        source_data, REPEATS_MEDIUM_OPS // 2 # This is slower, so fewer repeats
    )

    # --- Run Suite 2: Data Structures ---
    run_and_print_results(
        "2a. Data Structure Creation",
        [
            ("Dictionary", bench_create_dict),
            ("Tuple", bench_create_tuple),
            ("NamedTuple", bench_create_namedtuple),
            ("Standard Class", bench_create_class),
            ("DataClass", bench_create_dataclass),
        ],
        REPEATS_MEDIUM_OPS
    )
    
    # Pre-create data for access tests to be fair
    dict_data = [{'id': i, 'name': 'test', 'value': i * 3.14} for i in range(DATA_SIZE)]
    tuple_data = [(i, 'test', i * 3.14) for i in range(DATA_SIZE)]
    namedtuple_data = [NamedTuple(i, 'test', i * 3.14) for i in range(DATA_SIZE)]
    class_data = [MyClass(i, 'test', i * 3.14) for i in range(DATA_SIZE)]
    dataclass_data = [MyDataClass(i, 'test', i * 3.14) for i in range(DATA_SIZE)]
    np_struct_data = np.array(tuple_data, dtype=numpy_dtype)

    run_and_print_results(
        "2b. Data Structure Field Access",
        [
            ("Dictionary Access", lambda: bench_access_dict(dict_data, REPEATS_MEDIUM_OPS)),
            ("Tuple Access", lambda: bench_access_tuple(tuple_data, REPEATS_MEDIUM_OPS)),
            ("NamedTuple Access", lambda: bench_access_namedtuple(namedtuple_data, REPEATS_MEDIUM_OPS)),
            ("Standard Class Access", lambda: bench_access_class(class_data, REPEATS_MEDIUM_OPS)),
            ("DataClass Access", lambda: bench_access_dataclass(dataclass_data, REPEATS_MEDIUM_OPS)),
            ("Numpy Structured Array Access", lambda: bench_access_numpy(np_struct_data, REPEATS_MEDIUM_OPS)),
        ]
    )


    # --- Run Suite 3: Piecewise Function ---
    run_and_print_results(
        "3. Piecewise Function Evaluation",
        [
            ("Python Loop with if/else", bench_piecewise_python_loop),
            ("Numpy Boolean Masking", bench_piecewise_numpy_masking),
        ],
        numpy_data, REPEATS_SLOW_OPS
    )

    # --- Run Suite 4: Built-ins ---
    run_and_print_results(
        "4. Built-in Functions vs. Manual",
        [
            ("map() function", bench_map_function),
            ("map() via list comprehension", bench_map_comprehension),
            ("filter() function", bench_filter_function),
            ("filter() via list comprehension", bench_filter_comprehension),
            ("Built-in sum()", bench_builtin_sum),
            ("Manual sum() with for loop", bench_manual_sum),
        ],
        source_data, REPEATS_MEDIUM_OPS
    )
    # Sorting is very slow, so it gets its own run with fewer repeats
    data_for_sorting = random.sample(range(DATA_SIZE * 10), DATA_SIZE)
    run_and_print_results(
        "4b. Sorting",
        [
            ("Built-in sorted()", bench_builtin_sorted),
        ],
        data_for_sorting, REPEATS_SLOW_OPS
    )


    # --- Run Suite 5: Function Call Overhead ---
    run_and_print_results(
        "5. Function Call Overhead",
        [
            ("Operation in-line in loop", bench_inline_operation),
            ("Operation inside a called function", bench_function_call_operation),
        ],
        source_data, REPEATS_MEDIUM_OPS
    )