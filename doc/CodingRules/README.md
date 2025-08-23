# Coding Rules and Workflows

This directory contains the project-wide coding rules and a set of lightweight, task-focused workflows. Start with `global_rules.md` and then pick the workflow matching your task.

- Global rules: see [`global_rules.md`](./global_rules.md)
- Workflow template: see [`workflow_template.md`](./workflow_template.md)

## Workflow index

- Python
  - [`workflows/python/python_demo_workflow.md`](./workflows/python/python_demo_workflow.md)
  - [`workflows/python/python_ctypes_workflow.md`](./workflows/python/python_ctypes_workflow.md)
- C/C++
  - [`workflows/cpp/cpp_debugging_workflow.md`](./workflows/cpp/cpp_debugging_workflow.md)
  - [`workflows/cpp/cpp_prototype_workflow.md`](./workflows/cpp/cpp_prototype_workflow.md)
- Documentation
  - [`workflows/documentation/extract_documentation_workflow.md`](./workflows/documentation/extract_documentation_workflow.md)
  - [`workflows/documentation/inline_documentation_workflow.md`](./workflows/documentation/inline_documentation_workflow.md)

## Quick usage guidance

- Python demo (small numerical scripts): use the Python Demo workflow. Keep code compact, vectorize with NumPy, and gate debug prints by verbosity.
- C/C++ prototyping and debugging: use the C++ Prototype or Debugging workflows. Build/run via project bash scripts (e.g., `C/build.sh` or `tests_bash/.../*.sh`), prefer `printf` for debug, and fail loudly.
- Python ↔ C/C++ via ctypes: use the Python ctypes workflow.
  - Typical structure in this repo:
    - Python wrapper: `python/pyMolecular/eFF.py`
    - Auto-build + loader: `python/pyMeta/cpp_utils.py`
    - C ABI bridge (extern "C"): `cpp/libs/Molecular/eFF_lib.cpp`
    - C++ implementation: `cpp/common/molecular/eFF.h`
  - Set `cpp_utils.BUILD_PATH` to the library build dir and call `cpp_utils.loadLib('<name>')` to auto-recompile and load. Define `argtypes`/`restype` and pass contiguous NumPy arrays.
- Documentation:
  - Extract human-readable docs from code: use the Extract Documentation workflow.
  - Insert concise inline docs (Doxygen or 1-line Python docstrings): use the Inline Documentation workflow. Use safe update markers to allow regeneration.

Notes:
- Follow the global rules for style, debugging-first, and performance.
- A legacy notes file exists as `wrokflows.md`; prefer the curated workflows above.
