# user task:

create me description of coding rules and work flow for windsurf

read this tutorial how to do it before you start:

https://www.obviousworks.ch/en/vibe-coding-with-windsurf-ide/

also read some broder backround here:
https://blog.jetbrains.com/idea/2025/05/coding-guidelines-for-your-ai-agents/
https://devin.ai/agents101#introduction
https://www.anthropic.com/engineering/building-effective-agents

Now I want to crerate workflow specific for sketching demo how to solve some numerical scientific algorithm ( typically form domain physics, chemistry, computational geometry, graphics ). I haavily use numpy and matplotlib, I care about performance, 

there are some rules I already wrote

```Markdown
# Python numerical coding rules (NumPy/Matplotlib)

Concise rules tailored for developing numerical algorithms in python.

## Core principles
- Minimal dependencies: standard library + NumPy/Matplotlib, unless user ask otherwise
- Keep code concise and minimal; prefer small, testable changes, step-by-step, divide-and-conquer
- Reuse existing functions; avoid duplication; 
   - check if usefull functions exist and report if they need modification
- Prefer pure, data-oriented functions with explicit inputs/outputs
   - use default named arguments to avoit exesively long call-strings
- Fail loudly: avoid broad try/except; use assertions for invariants
- Comment out deprecated code when experimenting (do not delete it); mark with TODO/DEBUG

## Performance (NumPy-first)
- Avoid Python loops and branches in hot paths, they are slow
   - prefer pre-allocated numpy arrays over lists
   - prefer using numpy boolen array masks and advanced indexing over if branches 
   - if you need grow datastructute use list/dict/set comprehension
- prefer tuples and data-classes over dictionaries
- be mindful of memory layout and dtype (default float64, for pyOpenCL use float32)
- prefer importain library at the start of script/module rather then in functions

## Debugging/logging, test-driven development, (debuggability > UX)
- When you desing algorithm thing about how to debug it and test it from the start
  - for LLM we need information in terminal output in text-form
    - make sure debug prints are easy to read not exesively long, like a table
  - for user the plots and visializations are preferable
  - Add quick assertions and end-of-run sanity checks (area/mass/energy)
- at initial stages of development add lot of debug-prints which help follow the algorithm and spot errors
   - be carefull adding too many debug prints inside hot loops  
- Gate all debug prints by a verbosity level; keep prints as one-liners
- If multi-line diagnostics are needed, wrap them in a small helper (to avoid bloating the code with debug-prints)
- Do not store debug-only data in core objects by default
- Prefer simple prints to heavy logging frameworks
- for key numpy arrays print overall debug information (e.g. min,max,shape) which help use spot strange values (nan, inf, etc.) 

## Plotting (Matplotlib)
- Separate compute and visualization; no plotting in core algorithms
- Plotting should be optional: add flags like `--noPlot`, `--saveFig`
- Call `plt.show()` only in CLI/main, never in library code
- prefer to reuse general plotting functions like `plot1d`, `plot2d` from `plot_utils.py` or other such modules present in the repo
- make plotting code concise
   - unless asked otrwise, fancy titples, labels and formating is not necessary, do not spend many lines of code on that 

## API/CLI
- Add `if __name__ == "__main__":` guard to CLI code instead of `main()` function
- Use argparse for CLI arguments with most global parameters unless too long 
- Provide `--verbosity` to control diagnostics; default to quiet

## Style
- Short, clear variable names, idealy representing using physical or mathematicall symbols (`E` for Energy, `T` for temparature etc.) 
- Compact code, avoiding empty lines, prefering one-liners
- Avoid line wrapping when it harms readability of expressions (assume wide editor)
   - at the start of all files set numpy print options allow infinite line-length `np.set_printoptions(linewidth=np.inf)` 
- Inline comments behind code line 
   - reserve header comments only for high-level structures (functions, classes, modules)
```

There are some more

```Marksdown


## NEVER DO THIS

- Never delete, overide or rearange existing code without explicit permission (unless it is explicitly requested in the user prompt)
  - if you think it should be removed comment it out
- Never do random edits of existing code just for esthetic reasons (e.g. to stick to coding style or other good practices)
- Never try to fix problem by patching the code by some "quick-fix" which does not address the root cause. Something like hard-coding expected results to pass the test is strictly forbidden.

## Debugging first

- Focus on debuggability, not user experience => do not try to hide issues.  
- Add debug prints to track variables and flow (we remove them later).  
- Never silently handle errors => fail loudly.  
- Avoid try-except; crashes with stack trace are helpful.  
- Comment out old code instead of deleting it.  
- Make small, testable changes and run often.  
- Log function entries and key conditions when helpful.  
- Mark experimental or unfinished code clearly (e.g. # TODO, # DEBUG).

## Guidelines 

- Be concise when editing code. Focus on the task and do minimum changes in the code necessary to solve the task.
- When writing documentation, notes or plans be rather detailed and concrete, citing names of files, functions, variables, etc.
- Check in the existing code first if there are already usefull functions implemnting similar task and try to reuse them to minimize code duplication
- If you have any doubts, better ask how to proceed
- Strive to write elegant, concise (compact) and moduar code, reuse existing functions and classes, avoid code duplication
- Strive to write performant data oriented code that is cache friendly and avoiding overheads - we are oriented towards performance scientific high performance computing 

## Style

- prefer concise and compact code.
- short code is easier to read and understand and see higher level structure.
- avoid bloated and unnecessarily complicatied code
- prefer short variable names, even single letter from math and physics (e.g. E for energy, m for mass, etc.)
- prefer one-liner expressions.
- prefer explanatory comments behind the code line (i.e. inline comments), rather then above
- do not break lines, especially in function calls (assume infinite line-width)
- For Doxygen documentation use `///` comment instead of `/*`, avoid using `/*` and `*/` for any comments altogther
- in C++ we use `printf` instead of `std::cout` for printing and debugging.
- in C++ we use plain C arrays like `double*` instead of `std::vector` especially when passing as function arguments. 
- for vector math w in C/C++ we use `Vec3.h` and `Vec2.h` headers with `Vec3d` and `Vec2d` structs defining common operators like +, -, *, /,  and functions like dot(), norm(),cross(), etc.
```

Please restructure it to the recomented format

We should think how to sparate it between general rules (global rules)
and task-specific rules / guidelines / workflows

Generally I can imagine to have several workflow

1) python demo workflow
2) python debugging workflow
3) C/C++ initial implementation workflow (minimal prototype)
4) C/C++ debugging 
 - Note for C/C++ I typically use bash-scripts which directly compile and run the program, and produce log file. Agent should be aware of this.
5) Python which calls C/C++ library over ctypes interface (I need to keep track / synchronize changes in the C/C++ library, with the update of python API, and update the test-script)
6) PyOpenCL workflow. Here I need agent to care about both the kenrells in .cl files and the python interface and the test scripts. I have my own python class OpenCLBase  which builds upon pyOpenCL and and some more convenience functions/abstractions like crating buffest from dictionary of buffer names and sizes, parsing kernel headers from .cl file, metaprograming where certain makro-markets in .cl file are textually replaced by some pieces of code from liberary etc. I will provide you more details later. 

Now do not go into all details, give me recomandation for the overall structure, where you left some visible markse where we need to fulfill details later. The main thin is that you shoud recomand me how many files / rules / workflow to clreate with what contnet, how to separate it  



---

# LLM answer: Google Gemini 2.5 Pro

Based on my analysis of the provided articles and your existing rules, I recommend structuring your AI agent's instructions using the "Vibe Coding" methodology. This approach treats the AI as a junior partner that needs a clear, consistent "vibe" or context to perform effectively. We will separate this into three main parts:

1.  **Global Rules:** A high-level "constitution" that governs the agent's behavior on all projects.
2.  **Project-Specific Rules:** A detailed file for each project that specifies the tech stack, architecture, and language-specific conventions.
3.  **Development Workflows:** Step-by-step guides for common tasks, ensuring the agent follows a predictable and effective process.

Here is a proposed restructuring of your rules into this format.

***

## Recommended Structure for AI Coding Rules & Workflows

### Part 1: Global Rules (`global_rules.md`)

This file acts as the universal set of principles for the AI agent. It should be loaded for every task to ensure core behaviors are always followed.

```markdown
# WindSurf Global Coding Rules

This document defines the core principles and absolute constraints for all development tasks.

## 1. Core Philosophy: Debugging First

The primary goal is to produce code that is easy to understand, test, and debug.

- **[DF] Debuggability > User Experience:** Do not hide issues. The top priority is providing clear information for debugging.
- **[FL] Fail Loudly:** Never silently handle errors with broad `try-except` blocks. Crashes with a full stack trace are valuable for debugging.
- **[AC] Atomic Changes:** Make small, self-contained, and testable changes. Run tests and checks often.
- **[CDD] Comment, Don't Delete:** When experimenting or refactoring, comment out old code instead of deleting it. Mark it with `// TODO:` or `// DEBUG:` to indicate its status.

## 2. Absolute Constraints (NEVER DO THIS)

These are strict rules that must never be broken.

- **[NOC] No Unauthorized Code Changes:** Never delete, override, or rearrange existing code without explicit permission in the prompt.
- **[NEA] No Aesthetic Edits:** Do not perform random edits for purely stylistic reasons. Focus on the assigned task.
- **[NQF] No Quick-Fixes:** Never patch a problem with a "quick-fix" that doesn't address the root cause. Hard-coding results to pass a test is strictly forbidden.

## 3. General Guidelines

- **[C] Concise Code:** Keep code concise and minimal. Focus on the task and implement the minimum changes necessary.
- **[RF] Reuse Functions:** Before writing new code, check for existing functions that perform a similar task. Prioritize reuse to minimize duplication.
- **[DP] Detailed Planning:** When asked for a plan, be detailed and concrete. Cite specific file names, functions, and variables.
- **[AI] Ask if in Doubt:** If requirements are ambiguous, ask for clarification before proceeding.
```

### Part 2: Project-Specific Rules (Template for `.windsurfrules`)

This file provides the specific context for a given project, including the tech stack, directory structure, and language-specific conventions.

```markdown
# Project Rules for: [PROJECT NAME]

## 2.1. Project Context

- **Project Description:** // TODO: A brief one-sentence description of the project's goal.
- **Tech Stack:**
    - **Primary:** Python 3.10, NumPy, Matplotlib
    - **Allowed:** SciPy, pandas (only if requested)
    - **Forbidden:** // TODO: List any libraries to explicitly avoid.
- **Directory Structure:**
    - ` /src`: Core algorithm implementations.
    - ` /tests`: Test scripts.
    - ` /scripts`: Main executable scripts.
    - ` /data`: Input/output data.

## 2.2. Python-Specific Rules (Numerical Science)

### Performance (NumPy-first)
- Avoid Python loops in hot paths; prefer NumPy vectorization.
- Use pre-allocated NumPy arrays over appending to Python lists.
- Use boolean array masks for conditional logic instead of `if` branches in loops.
- Be mindful of memory layout and `dtype` (default `float64`).
- Import libraries at the start of the module, not inside functions.

### Debugging & Logging
- Gate all debug prints by a verbosity level (e.g., `if verbosity > 1:`).
- Keep debug prints as concise one-liners. For key NumPy arrays, print shape, min, and max to spot `nan`/`inf` values.
- Use assertions for invariants and end-of-run sanity checks.

### Plotting (Matplotlib)
- Separate computation from visualization. Core algorithm functions should not create plots.
- Make plotting optional with flags like `--noPlot` or `--saveFig`.
- Call `plt.show()` only from within an `if __name__ == "__main__"` block.
- Keep plotting code simple; fancy titles and labels are not necessary unless requested.

### API/CLI
- Guard all executable code with `if __name__ == "__main__":`.
- Use `argparse` for command-line arguments.
- Provide a `--verbosity` argument to control diagnostic output.

### Style
- Use short, clear variable names, often matching mathematical symbols (e.g., `E` for Energy).
- Write compact code, avoiding excessive empty lines.
- Do not wrap long lines, especially for mathematical expressions (assume a wide editor). Set `np.set_printoptions(linewidth=np.inf)` at the start of scripts.
- Use inline comments (`#`) placed after the code on the same line.

## 2.3. C/C++ Specific Rules (If applicable)

- Use `printf` for debugging instead of `std::cout`.
- Use plain C-style arrays (e.g., `double*`) for function arguments instead of `std::vector`.
- For vector math, use the `Vec3d` and `Vec2d` structs from `Vec3.h` and `Vec2.h`.
- Use `///` for Doxygen documentation; avoid `/* ... */` style comments.

## 2.4. // TODO: PyOpenCL-Specific Rules
// Details regarding the OpenCLBase class, kernel metaprogramming, and buffer management will be added here.
```

### Part 3: Development Workflows

This section outlines the high-level, step-by-step processes for different types of tasks. This provides the agent with a predictable plan of action.

```markdown
# Development Workflows

## 3.1. Workflow: Python Demo Sketching

**Goal:** Create a self-contained Python script to demonstrate a scientific algorithm.

1.  **Plan:** Prompt the agent to outline the core mathematical steps, required inputs, and expected outputs (data, plots).
2.  **Scaffold:** Create a single Python file. Add imports, `np.set_printoptions`, and an `if __name__ == "__main__"` guard.
3.  **Implement Core Logic:** Write the numerical algorithm using NumPy. Add assertions to check array shapes and value ranges at critical points.
4.  **Add Visualization:** Implement a separate function for plotting. Call this function from the main block, controlled by a command-line flag.
5.  **Add CLI:** Use `argparse` to define parameters like input data paths, algorithm settings, and verbosity.
6.  **Refine & Verify:** Run the script. Add debug prints to trace the algorithm's state and verify the results.

## 3.2. Workflow: Python Debugging

**Goal:** Identify and fix a bug in an existing Python script.

- // TODO: Define the systematic process for debugging.
- // Example steps could include: 1. Replicate the error. 2. Add targeted debug prints. 3. Analyze the output to isolate the fault. 4. Propose a fix. 5. Verify the fix.

## 3.3. Workflow: C/C++ Initial Implementation

**Goal:** Create a minimal, compilable C/C++ prototype.

- // TODO: Define the workflow for C/C++ development.
- // Note: Agent must be aware of the bash script used for compilation and execution (e.g., `compile_and_run.sh`). The workflow should include modifying the C++ source, running the script, and checking the log output.

## 3.4. Workflow: C/C++ Debugging

**Goal:** Debug a C/C++ program using logs and compilation scripts.

- // TODO: Define the C/C++ debugging cycle, emphasizing interaction with bash scripts and parsing text-based log files.

## 3.5. Workflow: Python with CTypes Interface

**Goal:** Update a C/C++ library and its Python `ctypes` wrapper.

- // TODO: Define the workflow for keeping the C/C++ library, Python API, and test scripts synchronized.
- // It should track changes in the C/C++ function signatures and update the Python `argtypes` and `restype` accordingly.

## 3.6. Workflow: PyOpenCL Development

**Goal:** Develop and test a PyOpenCL kernel and its host code.

- // TODO: Detail the specific workflow for PyOpenCL.
- // This should include steps for: 1. Modifying the `.cl` kernel file. 2. Updating the Python host code that uses the `OpenCLBase` class. 3. Handling the metaprogramming replacement markers. 4. Running the Python test script to verify kernel behavior.
```


# LLM answer: GPT-5

---

Got it — you want to evolve your existing rules into a **structured, Windsurf-compatible workflow & coding guidelines system**, organized so that (a) the global rules are clear, (b) task-specific workflows are modular, and (c) it’s easy for you (and an AI agent) to pick the right workflow depending on context.

Here’s a **recommended overall structure** (with placeholders where details can be filled later). I’m basing it both on the Windsurf tutorial style and the background agent guidelines.

---

# 📁 Proposed Structure

### 1. **Global Coding Rules**

File: `global_rules.md`
Purpose: Applies to all languages (Python, C/C++, OpenCL), all workflows. This is the “constitution” — the strict rules no workflow may break.

Content sections:

* **Core principles** (minimal dependencies, concise code, explicit inputs/outputs, fail loudly, comment-out vs delete, etc.)
* **Never do this** (no silent quick-fixes, no random edits, no deleting)
* **Debugging first** (crash > silent fail, debug prints gated by verbosity, etc.)
* **Style (cross-language)** (short names, inline comments, compactness, wide editor assumption, Doxygen conventions, printf vs std::cout, arrays vs std::vector, Vec2/Vec3 etc.)

👉 This file becomes the *shared root reference*. Every workflow inherits from it.

---

### 2. **Python Workflows**

Folder: `/workflows/python/`

Each workflow is one file. Suggested:

* `python_demo_workflow.md`
  For quick algorithm sketching (NumPy/Matplotlib). Minimal, fast iteration, clarity for presentation.

* `python_debugging_workflow.md`
  For bug-hunting or sanity checks. Verbose debug prints, heavy use of assertions, plots for diagnostics, gated verbosity.

* `python_ctypes_workflow.md`
  For interfacing with C/C++. Emphasis on synchronization between Python API, test scripts, and C/C++ library changes.

* `python_pyopencl_workflow.md`
  For GPU development. Covers `.cl` kernels, `OpenCLBase` helpers, buffers, macro substitutions. (To be detailed later when you provide more specifics.)

---

### 3. **C/C++ Workflows**

Folder: `/workflows/cpp/`

* `cpp_prototype_workflow.md`
  For minimal initial implementation (direct arrays, small test harness, bash build scripts).

* `cpp_debugging_workflow.md`
  For performance profiling and bug isolation. Includes: `printf` debugging, log file conventions, bash automation, sanity checks.

---

### 4. **Workflow Templates**

File: `/workflow_template.md`
Reusable schema for new workflows. Each workflow file should follow the same format, e.g.:

```markdown
# Workflow: [Name]

## Purpose
[When to use this workflow, what stage of development]

## Inheritance
- Inherits all [Global Rules]
- Overrides / adds rules for this context

## Workflow Steps
1. Setup (imports, build commands, data files, etc.)
2. Coding strategy (structure, language-specific constraints)
3. Debugging strategy
4. Testing & validation
5. Performance considerations
6. Visualization / reporting

## Notes
[Edge cases, special warnings, known pitfalls]
```

---

### 5. **Utilities**

Folder: `/utils/`

* `plot_utils.py` (shared plotting helpers)
* `debug_utils.py` (print helpers, array sanity checks, verbosity gates)
* Possibly: `test_utils.py` (sanity assertions, test harness)

---

# 🗂️ Overview of Files

```
coding_rules/
│
├── global_rules.md
├── workflow_template.md
│
├── workflows/
│   ├── python/
│   │   ├── python_demo_workflow.md
│   │   ├── python_debugging_workflow.md
│   │   ├── python_ctypes_workflow.md
│   │   └── python_pyopencl_workflow.md
│   │
│   └── cpp/
│       ├── cpp_prototype_workflow.md
│       └── cpp_debugging_workflow.md
│
└── utils/
    ├── plot_utils.py
    ├── debug_utils.py
    └── test_utils.py
```

---

# ✅ Why this structure?

* **Separation of concerns**: Global rules are central, workflows are contextual.
* **Agent-friendly**: Clear mapping: “I’m doing X → I open Y workflow.”
* **Human-friendly**: You don’t have to re-read the whole constitution when working — just open the relevant workflow.
* **Extensible**: You can later add e.g. “workflow for HPC cluster runs” or “workflow for visualization notebooks” without changing the base rules.
* **Reuses Windsurf conventions**: Each workflow = modular “vibe coding environment.”

---

Would you like me to **draft the `global_rules.md` and one example workflow (`python_demo_workflow.md`)** in this structure, so you can see the tone/format before we expand all others?


---

Perfect. Let’s fill in **all the workflow files** using the structured format we discussed, while making them compact, practical, and agent-friendly. I’ll base them on your existing rules and add context-specific adjustments.

---

# 📄 `global_rules.md`

```markdown
# Global Coding Rules

Applies to all workflows (Python, C/C++, OpenCL). These are the non-negotiable baseline rules.

## Core principles
- Minimal dependencies: standard library + NumPy/Matplotlib unless explicitly requested.
- Concise, modular code: prefer small, testable changes; divide-and-conquer.
- Reuse existing functions; avoid duplication.
- Pure, data-oriented functions with explicit inputs/outputs.
- Fail loudly: no silent handling; assertions for invariants.
- Comment out deprecated/experimental code instead of deleting; mark with TODO/DEBUG.

## NEVER DO THIS
- Never delete or override existing code without explicit permission.
- Never apply random edits for aesthetics or style.
- Never quick-fix by hardcoding results; always address root cause.

## Debugging First
- Debuggability > user experience. Do not hide issues.
- Add debug prints for key variables and flow (remove later).
- No broad try/except; crashes with stack traces are informative.
- Mark unfinished or experimental code clearly.
- Make small, testable changes; run often.
- Keep debug prints concise, preferably one-liners; use verbosity gates.

## Performance
- Prefer NumPy/vectorized ops (in Python).
- Avoid Python loops/branches in hot paths.
- Pre-allocate arrays; avoid dynamic growth in performance-critical sections.
- Be mindful of memory layout, dtype, cache-friendliness.
- In C/C++ prefer raw arrays (`double*`) over containers in hot paths.
- For vector math in C/C++ use `Vec2/Vec3` headers.

## Style
- Short, clear variable names (E for energy, T for temperature).
- Compact, minimal code. Avoid empty lines.
- Prefer inline comments after code lines; reserve headers for modules/classes/functions.
- Assume wide editor: avoid line wrapping in expressions.
- For Doxygen: use `///`, not `/* */`.
- In C++: prefer `printf` to `std::cout`.

## Visualization
- Separate compute vs plotting; no plotting in core algorithms.
- Plotting optional via flags (`--noPlot`, `--saveFig`).
- `plt.show()` only in CLI/main, never in library code.
- Reuse plotting utilities (`plot_utils.py`) where possible.
```

---

# 📄 `workflow_template.md`

```markdown
# Workflow: [Name]

## Purpose
When to use this workflow, what stage of development it supports.

## Inheritance
- Inherits all [Global Rules]
- Adds or overrides rules for this context

## Workflow Steps
1. Setup (imports, build commands, configs)
2. Coding strategy
3. Debugging strategy
4. Testing & validation
5. Performance considerations
6. Visualization/reporting

## Notes
Edge cases, warnings, known pitfalls.
```

---

# 📄 `/workflows/python/python_demo_workflow.md`

```markdown
# Workflow: Python Demo

## Purpose
Quick sketch of a scientific/numerical algorithm (physics, chemistry, geometry, graphics).
Optimized for clarity, minimal code, and simple visualization.

## Inheritance
- Inherits [Global Rules]
- Adds demo-oriented plotting and presentation rules

## Workflow Steps
1. Setup
   - Imports: `import numpy as np`, `import matplotlib.pyplot as plt`
   - Set NumPy print options: `np.set_printoptions(linewidth=np.inf)`

2. Coding Strategy
   - Keep functions small and self-contained.
   - Default args for key parameters (e.g. dt=0.01).
   - No plotting in core functions; only in `__main__`.

3. Debugging Strategy
   - Print shapes, min/max, `np.isnan` checks.
   - Assertions for physics invariants (energy, mass, etc.).

4. Testing & Validation
   - Compare against analytical solution or simpler case.
   - End-of-run sanity checks.

5. Performance
   - Use NumPy vectorization.
   - Avoid loops unless prototype requires clarity.

6. Visualization/Reporting
   - Provide 1D/2D plots using `plot_utils.py`.
   - Compact plot code; only basic labels.
   - `plt.show()` at the very end.

## Notes
- Do not overcomplicate formatting.
- Keep demo code under ~50–100 LOC.
```

---

# 📄 `/workflows/python/python_debugging_workflow.md`

```markdown
# Workflow: Python Debugging

## Purpose
Diagnosing incorrect behavior in scientific/numerical algorithms.

## Inheritance
- Inherits [Global Rules]
- Adds verbose debugging conventions

## Workflow Steps
1. Setup
   - Add `--verbosity` CLI flag.
   - Debug prints gated by verbosity.

2. Coding Strategy
   - Instrument functions with entry/exit prints.
   - Use helper functions in `debug_utils.py` for array stats (min, max, shape, NaN count).
   - Comment out old/experimental code, never delete.

3. Debugging Strategy
   - Print intermediate results (e.g. first few values of arrays).
   - Add quick plots for sanity checks (optional).
   - Assertions for invariants at every step.

4. Testing & Validation
   - Use reduced-size input to isolate bugs.
   - Add regression tests after bug is fixed.

5. Performance
   - Performance is secondary during debugging.
   - Prefer clarity, even loops, if it aids inspection.

6. Visualization/Reporting
   - Scatter/line plots for diagnostics.
   - Mark unusual values clearly (NaN, Inf).
```

---

# 📄 `/workflows/python/python_ctypes_workflow.md`

```markdown
# Workflow: Python + C/C++ (ctypes)

## Purpose
Glue layer between Python scripts and C/C++ libraries. Synchronize code changes across layers.

## Inheritance
- Inherits [Global Rules]
- Adds synchronization rules for mixed Python/C++

## Workflow Steps
1. Setup
   - Shared library built with bash script (e.g. `make.sh`).
   - Python loads via `ctypes.CDLL`.

2. Coding Strategy
   - Mirror function signatures in Python with `ctypes`.
   - Define argument/return types explicitly.
   - Keep wrapper functions minimal.

3. Debugging Strategy
   - On mismatch, print raw pointer values and array contents.
   - If crash: confirm array dtype/contiguity matches C++ side.

4. Testing & Validation
   - Always update Python test script when C++ signature changes.
   - Write small test harness comparing Python vs C++ outputs.

5. Performance
   - Pass NumPy arrays as raw pointers.
   - Avoid copying data between Python and C++ unnecessarily.

6. Visualization/Reporting
   - Use Python for plotting results of C++ computations.
   - Debug prints in C++ go to log file; Python parses/prints if needed.

## Notes
- Maintain sync: C++ header ↔ ctypes wrapper ↔ test script.
```

---

# 📄 `/workflows/python/python_pyopencl_workflow.md`

```markdown
# Workflow: Python + PyOpenCL

## Purpose
Numerical kernels running on GPU via `.cl` files, interfaced from Python. Handles buffer management and kernel compilation.

## Inheritance
- Inherits [Global Rules]
- Adds GPU-specific considerations

## Workflow Steps
1. Setup
   - Use `OpenCLBase` class for context, queue, and buffers.
   - Keep `.cl` kernels in `kernels/` folder.

2. Coding Strategy
   - Define kernel arguments clearly.
   - Parse kernel headers automatically where possible.
   - Use meta-programming macros for boilerplate expansions.

3. Debugging Strategy
   - Start with small problem sizes.
   - Check buffer contents with `cl.enqueue_copy`.
   - Print kernel build logs when compilation fails.

4. Testing & Validation
   - Cross-check GPU vs CPU NumPy reference implementation.
   - Add regression tests for correctness.

5. Performance
   - Prefer `float32` for GPU.
   - Align memory sizes to work-group size.
   - Minimize host-device transfers.

6. Visualization/Reporting
   - Use Python plotting for GPU results.
   - Plot timings and performance scaling.

## Notes
- Keep kernel code simple first; optimize later.
- Clearly mark macros and auto-generated code blocks.
```

---

# 📄 `/workflows/cpp/cpp_prototype_workflow.md`

```markdown
# Workflow: C/C++ Prototype

## Purpose
Minimal first-pass implementation of numerical algorithms in C/C++. Focus on clarity, not optimization.

## Inheritance
- Inherits [Global Rules]

## Workflow Steps
1. Setup
   - Use bash script to compile & run (`build.sh`).
   - All debug output goes to terminal or log file.

2. Coding Strategy
   - Use raw arrays (`double*`) for clarity.
   - Minimal structs/classes; only what is essential.
   - Prefer short, flat functions.

3. Debugging Strategy
   - Add `printf` at function entries and key computations.
   - Use assertions for invariants.
   - Comment out alternative implementations.

4. Testing & Validation
   - Compare outputs against Python prototype.
   - Sanity checks: conservation of mass/energy.

5. Performance
   - Initial version not optimized; correctness first.
   - Note hotspots for later profiling.

6. Visualization/Reporting
   - No plotting in C++; produce log files for Python plotting.
```

---

# 📄 `/workflows/cpp/cpp_debugging_workflow.md`

```markdown
# Workflow: C/C++ Debugging

## Purpose
Finding and fixing bugs in C/C++ numerical algorithms. Debugging performance and correctness.

## Inheritance
- Inherits [Global Rules]
- Adds debugging-specific practices

## Workflow Steps
1. Setup
   - Compile with `-g -O0` flags for debugging.
   - Use bash scripts to recompile and run tests automatically.

2. Coding Strategy
   - Keep changes small and localized.
   - Instrument with `printf` rather than debuggers when possible.

3. Debugging Strategy
   - Track memory allocations manually.
   - Print array sizes, key values, invariants.
   - Use log files for long runs.

4. Testing & Validation
   - Cross-check against Python prototype.
   - Write temporary "mini-tests" directly in `main()`.

5. Performance
   - Focus on correctness before optimizing.
   - Only optimize after bug-free.

6. Visualization/Reporting
   - Parse log files in Python; use matplotlib for diagnostics.
```

---

# 📁 `/utils/`

* `plot_utils.py` – reusable simple plotting (`plot1d`, `plot2d`, quick scatter).
* `debug_utils.py` – array sanity check (`min/max/shape/nans`), gated debug prints.
* `test_utils.py` – small harness for regression tests & sanity checks.

---

✅ This gives you **complete content for each workflow file + global rules + utilities plan**.
Would you like me to now **draft the `utils/plot_utils.py` and `debug_utils.py` skeletons** (minimal but reusable), so you have a working base?


---

