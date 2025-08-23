# Workflow: C/C++ Prototype

## Purpose
Minimal, correct first implementation of a numerical algorithm in C/C++. Emphasis on clarity and testability.

## Inheritance
- Inherits `doc/CodingRules/global_rules.md`

## Respect bash scripts 
- run via project bash scripts (typically `run.sh`),  it will automatically:
  - recompile the program
  - setup paths to inputs and libraries
  - run the program with proper arguments
  - ouput stdout to logfile
- you can edit the `run.sh` or other bash script if needed
- DO NOT try to compile it yourself using `make`, it will not work

## Workflow Steps
1. Setup
   - Use project bash scripts to build/run (see `tests_bash/`); capture logs
   - Produce deterministic, text‑based output for inspection
2. Coding strategy
   - Short, flat functions; avoid premature abstraction
   - produce numerical C-like code, from C++ use mostly just templates and classes
   - Prefer raw C arrays (`double*`) in hot paths; small POD structs
   - prefer passing pointers over deep copying when possible
   - avoid exesive used of std::algorithms unless user ask otherwise
   - use std::unordered_map and std::unordered_set for associative lookups where array is not suitable
   - use stack alocated arrays when possible (unless too long)
3. Debugging strategy
   - Add `printf` at entries and key computations; assertions for invariants
   - Comment out alternatives; mark TODO/DEBUG blocks
   - calculate min,max of resulting arrays to spot NaN, Inf or other strange values
4. Testing & validation
   - Add end‑of‑run sanity summaries (mass/energy/area)
5. Performance considerations
   - Mind cache locality; 
   - avoid unnecessary allocations and deallocations inside hot paths
6. Visualization/reporting
   - Leave plotting to Python; write simple logs or CSV for consumption

## Notes
- For vector math use `Vec2.h`/`Vec3.h` and `quaternion.h` helpers already in repo
