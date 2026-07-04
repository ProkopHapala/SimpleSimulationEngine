We develop rigorous scientific software where debuggability, physical consistency, numerical correctness and stability, as well as performance and simplicity are paramount. Follow these principles:

## Core Principles

- **KISS** (Keep It Simple), Simplest solution that works. One-liner > ten-liner.
- **AHA** (Avoid Hasty Abstractions), avoid boilerplate
- **YAGNI** : **Surgical Edits** — Touch only what's needed. No unrelated cleanup. Comment out, don't delete. Ask if ambiguous.
- **DRY** : Inventory existing code before writing new. Generalize rather than duplicate. See skill:`code-reuse`.
- **SoC** (Separation of Concerns), separate module for Compute, plotting, Backend, CLI, GUI. Thin test scripts call general workhorse function from shared modules.
- **SSOT** : Authoritative single source of truth must be defined to avoid ambiguity and confusion
- **TDD** : Define verification before coding. Parity checks vs reference/analytical/physical invariants. Run tests after every change.
- **Fail Fast** : No silent fallbacks (try-catch). Crashes with stack traces > masked bugs. Look for root cause, not symptoms.
- Compact code, unlimited line lengh (function call must be one line).  Short names for math symbols (`E_tot`, `T_ij`).

## Never Do This

- Never delete or rearrange existing code without explicit permission
- Never perform unrelated aesthetic/style edits
- Never apply quick-fixes that hide root causes (e.g. hard-coded outputs)
- Never reinvent functionality already implemented. Inventory fist, check the provided examples and references, base classs etc.

## Debugging & Testing

Fail loud — crashes with stack traces > masked bugs. Debug prints gated by verbosity. Tests in `tests/` with ref_data regression. Numerical correctness via parity checks vs reference/analytical/physical invariants. See skill:`test-runner`, skill:`visual-debug`, skill:`gpu-debug`, skill:`numerical-parity`,

## Performance & Languages

* Minimize Python orchestration; push compute to C/C++/OpenCL kernels. Flat arrays, cache-aware, preallocate. See skill:`python-perf`, skill:`port-to-opencl`.
* ctypes auto-build, zero-copy buffer wrapping, Fortran/C memory layout. See skill:`ctypes-bindings`.
* C/C++ Build & Memory : Build via cmake with ASAN/opt symlink switching. Memory ownership: one owner per heap array, borrowed aliases, `bOwn*` flags. See skill:`cpp-build`, skill:`cpp-memory`.
* GPU/OpenCL : memory latency, gather over scatter, local memory, minimize host-device transfers. See skill:`gpu-optimize`.
* CPU/C++ : cache hierarchy, SIMD vectorization, data-oriented design, loop optimization. See skill:`cpu-perf`.

## Code Reusability & Shared Libraries

- **Placement test**: "Could another app use this?" → if yes, `cpp/common/` or `cpp/common_SDL/` (global include path); if no, `cpp/apps/<AppName>/`.
- **Never copy-paste** between apps — extract to shared lib and include. Search `common/` and `common_SDL/` before writing any new function.

## Documentation & Navigation

- Before writing: search existing implementations (skill:`doc-read-navigate`)
- After implementing: update README.md + topical audits (skill:`doc-task-summary`)
- Dedicated doc work: OKF format, extract/inline workflows (skill:`doc-audit`)
- `doc/TopicalAudit/` — cross-language implementation maps
- `doc/AGENTs/skills/` — all skills; `doc/AGENTs/protocols/` — domain protocols
- `README.md` per folder — local index; `CODEMAP.md` — repo structure
- Visualization: separate compute from plotting; `plt.show()` only in CLI/main
