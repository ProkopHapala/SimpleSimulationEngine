---
name: cpp-build
description: Building C++ (FireCore/SPFF) — build dirs, cmake targets, ASAN vs opt modes
---

# FireCore C++ build (SPFF / pyBall)

## Build directory (mandatory)

We always build in `$REPO_ROOT/cpp/Build/`, which is a **symlink**.

Switch it to point to either:

- `cpp/Build-opt/` — production / fast runs
- `cpp/Build-asan/` — AddressSanitizer instrumentation

```bash
cd $REPO_ROOT/cpp
ln -sfn Build-opt Build    # or Build-asan
```

**Python loads only through the symlink**, not `Build-opt` or `Build-asan` directly:

- `pyBall/SPFF.py` → `cpp/Build/libs/Molecular/libSPFF_lib.so`
- `pyBall/cpp_utils_.py` → `cpp/Build/libs/`

Optional override: `export CPP_BUILD_PATH=/abs/path/to/cpp/Build/libs/`

Prefer `tests/tSPFF/run.sh` for SPFF tests — it recompiles fresh code and (when used for ASAN) sets sanitizer library paths. **Ensure the script’s build tree matches the run mode** (see ASAN section below).

## Compile SPFF_lib

After switching `Build` symlink, rebuild **in that tree** (mixed `.o` files across trees cause ASAN symbol errors or silent corruption):

```bash
cd $REPO_ROOT/cpp/Build
cmake --build . --target DynamicOpt SPFF_lib
# or legacy per-target make (as in tests/tSPFF/run.sh):
# cd Build/libs/Molecular && make SPFF_lib
```

Verify opt build is not linked against ASAN objects:

```bash
nm -D $REPO_ROOT/cpp/Build/libs/Molecular/libSPFF_lib.so | grep -i asan
# must print nothing for Build-opt
```

## Run modes

| Mode | `Build` → | Before Python | `LD_PRELOAD` |
|------|-----------|---------------|--------------|
| **opt** | `Build-opt` | `unset LD_PRELOAD ASAN_OPTIONS LSAN_OPTIONS` | off |
| **asan** | `Build-asan` | rebuild full dep chain in Build-asan | `libasan.so` required |

ASAN run (Hessian/debug only — **no matplotlib** in same process):

```bash
export LD_PRELOAD=$(g++ -print-file-name=libasan.so)
export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1:symbolize=1
export LSAN_OPTIONS=detect_leaks=0
```

## Common failures

| Error | Cause | Fix |
|-------|-------|-----|
| `undefined symbol: __asan_*` | opt `.so` or mixed asan/opt `.o` | `ln -sfn Build-asan Build`, rebuild all targets; or `unset LD_PRELOAD` for opt |
| ASAN + matplotlib `ft2font` crash | interceptor conflict | Hessian scripts without `import matplotlib`; plot in separate process |
| Stale library | `Build` symlink changed but `.so` not rebuilt | `cmake --build cpp/Build --target SPFF_lib` |

## `tests/tSPFF/run.sh` caveat

The script currently hardcodes `cpp/Build-opt` for `make` but always exports `LD_PRELOAD=libasan`. For trustworthy ASAN runs: point `Build` → `Build-asan`, rebuild there, preload ASAN, and run **without** mixing opt artifacts. For production: `Build` → `Build-opt`, **unset** `LD_PRELOAD`.

## Minimal ASAN test (no matplotlib)

```bash
cd tests/tSPFF
bash run_asan_minimal.sh asan    # or opt
```

Uses `test_asan_minimal.py` (numpy only). Default `os._exit(0)` after PASS avoids teardown false negatives.

## Related docs

- Memory ownership / double-free: `cpp-memory` skill → `doc/dev_notes/SPFF/Memory_Ownership_and_Deallocation.md`
- ASAN debugging plan: `doc/Topics/FTIR_Nanocrystals/Debug_ASan.md`
- Phonon physics (separate): `doc/Topics/FTIR_Nanocrystals/Debug_negative_phonon_freqs.md`
