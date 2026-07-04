---
name: ctypes-bindings
description: Python↔C/C++ ctypes bindings — auto-build, zero-copy buffers, memory layout
trigger:
  glob:
    - "**/pyBall/**/*.py"
    - "**/cpp/**/*.cpp"
---

## Setup

- Prefer auto-build via `pyMeta.cpp_utils`: set `BUILD_PATH` to lib build dir, `lib = cpp_utils.loadLib('<name>')` (auto-recompile if needed).
- Load functions and set `argtypes`/`restype` explicitly.

## Memory Layout: Fortran vs C/C++/Python

**Fortran is column-major** (fastest index first): `rho[imu,inu,ineigh,iatom]`
**C/C++/Python is row-major** (fastest index last): `rho[iatom,ineigh,inu,imu]`

**Do not manually reorder nor use `order='F'` in NumPy.** Accept the natural memory layout for numerical efficiency. Pass raw pointers transparently without copying/rearranging. The interface should be transparent and overhead-free.

## Buffer Wrapping (C++ → NumPy zero-copy)

This is the core pattern used across FireCore. C++ holds memory; Python gets live views.

**C++ side:** Register named pointers in a map, expose via `extern "C"`:
```cpp
std::unordered_map<std::string,double*> buffers;
double* getBuff(const char* name){ return buffers[name]; }
```

**Python side:** Wrap raw pointer into `np.ctypeslib.as_array`:
```python
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = c_double_p
def getBuff(name, sh):
    ptr = lib.getBuff(name.encode('utf8'))
    return np.ctypeslib.as_array(ptr, shape=sh)

# Usage: live view into C++ memory
apos  = getBuff("apos",  (natoms,3))
fapos = getBuff("fapos", (natoms,3))
Es    = getBuff("Es",    (6,))      # energy terms
```

Arrays must be contiguous with matching dtype/shape. Ownership stays on C++ side.

## Existing Bindings Examples

**Switches (tri-state: -1 disable, 0 keep, +1 enable):**
```python
lib.setSwitchesUFF.argtypes = [c_int]*7
lib.setSwitchesUFF.restype  = None
def setSwitchesUFF(DoBond=0, DoAngle=0, DoDihedral=0, DoInversion=0, DoAssemble=0, SubtractBondNonBond=0, ClampNonBonded=0):
    return lib.setSwitchesUFF(DoBond, DoAngle, DoDihedral, DoInversion, DoAssemble, SubtractBondNonBond, ClampNonBonded)

# Disable nonbonded for parity testing
setSwitchesUFF(DoBond=+1, DoAngle=+1, SubtractBondNonBond=-1)
```

**Initialize and fetch buffers:**
```python
lib.init_buffers.argtypes = [c_bool]   # bUFF
lib.init_buffers.restype  = None
init_buffers(bUFF=True)

# Then read ndims to know shapes, then getBuff()
ndims = getIBuff("ndims", (10,))   # [_natoms, nbonds, nangles, ...]
natoms = ndims[0]
apos   = getBuff("apos", (natoms,3))
```

## Testing

- Maintain small Python harness comparing CPU NumPy vs C/C++ outputs.
- Update on signature change. Keep sync triad: C/C++ header ↔ Python wrapper ↔ test script.
