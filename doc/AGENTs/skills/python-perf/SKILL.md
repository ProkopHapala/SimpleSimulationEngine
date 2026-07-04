---
name: python-perf
description: Performant Python for scientific computing — vectorization, NumPy anti-patterns, preallocation
trigger:
  glob:
    - "**/*.py"
    - "**/tests/**/*"
    - "**/*bench*.py"
    - "**/*perf*.py"
---

## Core Principle

**Python is the harness, not the engine.**  Call NumPy or OpenCL kernels on as many items at once as possible to minimize harness overhead. Use advanced array slicing. NEVER write low-level hot loops in Python.

## Rules

### 1. Batch Everything
- Minimize number of loop iterations in python.
- Process entire arrays in single NumPy/OpenCL calls
- Minimize Python function call overhead
- One kernel call for 1M points > 1M kernel calls for 1 point each

### 2. Vectorized Operations Only
- **Forbidden:** `for i in range(len(array)): result[i] = array[i] * 2`
- **Required:** `result = array * 2`
- **Forbidden:** Nested loops over 2D/3D grids
- **Required:** `np.meshgrid`, slicing, broadcasting, or OpenCL kernels

### 3. Minimal Allocation
- Preallocate buffers, reuse in hot paths
- Use in-place operations (`array *= 2` vs `array = array * 2`)
- Avoid intermediate arrays where possible

### 4. Advanced Slicing
- `array[mask]` for conditional selection
- `array[:, :, 0]` for channel extraction
- `array.reshape(-1)` for flattening

## Common Anti-Patterns

### ❌ Grid Iteration
```python
# WRONG: 1M Python iterations
for i in range(nx):
    for j in range(ny):
        result[i,j] = grid[i,j] * weight
```

### ✅ Vectorized
```python
# CORRECT: Single NumPy operation
result = grid * weight
```

### ❌ Position Generation
```python
# WRONG: Nested loops + stack overhead
positions = []
for ix in range(nx):
    for iy in range(ny):
        positions.append([x[ix], y[iy], z])
positions = np.array(positions)
```

### ✅ Meshgrid (Simple)
```python
# CORRECT: 3D grid directly, minimum allocation
xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
# Use xx, yy, zz directly as 3D arrays
```

### ❌ Conditional Loop
```python
# WRONG: Python loop with condition
for i in range(len(data)):
    if data[i] > threshold:
        result[i] = data[i] * scale
```

### ✅ Boolean Masking
```python
# CORRECT: Vectorized
mask = data > threshold
result[mask] = data[mask] * scale
```

## When Python Loops ARE Acceptable

- Orchestration: file I/O, argument parsing
- Small constants: <100 iterations, not in hot path
- Debugging: temporary prints (remove before commit)

