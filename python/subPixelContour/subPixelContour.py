import numpy as np
import matplotlib.pyplot as plt
import argparse

'''
In many situations we wanto to describe objects with sharp boundary using rather low-resolution texture (e.g. lakes or castle walls on the maps, colision objects on grid, text font or glyphts, area with one component in multi-component fluid simulation or continuum mechanics). This provides representation which is both efficient to evaluate (as position of relevant node points can be quickly found by indexing the grid) and memory efficient (thanks to low resolution).

To do so we can express the contour of the object as implicit surface of 2D function F(x,y)=0 and express it as linar combination of basis functions positioned in each point F(x,y)= sum_i * c_i * B(x-x_i,y-y_i). In this sum we consider just a few nearest neighto grid points p_i=(x_i,y_j) because our basis function B(x,y) hav finite support. E.g. bicubic or bilinear interpolation.

Now write code which find optimal coefficients to describe certain shape using bilinear intepolation. Lets do it just on small patch of 4x4 points. The reference function describing the shape should be provided as analytical binary function (returs 1 if we are in, and -1 if we are out). The function is sampled on bunch of points with arbitrary position (typically user will provide points which are near the boundary to describe it precisely). In the end we do not have to care how these points are distributed. We only solve the fitting probelm that our linear combination F(x,y) must reproduce F_ref(x,y) at all the points with minimum error. So it is just linar regression problem.

---

Implementation notes:

1) Do not make hot-loop in python. Use efficient numpy operations (arry slicing, advanced indexing). Instead of if else use numpy boolean mask for indexing the arryay especially in evalutation of basis function B(x,y) 

2) the sample point should be rather uniformaly distributed on regular recantugalr grid rather than randomly 3) I see you have some error there

3) for debug plot use imshow. On that imshow plot the sub-pixel sampling points (with color black-or-white do illustrate if they are in or not) with nearst interpolation. Also plot there the bigger pints for the position of basis points. There should be 3 plots - reference function, model function, and the error. Print the final error (sum of square error over the whole image)

4) Print the expansion coefs at node points so you can debug it. 

5) Importaint is to consider that we should not fit the whole function, but we are interested just in the zero-line F(x,y)=0 (and its neighborhood). To promote this we need to calculate the fitting error with clamping the model fuction F(x,y) to interval (-1,1). This m,akes evaluation of the error non-linear and I'm not sure if lienar albebra like np.linalg.lstsq will work in that case. So perhaps we need to evaluate derivatives and do gradinet descent (?)
'''

# unlimited line lengh for numpy print
np.set_printoptions(linewidth=np.inf)

# Global exponent for bilinear_power; set from argparse in main
P_POWER = 2.0


def points_with_radius(x, y, radius=1.2):
    """Return boolean mask for points (x,y) lying within given radius from origin (0,0)."""
    return (x*x + y*y) <= (radius*radius)

def reference_shape(x, y, center=(1.5, 1.5), radius=1.2):
    """Analytical binary circle: returns +1 inside, -1 outside."""
    inside = (x - center[0])**2 + (y - center[1])**2 < radius**2
    return np.where(inside, 1.0, -1.0)

def evaluate_bilinear_bases(points, basis_points):
    """
    Evaluates the bilinear basis functions for a set of points.
    
    Args:
        points (np.ndarray): Array of points to evaluate at, shape (n_points, 2).
        basis_points (np.ndarray): Array of basis point locations, shape (n_basis, 2).
        
    Returns:
        np.ndarray: Matrix of basis function values, shape (n_points, n_basis).
    """
    n_points = points.shape[0]
    n_basis = basis_points.shape[0]
    A = np.zeros((n_points, n_basis))
    for j, (bx, by) in enumerate(basis_points):
        dx = points[:, 0] - bx
        dy = points[:, 1] - by
        # Bilinear basis function B(dx, dy) = max(0, 1 - |dx|) * max(0, 1 - |dy|)
        B_vals = np.maximum(0, 1 - np.abs(dx)) * np.maximum(0, 1 - np.abs(dy))
        A[:, j] = B_vals
    return A

def cubic_bspline_1d(t):
    """Cubic B-spline basis B(t) with support |t|<2 (vectorized)."""
    a = np.abs(t)
    out = np.zeros_like(a)
    m1 = a < 1.0
    out[m1] = (2.0/3.0) - a[m1]**2 + 0.5*a[m1]**3
    m2 = (a >= 1.0) & (a < 2.0)
    out[m2] = (1.0/6.0) * (2.0 - a[m2])**3
    return out

def evaluate_bspline_bases(points, basis_points):
    """Evaluate 2D cubic B-spline bases (Cartesian product of 1D splines)."""
    dx = points[:, 0][:, None] - basis_points[None, :, 0]
    dy = points[:, 1][:, None] - basis_points[None, :, 1]
    return cubic_bspline_1d(dx) * cubic_bspline_1d(dy)

def evaluate_rbf2x2_bases(points, basis_points, h=None, kernel='r2_linear'):
    """Evaluate compact 2x2-cell radial bases on a grid using vectorized numpy.

    Kernels (support r<1):
      - 'r2_linear'   : B = max(0, 1 - r^2)
      - 'r2_quadratic': B = max(0, 1 - r^2)^2
      - 'wendland_c2' : B = (1 - r)^4_+ * (4 r + 1)  [C2 Wendland, 2D]

    r^2 = ((x-x_i)/h)^2 + ((y-y_j)/h)^2, with h the grid spacing (deduced if None).
    Only four corners of the containing cell contribute.
    """
    xs = np.unique(basis_points[:, 0]); ys = np.unique(basis_points[:, 1])
    nx, ny = xs.size, ys.size
    if h is None:
        hx = np.min(np.diff(xs)) if nx > 1 else 1.0
        hy = np.min(np.diff(ys)) if ny > 1 else 1.0
        h = min(hx, hy)

    n = points.shape[0]
    A = np.zeros((n, basis_points.shape[0]))

    # integer cell corner indices for each point (floors)
    i0 = np.floor(points[:, 0]).astype(int)
    j0 = np.floor(points[:, 1]).astype(int)

    # 4 offsets for the cell corners, vectorized shape (4, n)
    ox = np.array([0,1,0,1])[:, None]
    oy = np.array([0,0,1,1])[:, None]
    ii = i0[None, :] + ox
    jj = j0[None, :] + oy

    # valid mask inside grid
    m = (ii >= 0) & (ii < nx) & (jj >= 0) & (jj < ny)
    if not np.any(m):
        return A

    # columns (basis indices) and rows (point indices) for scattered assignment
    cols = (jj * nx + ii)
    rows = np.broadcast_to(np.arange(n), ii.shape)

    # basis node coordinates for those (ii,jj); clip to stay index-safe, we'll mask later
    ii_s = np.clip(ii, 0, nx-1)
    jj_s = np.clip(jj, 0, ny-1)
    bx = xs[ii_s]
    by = ys[jj_s]

    # distances (4,n)
    dx = points[None, :, 0] - bx
    dy = points[None, :, 1] - by
    r2 = (dx / h)**2 + (dy / h)**2

    # evaluate kernel
    if kernel == 'r2_linear':
        vals = np.maximum(0.0, 1.0 - r2)
    elif kernel == 'r2_quadratic':
        vals = np.maximum(0.0, 1.0 - r2)
        vals = vals * vals
    elif kernel in ('wendland', 'wendland_c2'):
        r = np.sqrt(r2)
        t = np.maximum(0.0, 1.0 - r)
        vals = (t**4) * (4.0*r + 1.0)
    else:
        raise ValueError(f"Unknown rbf2x2 kernel: {kernel}")

    # scatter-add into A using mask
    A[rows[m], cols[m]] = vals[m]

    return A

def solve_with_gradient_descent(A, b, learning_rate=0.1, n_iterations=200, verbose=True, log=None, ftol=None):
    """
    Solves for the coefficients using gradient descent with clamping.
    """
    n_basis = A.shape[1]
    # Initialize coefficients to zero
    coeffs = np.zeros(n_basis)
    for i in range(n_iterations):
        F          = A @ coeffs
        F_clamped  = np.clip(F, -1, 1)
        error      = F_clamped - b
        d_clamp_dF = np.where((F > -1) & (F < 1), 1, 0)
        gradient   = A.T @ (d_clamp_dF * error)
        coeffs    -= learning_rate * gradient
        # logging
        if log is not None:
            log.setdefault('cost', []).append(float(np.sum(error*error)))
            log.setdefault('gnorm_c', []).append(float(np.linalg.norm(gradient)))
        # Optional: Print cost for debugging
        if verbose and (i % 10 == 0):
            cost = np.sum(error**2)
            print(f"Iteration {i}, Cost: {cost}")
        # early stop on force
        if ftol is not None:
            if np.linalg.norm(gradient) < ftol:
                if verbose:
                    print(f"[GD] Early stop at iter {i}: ||grad||={np.linalg.norm(gradient):.3e} < ftol={ftol}")
                break
            
    return coeffs

def solve_with_damped_dynamics(A, b, learning_rate=0.1, mu=0.9, n_iterations=200, verbose=True, log=None, ftol=None):
    """Damped-dynamics (with inertia) solver for linear model F = A @ c, with clamped error.

    v is the velocity in coefficient space. Force is -grad. If force opposes velocity (f·v < 0), reset v=0.
    """
    n = A.shape[1]
    c = np.zeros(n)
    v = np.zeros_like(c)
    for i in range(n_iterations):
        F = A @ c
        F_c = np.clip(F, -1.0, 1.0)
        e = F_c - b
        dphi = ((F > -1.0) & (F < 1.0)).astype(float)
        g = A.T @ (dphi * e)
        f = -g
        if np.dot(f, v) < 0.0:
            v[:] = 0.0
        v = mu * v + learning_rate * f
        c += v
        # logging
        if log is not None:
            log.setdefault('cost', []).append(float(np.sum(e*e)))
            log.setdefault('gnorm_c', []).append(float(np.linalg.norm(g)))
        if verbose and (i % 10 == 0):
            print(f"[Dyn] Iteration {i}, Cost: {np.sum(e*e):.6f}")
        # early stop
        if ftol is not None:
            if np.linalg.norm(g) < ftol:
                if verbose:
                    print(f"[Dyn] Early stop at iter {i}: ||grad||={np.linalg.norm(g):.3e} < ftol={ftol}")
                break
    return c

def predict_rational(A, c, w=None, clamp=False, eps=1e-12):
    """Predict F for a rational (normalized) linear combination.

    Definition (fixed weights w):
      Let A[p,i] = B_i(x_p) be basis samples. Define A_w = A diag(w), with w_i>=0 (default w_i=1).
      Numerator:   N_p = sum_i c_i (A_w)[p,i] = (A @ (w*c))_p
      Denominator: D_p = sum_i       (A_w)[p,i] = (A @ w)_p
      Prediction:  F_p = N_p / D_p

    Notes:
    - Locality: with compact support bases A is sparse per-row; D_p = A[p,:] @ w reduces to the local 4 neighbors for 2x2 support automatically.
    - Efficiency: compute N = A @ (w*c) and D = A @ w without forming A_w explicitly.
    - Clamp: optional clamping to [-1,1] is for diagnostics; do not use clamp here if you need the analytical derivative of F.
    """
    n = A.shape[0]
    if w is None: w = np.ones(A.shape[1])
    wc = w * c
    N = A @ wc
    D = A @ w
    F = N / np.maximum(D, eps)
    if clamp:
        F = np.clip(F, -1.0, 1.0)
    return F, N, D

def grad_rational_coeffs(A, c, b, w=None, clamp=True, eps=1e-12):
    """Gradient of squared error for rational model w.r.t. coefficients c.

    Model:
      F_p = N_p / D_p,  where  N_p = sum_i c_i w_i B_i(x_p) = (A @ (w*c))_p,  D_p = sum_i w_i B_i(x_p) = (A @ w)_p.
      Use cost J = 1/2 sum_p (phi(F_p) - b_p)^2, with phi the clamp to [-1,1] (optional).

    Derivation (w fixed, optimizing only c):
      ∂F_p/∂c_k = (w_k B_k(x_p)) / D_p = A[p,k] w_k / D_p.
      Let e_p = phi(F_p) - b_p and dphi_p = d phi(F_p) / dF_p = 1 if F_p in (-1,1), else 0.
      Then ∂J/∂c_k = sum_p e_p dphi_p ∂F_p/∂c_k = sum_p e_p dphi_p (A[p,k] w_k / D_p).
      Vector form: define t_p = (e_p dphi_p) / D_p; then g_k = w_k (A[:,k]^T t).
      Hence g = (A^T t) ∘ w, where ∘ is Hadamard product.

    Implementation:
      F, N, D = predict_rational(A,c,w,clamp=False). If clamp is True, use phi=clip and dphi=1 on (-1,1), 0 outside; else phi(F)=F and dphi=1.
      t = (e * dphi) / max(D,eps); g = (A.T @ t) * w; return g.
    """
    if w is None: w = np.ones(A.shape[1])
    F, N, D = predict_rational(A, c, w=w, clamp=False, eps=eps)
    if clamp:
        F_c = np.clip(F, -1.0, 1.0)
        e = F_c - b
        dphi = ((F > -1.0) & (F < 1.0)).astype(float)
    else:
        e = F - b
        dphi = 1.0
    t = (e * dphi) / np.maximum(D, eps)
    g = (A.T @ t) * w
    return g

def grad_rational_weights(A, c, b, w=None, clamp=True, eps=1e-12):
    """Gradient of squared error for rational model w.r.t. non-uniform weights w.

    Model: F = N/D, N = A @ (w∘c), D = A @ w.
    Derivative: ∂F_p/∂w_k = A[p,k] (c_k - F_p) / D_p.
    With e = phi(F) - b, dphi = 1 on (-1,1) else 0 when clamp=True.
    Thus: ∂J/∂w_k = ∑_p e_p dphi_p A[p,k] (c_k - F_p) / D_p.
    Vector form: let t = (e∘dphi)/D, then
      g_w = (A^T @ t) ∘ c - A^T @ (t ∘ F)
    """
    if w is None: w = np.ones(A.shape[1])
    F, N, D = predict_rational(A, c, w=w, clamp=False, eps=eps)
    if clamp:
        F_c = np.clip(F, -1.0, 1.0)
        e = F_c - b
        dphi = ((F > -1.0) & (F < 1.0)).astype(float)
    else:
        e = F - b
        dphi = 1.0
    t = (e * dphi) / np.maximum(D, eps)
    At_t = A.T @ t
    At_tF = A.T @ (t * F)
    g_w = (At_t * c) - At_tF
    return g_w

def solve_with_gradient_descent_rational(A, b, w=None, learning_rate=0.1, n_iterations=200, verbose=True, eps=1e-12, train_w=False, w_learning_rate=None, w_clip_min=1e-8, w_clip_max=None, log=None, ftol=None, wftol=None):
    """Gradient descent solver for rational model F = (A@(w*c)) / (A@w).

    Notes:
    - Keeps existing clamped-error strategy for robustness near the target contour.
    - If you set w=None it reduces to uniform normalization with D_p = sum_i B_i(x_p).
    - If train_w=True, simultaneously optimize non-uniform weights w (initialized to 1).
    """
    if w is None: w = np.ones(A.shape[1])
    if w_learning_rate is None: w_learning_rate = 0.5 * learning_rate
    c = np.zeros(A.shape[1])
    for i in range(n_iterations):
        # gradients at current state
        g = grad_rational_coeffs(A, c, b, w=w, clamp=True, eps=eps)
        if train_w:
            gw = grad_rational_weights(A, c, b, w=w, clamp=True, eps=eps)
        # logging and early-stop checks based on current gradients
        F, _, _ = predict_rational(A, c, w=w, clamp=False, eps=eps)
        e = np.clip(F, -1.0, 1.0) - b
        if log is not None:
            log.setdefault('cost', []).append(float(np.sum(e*e)))
            log.setdefault('gnorm_c', []).append(float(np.linalg.norm(g)))
            if train_w:
                log.setdefault('gnorm_w', []).append(float(np.linalg.norm(gw)))
        if ftol is not None:
            cond_c = (np.linalg.norm(g) < ftol)
        else:
            cond_c = False
        cond_w = False
        if train_w and (wftol is not None):
            cond_w = (np.linalg.norm(gw) < wftol)
        if cond_c and (not train_w or cond_w):
            if verbose:
                tag = "+NU" if train_w else ""
                print(f"[Rational{tag}] Early stop at iter {i}: ||gc||={np.linalg.norm(g):.3e}{' ||gw||='+format(np.linalg.norm(gw),'.3e') if train_w else ''}")
            break
        # updates
        c -= learning_rate * g
        if train_w:
            w -= w_learning_rate * gw
            if w_clip_min is not None: w = np.maximum(w, w_clip_min)
            if w_clip_max is not None: w = np.minimum(w, w_clip_max)
        if verbose and (i % 10 == 0):
            if train_w:
                print(f"[Rational+NU] Iteration {i}, Cost: {np.sum(e*e):.6f}")
            else:
                print(f"[Rational] Iteration {i}, Cost: {np.sum(e*e):.6f}")
    return (c, w) if train_w else c

def solve_with_damped_dynamics_rational(A, b, w=None, learning_rate=0.1, damping=0.1, n_iterations=200, verbose=True, eps=1e-12, train_w=False, w_clip_min=1e-8, w_clip_max=None, log=None, ftol=None, wftol=None):
    """Damped-dynamics (with inertia) for rational model F = (A@(w*c))/(A@w).

    - c update uses velocity v_c. Force f_c = -grad_c.
    - optional w update with velocity v_w. Force f_w = -grad_w. If f·v < 0 then reset the respective velocity to zero.
    - w is clipped to [w_clip_min, w_clip_max] if provided.
    """
    if w is None: w = np.ones(A.shape[1])
    c = np.zeros(A.shape[1])
    v_c = np.zeros_like(c)
    v_w = np.zeros_like(w) if train_w else None
    for i in range(n_iterations):
        # coefficients gradient and state metrics
        g_c = grad_rational_coeffs(A, c, b, w=w, clamp=True, eps=eps)
        # optional weights gradient computed after c update decision? compute now for logging & early-stop coherence
        g_w = grad_rational_weights(A, c, b, w=w, clamp=True, eps=eps) if train_w else None
        F, _, _ = predict_rational(A, c, w=w, clamp=False, eps=eps)
        e = np.clip(F, -1.0, 1.0) - b
        # logging
        if log is not None:
            log.setdefault('cost', []).append(float(np.sum(e*e)))
            log.setdefault('gnorm_c', []).append(float(np.linalg.norm(g_c)))
            if train_w:
                log.setdefault('gnorm_w', []).append(float(np.linalg.norm(g_w)))
        # early stop
        stop_c = (ftol is not None) and (np.linalg.norm(g_c) < ftol)
        stop_w = (train_w and (wftol is not None) and (np.linalg.norm(g_w) < wftol))
        if stop_c and (not train_w or stop_w):
            if verbose:
                print(f"[Dyn{' +NU' if train_w else ''}] Early stop at iter {i}: ||gc||={np.linalg.norm(g_c):.3e}{' ||gw||='+format(np.linalg.norm(g_w),'.3e') if train_w else ''}")
            break
        # coefficients step
        f_c = -g_c
        if np.dot(f_c, v_c) < 0.0:
            v_c[:] = 0.0
        v_c = mu * v_c + learning_rate * f_c
        c += v_c

        # weights step (optional)
        if train_w:
            f_w = -g_w
            if np.dot(f_w, v_w) < 0.0:
                v_w[:] = 0.0
            v_w = w_mu * v_w + w_learning_rate * f_w
            w += v_w
            if w_clip_min is not None: w = np.maximum(w, w_clip_min)
            if w_clip_max is not None: w = np.minimum(w, w_clip_max)

        if verbose and (i % 10 == 0):
            print(f"[Dyn{' +NU' if train_w else ''}] Iteration {i}, Cost: {np.sum(e*e):.6f}")
    return (c, w) if train_w else c

def generate_basis_points(nx=4, ny=4):
    """Generate positions of bilinear basis points on an integer grid nx×ny; returns (n,2)."""
    gx, gy = np.meshgrid(np.arange(nx), np.arange(ny))
    return np.vstack([gx.ravel(), gy.ravel()]).T

def generate_grid_points(xmin=-0.5, xmax=3.5, ymin=-0.5, ymax=3.5, res=30):
    """Generate a regular grid and stacked point list; returns (X,Y,points)."""
    x, y = np.meshgrid(np.linspace(xmin, xmax, res), np.linspace(ymin, ymax, res))
    pts = np.vstack([x.ravel(), y.ravel()]).T
    return x, y, pts

def _normalize_circles(circles):
    """Normalize circles spec to (centers (m,2), radii (m,)). Accepts list[(cx,cy,r)] or array (m,3)."""
    if circles is None: return None, None
    arr = np.asarray(circles, dtype=float)
    if arr.ndim == 1 and arr.size == 3: arr = arr.reshape(1, 3)
    assert arr.shape[1] == 3, "circles must be (m,3) of (cx,cy,r)"
    centers, radii = arr[:, :2], arr[:, 2]
    return centers, radii

def points_in_any_circle(points, circles):
    """Return boolean mask of points inside any of the provided circles."""
    centers, radii = _normalize_circles(circles)
    if centers is None or len(centers) == 0: return np.zeros(points.shape[0], dtype=bool)
    d = points[:, None, :] - centers[None, :, :]
    d2 = (d[:, :, 0]**2 + d[:, :, 1]**2)
    inside = d2 < (radii[None, :]**2)
    return np.any(inside, axis=1)

def points_inside_convex_polygon(points, vertices):
    """Return boolean mask of points inside a convex polygon defined CCW by vertices (k,2)."""
    if vertices is None: return np.zeros(points.shape[0], dtype=bool)
    V = np.asarray(vertices, dtype=float)
    k = V.shape[0]
    Vi = V
    Vj = np.roll(V, -1, axis=0)
    e = Vj - Vi
    # cross(e, p-vi) >= 0 for left side (assuming CCW)
    pvi_x = points[:, 0][:, None] - Vi[None, :, 0]
    pvi_y = points[:, 1][:, None] - Vi[None, :, 1]
    cross = e[None, :, 0] * pvi_y - e[None, :, 1] * pvi_x
    left = cross >= 0.0
    return np.all(left, axis=1)

def reference_from_shapes(points, circles=None, polygon=None):
    """Combine shapes: inside if in any circle OR outside polygon; map to +/-1 at the end."""
    n = points.shape[0]
    in_c = points_in_any_circle(points, circles) if circles is not None else np.zeros(n, dtype=bool)
    in_p = points_inside_convex_polygon(points, polygon) if polygon is not None else None
    if in_p is None:
        inside = in_c
    elif circles is None:
        inside = ~in_p
    else:
        inside = in_c | (~in_p)
    return np.where(inside, 1.0, -1.0)

def evaluate_reference(points, ref_fn=reference_shape, circles=None, polygon=None, **kwargs):
    """Evaluate reference implicit shape at points. If circles or polygon given, use boolean masks then map to +/-1."""
    if (circles is not None) or (polygon is not None):
        return reference_from_shapes(points, circles=circles, polygon=polygon)
    return ref_fn(points[:, 0], points[:, 1], **kwargs)

def build_A(points, basis_points, kind='bilinear'):
    """Build matrix A with chosen basis at points: kind in {'bilinear','bilinear_aug','bilinear_power','bilinear_nu','bspline','rbf2x2','rbf2x2_sq','wendland'}"""
    if kind == 'bilinear':
        return evaluate_bilinear_bases(points, basis_points)
    if kind in ('bilinear_aug',):
        # currently same as bilinear; placeholder for future combined augmentations
        return evaluate_bilinear_bases(points, basis_points)
    if kind in ('bilinear_power','bilinear_p','bilinear_pow'):
        # raise bilinear basis to power P_POWER
        A0 = evaluate_bilinear_bases(points, basis_points)
        return np.power(A0, P_POWER)
    if kind in ('bilinear_nu','bilinear_nonuniform','bilinear_w'):
        # same bilinear A; non-uniform weights training is handled in solver path
        return evaluate_bilinear_bases(points, basis_points)
    elif kind in ('bspline', 'cubic', 'cubic_bspline', 'b_spline'):
        return evaluate_bspline_bases(points, basis_points)
    if kind in ('radial', 'rbf', 'rbf2x2', 'nn_rbf'):
        return evaluate_rbf2x2_bases(points, basis_points, kernel='r2_linear')
    if kind in ('rbf2x2_sq', 'rbf_sq', 'radial_sq'):
        return evaluate_rbf2x2_bases(points, basis_points, kernel='r2_quadratic')
    if kind in ('wendland', 'wendland_c2', 'rbf_wendland'):
        return evaluate_rbf2x2_bases(points, basis_points, kernel='wendland_c2')
    raise ValueError(f"Unknown basis kind: {kind}")

def fit_coeffs(A, b, learning_rate=0.1, n_iterations=200, log=None, ftol=None):
    """Fit coefficients with clamped error using gradient descent.

    Extended with optional logging and early-stop threshold on gradient norm.
    """
    return solve_with_gradient_descent(A, b, learning_rate=learning_rate, n_iterations=n_iterations, log=log, ftol=ftol)

# ================================
# Generic, model-agnostic relaxators
# ================================

def gradient_descent_relax(force_fn, u0, dt=0.1, n_iterations=200, ftol=None, wftol=None, log=None, verbose=True, project=None, mass=None):
    """Model-agnostic gradient descent on DOF vector u.

    force_fn(u) -> (f, E) where f is descent force (-grad-like) and E is scalar error magnitude (e.g., ||residual||).
    project(u)   optional projection/clipping applied after each update.
    Logging: if log is a list, append tuple (E, ||f||).
    """
    u = np.array(u0, dtype=float, copy=True)
    mass_vec = 1.0 if (mass is None) else np.array(mass, dtype=float, copy=False)
    for i in range(n_iterations):
        f, E = force_fn(u)
        Fn = float(np.linalg.norm(f))
        if isinstance(log, list):
            log.append((float(E), Fn))
        # early stop on force magnitude only (simple and general)
        if (ftol is not None) and (Fn < ftol):
            if verbose:
                print(f"[GD-relax] Early stop at iter {i}: |F|={Fn:.3e}, |E|={float(E):.6f}")
            break
        # update with mass scaling
        u = u + dt * (f / mass_vec)
        if project is not None:
            u = project(u)
        if verbose and (i % 10 == 0):
            print(f"[GD-relax] Iter {i}, |E|={float(E):.6f}, |F|={Fn:.3e}")
    return u

def dynamical_relaxation(force_fn, u0, dt=0.1, damping=0.1, n_iterations=200, ftol=None, wftol=None, log=None, verbose=True, project=None, mass=None):
    """Model-agnostic damped dynamics on DOF vector u with velocity.

    Velocity update: v = v * (1 - damping) + dt * f; reset v to 0 if f opposes v.
    Logging: if log is a list, append tuple (E, ||f||).
    """
    u = np.array(u0, dtype=float, copy=True)
    v = np.zeros_like(u)
    mass_vec = 1.0 if (mass is None) else np.array(mass, dtype=float, copy=False)
    cdamp = 1.0 - float(damping)
    for i in range(n_iterations):
        f, E = force_fn(u)
        Fn = float(np.linalg.norm(f))
        if isinstance(log, list):
            log.append((float(E), Fn))
        # early stop on force magnitude only
        if (ftol is not None) and (Fn < ftol):
            if verbose:
                print(f"[Dyn-relax] Early stop at iter {i}: |F|={Fn:.3e}, |E|={float(E):.6f}")
            break
        # dynamics update with mass scaling
        f_eff = f / mass_vec
        if np.dot(f_eff, v) < 0.0:
            v[:] = 0.0
        v = cdamp * v + dt * f_eff
        u = u + v
        if project is not None:
            u = project(u)
        if verbose and (i % 10 == 0):
            print(f"[Dyn-relax] Iter {i}, |E|={float(E):.6f}, |F|={Fn:.3e}")
    return u

# ================================
# Model-specific force/project factories
# ================================

def make_force_linear(A, b):
    """Return (force_fn, project, u_dim) for linear clamped-error model with DOFs u=c.

    force(u) returns (f, E) where E = ||e|| with e = clip(Au)-b.
    """
    n = A.shape[1]
    def force(u):
        F = A @ u
        F_c = np.clip(F, -1.0, 1.0)
        e = F_c - b
        dphi = ((F > -1.0) & (F < 1.0)).astype(float)
        g = A.T @ (dphi * e)
        f = -g
        E = float(np.linalg.norm(e))
        return f, E
    def project(u):
        return u
    return force, project, n

def make_force_rational(A, b, eps=1e-12, clamp=True, wclip_min=1e-8, wclip_max=None, train_w=True):
    """Return (force_fn, project, pack/unpack) for rational model.

    If train_w is True, DOFs u = [c(0..n-1), w(0..n-1)], else u=c and w is implicit ones.
    force(u) returns (f, E) where E = ||e|| with e = clip(F)-b and F = (A@(w*c))/(A@w).
    """
    n = A.shape[1]
    def unpack(u):
        if train_w:
            return u[:n], u[n:]
        else:
            return u, np.ones(n)
    def pack(c, w):
        return np.concatenate([c, w]) if train_w else np.array(c, copy=True)
    def force(u):
        c, w = unpack(u)
        # grads (with clamped error)
        g_c = grad_rational_coeffs(A, c, b, w=w, clamp=clamp, eps=eps)
        f_c = -g_c
        g_w = grad_rational_weights(A, c, b, w=w, clamp=clamp, eps=eps) if train_w else None
        f_w = -g_w if train_w else None
        # error magnitude
        F, _, _ = predict_rational(A, c, w=w, clamp=False, eps=eps)
        e = np.clip(F, -1.0, 1.0) - b
        E = float(np.linalg.norm(e))
        f = np.concatenate([f_c, f_w]) if train_w else f_c
        return f, E
    def project(u):
        if not train_w: return u
        c, w = unpack(u)
        if wclip_min is not None: w = np.maximum(w, wclip_min)
        if wclip_max is not None: w = np.minimum(w, wclip_max)
        return pack(c, w)
    return force, project, (n if not train_w else 2*n), unpack

def imshow_ax(ax, Z, extent, title, vlims=None, cmap='bwr', points=None):
    """Helper for imshow on an axis; returns the image handle."""
    if vlims is None: vlims = (-np.abs(Z).max(), np.abs(Z).max())
    im = ax.imshow(Z, extent=extent, origin='lower', cmap=cmap, vmin=vlims[0], vmax=vlims[1])
    ax.set_title(title)
    if points is not None:
        ax.scatter(points[:, 0], points[:, 1], c='yellow', s=1, marker='o', label='Basis Points')
    return im

def plot_results(ref_Z, model_Z, extent, basis_points, sample_points=None, sample_values=None, vmin=-1.0, vmax=1.0, cmap='bwr'):
    """Plot reference, model, and error using shared imshow helper."""
    fig, axes = plt.subplots(2,2, figsize=(12,12))

    err_Z = model_Z - ref_Z
    RMSE  = np.sqrt(np.mean(err_Z**2))
    print(f"\nFinal Root Mean Square Error (on plot grid): {RMSE:.4f}")

    im0 = imshow_ax(axes[0,0], ref_Z,   extent, "Reference Function", vlims=(vmin, vmax), cmap=cmap, points=basis_points)
    im3 = imshow_ax(axes[0,1], err_Z,   extent, "Error (Model - Reference)",  cmap=cmap, points=basis_points)
    im1 = imshow_ax(axes[1,0], model_Z, extent, "Model (clamped)",  vlims=(vmin, vmax), cmap=cmap, points=basis_points)
    im2 = imshow_ax(axes[1,1], model_Z, extent, "Model (full)",     cmap=cmap, points=basis_points)
    #axes[2].scatter(basis_points[:, 0], basis_points[:, 1], c='yellow', s=100, edgecolor='black', marker='o')

    fig.colorbar(im2, ax=axes.ravel().tolist(), orientation='vertical', fraction=0.02, pad=0.04)
    plt.suptitle("Shape Fitting using Bilinear Interpolation")
    #plt.show()

def plot_opt_log(log, bLog=True):
    """Plot optimization convergence.

    Preferred format: log is a list of tuples (|E|, |F|) per iteration.
    Fallback: if log is a dict with keys like 'cost' and 'gnorm'/'gnorm_c', convert minimally.
    """
    if log is None or len(log) == 0: return
    # preferred path: list of (E, F)
    if isinstance(log, list):
        E = np.array([t[0] for t in log], dtype=float)
        F = np.array([t[1] for t in log], dtype=float)
    else:
        # minimal fallback for old dict logs
        cost = np.asarray(log.get('cost', []), dtype=float)
        E = np.sqrt(cost) if cost.size > 0 else np.array([])
        F = np.asarray(log.get('gnorm', log.get('gnorm_c', [])), dtype=float)
    it = np.arange(max(len(E), len(F)))
    fig, axes = plt.subplots(1, 2, figsize=(10,4))
    if len(E) > 0:
        axes[0].plot(np.arange(len(E)), E, '-k')
    axes[0].set_title('|Error|'); axes[0].set_xlabel('iter'); axes[0].set_ylabel('||e||'); axes[0].grid(True, alpha=0.3)
    axes[0].set_yscale('log' if bLog else 'linear')
    if len(F) > 0:
        axes[1].plot(np.arange(len(F)), F, '-b')
    axes[1].set_title('|Force|'); axes[1].set_xlabel('iter'); axes[1].set_ylabel('||f||'); axes[1].grid(True, alpha=0.3)
    axes[1].set_yscale('log' if bLog else 'linear')
    fig.tight_layout()

def plot_single_basis_2d(basis_points, idx, extent, kind='bilinear', res=150):
    """Plot a single basis function as 2D map by taking column idx of A on a grid."""
    xmin, xmax, ymin, ymax = extent
    _, _, pts = generate_grid_points(xmin, xmax, ymin, ymax, res=res)
    A = build_A(pts, basis_points, kind=kind)
    Z = A[:, idx].reshape(res, res)
    plt.figure(figsize=(5,4))
    plt.imshow(Z, extent=extent, origin='lower', cmap='viridis')
    plt.colorbar(label='Basis value')
    plt.scatter(basis_points[:,0], basis_points[:,1], c='w', s=10)
    plt.title(f"Single basis 2D map (kind={kind}, idx={idx})")
    plt.tight_layout(); 
    #plt.show()

def plot_single_basis_1d(basis_points, idx, kind='bilinear', p0=None, p1=None, n=500):
    """Plot a single basis along a 1D line segment p(s)=p0+(p1-p0)*s, s in [0,1]."""
    if p0 is None or p1 is None:
        # default: horizontal mid-line across the grid span
        xs = np.unique(basis_points[:,0]); ys = np.unique(basis_points[:,1])
        x0, x1 = xs.min()-0.5, xs.max()+0.5
        y = 0.5*(ys.min()+ys.max())
        p0 = np.array([x0, y]); p1 = np.array([x1, y])
    t = np.linspace(0.0, 1.0, n)
    pts = (1.0-t)[:,None]*p0[None,:] + t[:,None]*p1[None,:]
    A = build_A(pts, basis_points, kind=kind)
    yv = A[:, idx]
    plt.figure(figsize=(6,3))
    plt.plot(t, yv, '-k')
    plt.xlabel('t along line'); plt.ylabel('basis value')
    plt.title(f"Single basis 1D slice (kind={kind}, idx={idx})")
    plt.grid(True, alpha=0.3)
    plt.tight_layout(); 
    #plt.show()

if __name__ == "__main__":

    # run like this:
    # python -u -m subPixelContour.subPixelContour

    # Parameters via argparse
    parser = argparse.ArgumentParser(description="Sub-pixel contour fitting demo")
    parser.add_argument('--size', type=int, default=None, help='number of node points per axis (sets nx=ny=size)')
    parser.add_argument('--nx', type=int, default=None, help='number of node points in x (overrides --size)')
    parser.add_argument('--ny', type=int, default=None, help='number of node points in y (overrides --size)')
    parser.add_argument('--res', type=int, default=100, help='sampling grid resolution per axis')
    parser.add_argument('--basis', type=str, default='bilinear', choices=['bilinear','bilinear_aug','bilinear_power','bilinear_nu','bspline','rbf2x2','rbf2x2_sq','wendland'], help='basis kind to use')
    parser.add_argument('--rational', action='store_true', help='use rational (normalized) model: F=(A@(w*c))/(A@w)')
    parser.add_argument('--p', type=float, default=2.0, help='power exponent for bilinear_power (p>0)')
    parser.add_argument('--nonuniform', action='store_true', help='optimize non-uniform weights w in rational model (learn w_i >= 0)')
    parser.add_argument('--test_const', action='store_true', help='diagnostic: verify rational constancy between two equal neighbors')
    parser.add_argument('--niter', type=int, default=200, help='number of iterations for gradient descent')
    # optimizer controls
    parser.add_argument('--opt', type=str, default='gd', choices=['gd','dyn'], help='optimizer: gd=gradient descent (default), dyn=damped dynamics with inertia')
    # generic relaxator params
    parser.add_argument('--dt', type=float, default=0.1, help='time step for relaxation (aka learning rate). If not set, falls back to --lr')
    parser.add_argument('--damping', type=float, default=0.1, help='damping in [0,1]; velocity multiplier is (1-damping). If not set, derived from legacy --mu as damping=1-mu')
    # legacy params (kept for compatibility)
    parser.add_argument('--lr', type=float, default=0.1, help='[deprecated] learning rate; use --dt instead')
    parser.add_argument('--wlr', type=float, default=None, help='learning rate for weights w (rational). Default: 0.5*lr')
    parser.add_argument('--mu', type=float, default=0.9, help='momentum (damping) coefficient for damped dynamics')
    parser.add_argument('--wmu', type=float, default=None, help='[deprecated] momentum for weights w; not used by generic relaxators')
    parser.add_argument('--wclip_min', type=float, default=1e-8, help='min clip for weights w when training non-uniform (to keep denominator positive)')
    parser.add_argument('--wclip_max', type=float, default=None, help='max clip for weights w when training non-uniform (optional)')
    # masses per DOF kind
    parser.add_argument('--mass_c', type=float, default=1.0, help='mass for coefficient DOFs (higher = slower updates)')
    parser.add_argument('--mass_w', type=float, default=1.0, help='mass for weight DOFs (higher = slower updates)')
    # early-stop thresholds and plotting
    parser.add_argument('--ftol', type=float, default=None, help='force/gradient norm threshold for coefficients (early stop)')
    parser.add_argument('--wftol', type=float, default=None, help='force/gradient norm threshold for weights w (early stop, rational NU)')
    parser.add_argument('--plot_opt_log', action='store_true', help='plot optimization log (|Error|, |Force|) after run')
    args = parser.parse_args()

    nx = args.nx if args.nx is not None else (args.size if args.size is not None else 10)
    ny = args.ny if args.ny is not None else (args.size if args.size is not None else 10)
    sample_res = args.res
    basis_kind = args.basis
    # set module-level power exponent
    P_POWER = float(args.p)
    assert P_POWER > 0.0, "--p must be > 0"

    xmin, xmax, ymin, ymax = -0.5, nx-0.5, -0.5, ny-0.5
    # xmin, xmax, ymin, ymax = -0.5, 3.5, -0.5, 3.5

    basis_points = generate_basis_points(nx, ny)
    sx, sy, sample_points = generate_grid_points(xmin, xmax, ymin, ymax, res=sample_res)

    # Example shapes
    circles = [
        (nx*0.35, ny*0.35, 1.8),
        (nx*0.65, ny*0.65, 1.2),
    ]
    polygon = np.array([
        [0.0,    0.0],
        [nx-2.0, 1.0],
        [nx-2.0, ny-2.0],
        [1.0,    ny-2.0],
    ], dtype=float)  # CCW

    b = evaluate_reference(sample_points, circles=circles, polygon=polygon)
    # For bilinear_power we fit on transformed targets: b_fit = sign(b)*|b|^p
    if basis_kind in ('bilinear_power','bilinear_p','bilinear_pow'):
        b_fit = np.sign(b) * (np.abs(b) ** P_POWER)
    else:
        b_fit = b
    # basis_kind set by argparse
    A = build_A(sample_points, basis_points, kind=basis_kind)
    # decide if we train non-uniform weights
    train_w = args.nonuniform or (basis_kind in ('bilinear_nu','bilinear_nonuniform','bilinear_w'))
    learned_w = None
    # If user chose 'bilinear_nu', enforce rational model
    use_rational = args.rational or (basis_kind in ('bilinear_nu','bilinear_nonuniform','bilinear_w'))
    # map generic relaxator params
    dt = float(args.dt) if (args.dt is not None) else float(args.lr)
    # cdamp legacy: mu ~ (1 - damping) => damping = 1 - mu
    mu_legacy = float(args.mu)
    damping = float(args.damping) if (args.damping is not None) else max(0.0, min(1.0, 1.0 - mu_legacy))
    # optimization log: prefer simple list of tuples (|E|, |F|) per iteration
    opt_log = []

    # choose model force and projectors
    if use_rational:
        force_fn, project, u_dim, unpack = make_force_rational(A, b_fit, clamp=True, wclip_min=args.wclip_min, wclip_max=args.wclip_max, train_w=train_w)
        if train_w:
            n = A.shape[1]
            u0 = np.concatenate([np.zeros(n), np.ones(n)])
        else:
            u0 = np.zeros(A.shape[1])
    else:
        force_fn, project, u_dim = make_force_linear(A, b_fit)
        u0 = np.zeros(u_dim)

    # build mass vector
    if use_rational and train_w:
        n = A.shape[1]
        mass = np.concatenate([np.full(n, args.mass_c), np.full(n, args.mass_w)])
    else:
        mass = np.full(A.shape[1], args.mass_c)

    # choose relaxator
    if args.opt == 'gd':
        u = gradient_descent_relax(force_fn, u0, dt=dt, n_iterations=args.niter, ftol=args.ftol, wftol=args.wftol, log=opt_log, verbose=True, project=project, mass=mass)
    else:
        u = dynamical_relaxation(force_fn, u0, dt=dt, damping=damping, n_iterations=args.niter, ftol=args.ftol, wftol=args.wftol, log=opt_log, verbose=True, project=project, mass=mass)

    # unpack solution
    learned_w = None
    if use_rational and train_w:
        c, w = unpack(u)
        coeffs, learned_w = c, w
    else:
        coeffs = u

    # 5) Report coefficients in grid form
    print("Optimal Expansion Coefficients:")
    print(coeffs.reshape(ny, nx))
    if learned_w is not None:
        print("\nLearned Non-Uniform Weights (w):")
        print(learned_w.reshape(ny, nx))

    extent = [xmin, xmax, ymin, ymax]


    #plot_res = 100
    #sx, sy, plot_points = generate_grid_points(xmin, xmax, ymin, ymax, res=plot_res)
    plot_points = sample_points; plot_res = sample_res

    ref_values_plot   = evaluate_reference(plot_points, circles=circles, polygon=polygon).reshape(plot_res, plot_res)
    A_plot            = build_A(plot_points, basis_points, kind=basis_kind)
    # helper: signed power
    def signed_pow(x, p):
        return np.sign(x) * (np.abs(x) ** p)
    # predict linear or rational on transformed basis
    if use_rational:
        F_lin = predict_rational(A_plot, coeffs, w=(learned_w if learned_w is not None else None), clamp=False)[0].reshape(plot_res, plot_res)
    else:
        F_lin = (A_plot @ coeffs).reshape(plot_res, plot_res)
    # invert the target transform for visualization if bilinear_power was used
    if basis_kind in ('bilinear_power','bilinear_p','bilinear_pow'):
        model_values_plot = signed_pow(F_lin, 1.0 / P_POWER)
    else:
        model_values_plot = F_lin

    plot_results(ref_values_plot, model_values_plot, extent, basis_points, sample_points=sample_points, sample_values=b)
    if args.plot_opt_log:
        plot_opt_log(opt_log)
        plt.show()

    # DEBUG: visualize a single basis function in 2D and 1D
    # pick center node by (ix,iy)
    ix = nx//2; iy = ny//2
    idx = iy*nx + ix
    print(f"\n# DEBUG single basis: kind={basis_kind}, node=({ix},{iy}), idx={idx}")
    plot_single_basis_2d(basis_points, idx, extent, kind=basis_kind, res=200)
    # 1D slice across the center horizontally
    xs = np.unique(basis_points[:,0]); 
    ys = np.unique(basis_points[:,1])
    p0 = np.array([xs.min()-0.5, iy])
    p1 = np.array([xs.max()+0.5, iy])
    plot_single_basis_1d(basis_points, idx, kind=basis_kind, p0=p0, p1=p1, n=800)

    plt.show()

    # Optional diagnostic: constancy along a line for two equal neighbors with rational model
    if args.test_const:
        if basis_kind not in ('rbf2x2','rbf2x2_sq','wendland'):
            print("[test_const] Skipped: enable with --basis rbf2x2 (or rbf2x2_sq, wendland)")
        else:
            # pick horizontal neighbors at row iy
            i0 = max(0, nx//2 - 1); i1 = i0 + 1; j = ny//2
            c = np.zeros(nx*ny)
            c[j*nx + i0] = 0.7
            c[j*nx + i1] = 0.7
            # sample points along y=j exactly
            xs_line = np.linspace(i0, i1, 200)
            pts_line = np.stack([xs_line, np.full_like(xs_line, j, dtype=float)], axis=1)
            A_line = build_A(pts_line, basis_points, kind=basis_kind)
            F_line, _, _ = predict_rational(A_line, c, w=None, clamp=False)
            print(f"[test_const] two equal neighbors (c=0.7) on x∈[{i0},{i1}] y={j}: min={F_line.min():.6f}, max={F_line.max():.6f}, std={F_line.std():.6e}")
