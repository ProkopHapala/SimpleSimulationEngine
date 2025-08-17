import numpy as np
import matplotlib.pyplot as plt

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

def solve_with_gradient_descent(A, b, learning_rate=0.1, n_iterations=200, verbose=True):
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
        # Optional: Print cost for debugging
        if verbose and (i % 10 == 0):
            cost = np.sum(error**2)
            print(f"Iteration {i}, Cost: {cost}")
            
    return coeffs

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
    """Build matrix A with chosen basis at points: kind in {'bilinear','bspline'}"""
    if kind == 'bilinear':
        return evaluate_bilinear_bases(points, basis_points)
    if kind in ('bspline', 'cubic', 'cubic_bspline', 'b_spline'):
        return evaluate_bspline_bases(points, basis_points)
    raise ValueError(f"Unknown basis kind: {kind}")

def fit_coeffs(A, b, learning_rate=0.1, n_iterations=200):
    """Fit coefficients with clamped error using gradient descent."""
    return solve_with_gradient_descent(A, b, learning_rate=learning_rate, n_iterations=n_iterations)

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
    plt.show()

if __name__ == "__main__":

    # run like this:
    # python -u -m subPixelContour.subPixelContour

    # Parameters
    nx = ny    = 10
    sample_res = 100
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
    basis_kind = 'bspline'  # 'bilinear' or 'bspline'
    A = build_A(sample_points, basis_points, kind=basis_kind)
    coeffs = fit_coeffs(A, b)

    # 5) Report coefficients in grid form
    print("Optimal Expansion Coefficients:")
    print(coeffs.reshape(ny, nx))

    extent = [xmin, xmax, ymin, ymax]


    #plot_res = 100
    #sx, sy, plot_points = generate_grid_points(xmin, xmax, ymin, ymax, res=plot_res)
    plot_points = sample_points; plot_res = sample_res

    ref_values_plot   = evaluate_reference(plot_points, circles=circles, polygon=polygon).reshape(plot_res, plot_res)
    A_plot            = build_A(plot_points, basis_points, kind=basis_kind)
    model_values_plot = (A_plot @ coeffs).reshape(plot_res, plot_res)

    plot_results(ref_values_plot, model_values_plot, extent, basis_points, sample_points=sample_points, sample_values=b)
