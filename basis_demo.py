import numpy as np
import matplotlib.pyplot as plt

# --- knobs ---
basis = 'bpline'           # 'gaussian' or 'bspline'
betas = np.array([0.03, 0.1, 0.13, 0.4, 0.9, 0.03]) # Needs to have at least 8 for degree 3
sigma = 35                  # gaussian width (ignored for bspline)
uniform_knots = False       # bspline: True = even spacing, False = quantile spacing

# Your actual sampled conditions (ms), denser near 0, both signs
mags = [4, 8, 13, 19, 27, 36, 48, 63, 80, 100, 125]
conditions = np.array(sorted(set(mags + [-m for m in mags])))

# ------------- Curve fitting
K = len(betas) # The number of functions

def gaussian_basis(conditions, K, sigma):
    '''
    Parameters:
    dmin - 
    dmax -
    K -
    '''
    dmin, dmax = conditions.min(), conditions.max()
    centers = np.linspace(dmin, dmax, K)
    d = np.linspace(dmin, dmax, 600)
    return np.stack([np.exp(-(d - c)**2 / (2 * sigma**2)) for c in centers], axis=1)

def bspline_basis(K, uniform_knots, conditions, degree):
    """
    Generate B-spline basis functions for given conditions.

    Args:
        K (int): The number of basis functions to generate.
        uniform_knots (bool): If True, knots are uniformly spaced. 
                              If False, knots are placed at quantiles of `conditions`.
        conditions (array): The sampled conditions.
        degree (int): The polynomial degree of the B-spline.

    Returns:
        numpy.ndarray: A 2D array of shape (len(d), K) where each column 
                       represents the evaluation of one basis function 
                       over the domain `d`.
    """
    from scipy.interpolate import BSpline

    # Create a grid for evaluating the basis functions
    dmin, dmax = conditions.min(), conditions.max()
    q = np.linspace(0, 1, K - degree + 1)
    d = np.linspace(dmin, dmax, 600)

    # Calculate interior knot positions
    interior = np.linspace(dmin, dmax, K - degree + 1) if uniform_knots else np.quantile(conditions, q)
    knots = np.concatenate([[dmin] * degree, interior, [dmax] * degree])
    
    # Construct the full knot vector
    # np.eye: returns a 2D array with ones on the diagonal and zeros elsewhere
    return np.stack([BSpline(knots, np.eye(K)[k], degree)(d) for k in range(K)], axis=1)

# Design matrix of basis functions
phi = gaussian_basis(conditions, K, sigma) if basis == 'gaussian' else bspline_basis(K, uniform_knots, conditions)

# Matrix multiplication
f = phi @ betas

# Plotting
plt.figure(figsize=(8, 4))
plt.plot(d, phi * betas, '--', lw=1, alpha=0.6)
plt.plot(d, f, 'b', lw=3, label='f(d) = sum beta_k phi_k(d)')
plt.plot(conditions, np.zeros_like(conditions), '|', color='k', ms=12, label='conditions')
plt.axhline(0, color='gray', lw=1)
plt.xlabel('d'); plt.title(f'{basis} basis, K={K}'); plt.legend()
plt.tight_layout(); plt.show()