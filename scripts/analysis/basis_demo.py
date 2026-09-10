#! /usr/bin/env python
# Time-stamp: <09-06-2026 a.tabas@bcbl.eu>
# Minor adaptations by m.utrosa@bcbl.eu

"""
The idea is to express the response function pe(delta_t) as a sum of basis functions.
    - Each basis function is like a bump that has a finite width
    - By adding more/less you can have different response curves
    - All of these curves are smooth: pe(delta_t + 1) ~ pe(delta_t )

The expression is:
pe(delta_t) = b_1 * basis_1(delta_t) + b_2 * basis_2(delta_t) + ... + b_K * basis_K(delta t).

The larger K, the more basis functions, the higher the resolution in delta_t.
The script produces the bumps (dashed lines) and final pe(delta_t) curve (blue line)
for different values of the b_1, b_2, ... (i.e.: betas = [ ] vector in line 53).

The way this works:
    - decide on the number of functions K
    - compute the values for the K basis functions `basis_k(delta_t)` for each of the conditions
    - build k regressors with the values of the bases functions
      (e.g.: for condition -4, plug in the values of basis_1(-4) ... basis_k(-4), in each of the regressors 1 ... k
    - run SPM with those regressors
    - get a set of betas (b_1, b_2, ... b_k) for each vector
    - by plugging in the betas in the pe(delta_t) function, you recover the response function for that voxel

How to decide K? 
    Run cross-validation to see what's the best value.
    For example, we fit several models using sessions 1 and 2 with different Ks,
    evaluate the models on the data for session 3, and take the K that best explains the data.

In addition, we could have two versions of the model.
    a.) The one below allows for asymmetric responses across positive and negative deltas.
        With a symmetric version we could afford lower Ks, which will give us more statistical power.
    b.) We could also make a version that just takes in the absolute value of delta_t and
        assumes pe(delta_t) = pe(-delta_t). With the non-symmtric version, we'll need a 
        larger K to have good response functions.

What are the basis functions?
The script has two versions of basis functions:
    a.) basis_k is a Gaussian centred in some value delta_t_k and with some variance sigma_k
    b.) basis_k are optimised for what we are trying to do with; bslines are functions shaped
        similarly to a Gaussian that have the ideal properties to use them to define smooth
        functions with a bunch of parameters

- play with the functions
- read online about the bsline functions
- for the analysis we should use b-splines
- try to implement the asymmetric version of the model
 
"""
import numpy as np
import matplotlib.pyplot as plt

# ------------- Knobs
basis = 'bspline'          # 'gaussian' or 'bspline'
betas = np.array([0.0, 0.1, 0.05, 0.8, -0.05, 0.1, -0.1, 0.0]) # Needs to have at least 8 for degree 3
sigma = 35                  # gaussian width (ignored for bspline)
uniform_knots = False       # bspline: True = even spacing, False = quantile spacing

# Your actual sampled conditions (ms), denser near 0, both signs
mags = [0, 4, 8, 13, 19, 27, 36, 48, 63, 80, 100, 125]
conditions = np.array(sorted(set(mags + [-m for m in mags])))

# ------------- Curve Fitting
K = len(betas)
dmin, dmax = conditions.min(), conditions.max()
d = np.linspace(dmin, dmax, 600)

def gaussian_basis(d, K, sigma):
    '''
        Args:
            d
            K (int): The number of basis functions (equal to no. of betas)
            sigma
    '''
    centers = np.linspace(dmin, dmax, K)
    return np.stack([np.exp(-(d - c)**2 / (2 * sigma**2)) for c in centers], axis=1)

def bspline_basis(d, K, p=3):
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
    q = np.linspace(0, 1, K - p + 1)

    # Calculate interior knot positions
    interior = np.linspace(dmin, dmax, K - p + 1) if uniform_knots else np.quantile(conditions, q)
    knots = np.concatenate([[dmin] * p, interior, [dmax] * p])

    # Construct the full knot vector
    # np.eye: returns a 2D array with ones on the diagonal and zeros elsewhere
    return np.stack([BSpline(knots, np.eye(K)[k], p)(d) for k in range(K)], axis=1)

# ------------- Design matrix of basis functions
phi = gaussian_basis(d, K, sigma) if basis == 'gaussian' else bspline_basis(d, K)

# ------------- Matrix multiplication
f = phi @ betas

# ------------- Plotting
plt.figure(figsize=(8, 4))
plt.plot(d, phi * betas, '--', lw=1, alpha=0.6)
plt.plot(d, f, 'b', lw=3, label='f(d) = sum beta_k phi_k(d)')
plt.plot(conditions, np.zeros_like(conditions), '|', color='k', ms=12, label='conditions')
plt.axhline(0, color='gray', lw=1)
plt.xlabel('d'); plt.title(f'{basis} basis, K={K}'); plt.legend()
plt.tight_layout(); plt.show()