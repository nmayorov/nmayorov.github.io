---
title: 'Ellipsoid fitting'
date: 2025-11-05
katex: true
plotly: true
---

In this note a problem of ellipsoid fitting to a set of points in arbitrary number of dimensions is considered.
The algorithm development was motivated by the task of magnetometer calibration and some example results for that are given at the end.

# Problem formulation

A hyper-ellipsoid in $m$ dimensions is characterized by a center $c \in R^m$ and a symmetric positive-definite shape matrix $P \in R^{m \times m}$ which encodes length and direction of its axes.
Coordinates of a point $x \in R^m$ lying on an ellipsoid satisfy the equation:
$$
(x - c)^T P^{-1} (x - c) = 1
$$

Given a set of points $\\{x_i, i = 1, 2 \ldots, n\\}$ supposedly belonging to an ellipsoid, but recorded with noise, the problem is to determine the parameters of such ellipsoid in some optimal sense.
One natural approach to mathematically formulate this problem is a nonlinear least-squares minimization problem:
$$
\text{minimize } F(c, P) = \sum_{i = 1}^n f_i^2(c, P) \text{ with respect to } c \text{ and positive-definite } P, \\\\
\text{where } f_i(c, P) = (x_i - c)^T P^{-1} (x_i - c) - 1
$$

The proposed innovation $f_i$ represents a deviation of point's "level" from the ellipsoid level of 1.
It is also formally a deviation from 1 of squared Mahalanobis distance between the ellipsoid center and a point using the ellipsoid shape matrix $P$.

For a circle it comes down to the deviation of the squared distance between the center and a point from the squared radius, normalized by the radius:

- $d^2 = (x - x_0)^2 + (y - y_0)^2 - r^2$ --- deviation of the squared distance from the squared radius
- $f = \dfrac{d^2}{r^2} = \dfrac{(x - x_0)^2}{r^2} + \dfrac{(y - y_0)^2}{r^2} - 1$ --- the normalized deviation

In general such normalized innovations make formulation agnostic to the ellipsoid scale, which is a positive thing for an optimization algorithm.

# Approach to the solution

To solve the proposed problem by an existing algorithm like [least_squares](https://scipy.github.io/devdocs/reference/generated/scipy.optimize.least_squares.html#scipy.optimize.least_squares) we need to clarify some issues.

First, to ensure that $P$ is positive-definite we parameterize it by lower-triangular Cholesky factors $P^{-1} = L^T L$.
The problem becomes unconstrained in terms of the matrix $L$.
It can be stated as finding a transform $z = L(x - c)$ to make points $z$ lie approximately on a unit sphere centered at origin.

Second, as seen in practice, when the points occupy only a small sector of the ellipsoid (especially if concentrated in one "hemisphere"), the problem becomes ill-conditioned where many solutions with comparable cost exist.
Moreover, such solutions often tend to drift towards large and skewed ellipsoids, which don't make sense in practice.
A remedy to that is to introduce a regularization term which will penalize ellipsoids which strongly deviate from a sphere of a given radius $r$ (algorithm parameter):
$$
P \approx r^2 I \\\\
r^2 P^{-1} \approx I \\\\
r L \approx I
$$
The constraint comes down to adding $m (m + 1) / 2$ residuals of the form: 
$$
f_r = r L - I
$$
Where the matrix is flattened to make a vector.

We also need to properly balance the fitting and tregularization residuals.
To do that, we divide the two kind of residuals by square root of their count and also introduce a factor $\alpha$ to control the extent of the regularization:
$$
f_i^\prime = \frac{1}{\sqrt{n}} f_i \\\\
f_r^\prime = \sqrt{\frac{2 \alpha}{m (m + 1)}} f_r
$$
The total expected cost function will be:
$$
\left< F^\prime \right> = \left<F_f\right> + \alpha \left<F_r\right>
$$
Where $F_f$ is a single point fit error cost and $F_r$ is the regularization error cost per matrix element.
The point is that in this form the error terms are agnostic to the number of points and dimensions, whereas $\alpha$ controls the level of regularization in a straightforward fashion.
The reasonable values of $\alpha$ are in the range from 0 to 1, with 0 obviously being the case without regularization.

## Jacobian formulas

It's preferable to provide an analytic Jacobian of residuals with respect to unknown parameters for least-squares minimization. 
For convenience in this section the index of a point is omitted and vector component indices are used instead.

The Jacobian (gradient vector) with respect to $c$ is a known result from matrix calculus:
$$
\frac{\partial f}{\partial c} = 2 L L^T (c - x)
$$

It is easier to obtain the Jacobian with respect to the elements of $L$ considering an example with 3 dimensions.
Let denote $y = x - c$, then
$$
z = L y = \begin{bmatrix}
L_{11} y_1 \\\\
L_{12} y_1 + L_{22} y_2 \\\\
L_{13} y_1 + L_{23} y_2 + L_{33} y_3
\end{bmatrix} \\\\
f = z^T z - 1 = (L_{11} y_1)^2 + (L_{12} y_1 + L_{22} y_2)^2 + (L_{13} y_1 + L_{23} y_2 + L_{33} y_3)^3 - 1
$$
Then it's easy to compute derivatives arranging them in a matrix for convenient representation:
$$
\frac{\partial f}{\partial L} = 2 \begin{bmatrix}
z_1 y_1 & & \\\\
z_2 y_1 & z_2 y_2 & \\\\
z_3 y_1 & z_3 y_2 & z_3 y_3
\end{bmatrix} \\\\
\frac{\partial f}{\partial L_{ij}} = 2 z_i y_j
$$

The regularization residual is linear in $L$ and thus the Jacobian is constant:
$$
\frac{\partial f_r}{\partial L_{ij}} = r
$$

## Parameter initialization

The center $c$ is initialized as a mean value of all given points (centroid).
For the scale matrix $L$ there are two options:

- If the regularization radius $r$ is given, $L$ is initialized from that
- If not given, $L$ is initialized as diagonal from a bounding box with sides along coordinate axes

In both cases no initial axes rotations are assumed.

# Implementation in Python

Here is a possible implementation of this algorithm relying on [least_squares](https://scipy.github.io/devdocs/reference/generated/scipy.optimize.least_squares.html#scipy.optimize.least_squares).
For simplicity the core optimization is run with default parameters.
The code should not be considered "library quality", but it is likely to work in most cases.

```Python
import numpy as np
import scipy


def fit_ellipsoid(X, alpha=0, expected_radius=None):
    """Fit ellipsoid to points.

    The function finds parameters of hyper-ellipsoid which in least-squares sense
    describes a given set of points. The ellipsoid characterized by its center `c` and
    positive-definite shape matrix `P = (L^T L)^-1`, where `L` is lower triangular
    matrix. The residuals in least-squares minimization are formed as normalized
    squared distances of the form `f = (x - c)^T L^T L (x - c) - 1`.

    The algorithm uses optional regularization which penalizes deviation of the
    ellipsoid from a sphere of a given radius.

    Parameters
    ----------
    X : array_like, shape (n_points, n_dimensions)
        Array of points. Each row represents a point.
    alpha : float, optional
        Ellipsoid shape regularization parameter. Default is 0.
    expected_radius : float or None, optional
        Expected radius of a sphere from which the ellipsoid shape must not deviate too
        much. If None, regularization is disabled.

    Returns
    -------
    c : ndarray, shape (n_dimensions,)
        Ellipsoid center.
    P : ndarray, shape (n_dimensions, n_dimensions)
        Positive-definite ellipsoid shape matrix encoding direction and length of axes.
    L : ndarray, shape (n_dimensions, n_dimensions)
        Lower-triangular transform matrix.
    L_inv : ndarray, shape (n_dimensions, n_dimensions)
        Inverse of L.
    opt_result : scipy.optimize.OptimizationResult
        Raw optimization result from scipy.optimize.least_squares.
    """
    X = np.asarray(X)
    n, m = X.shape

    row_ind, col_ind = np.tril_indices(m)
    I_flat = np.eye(m)[row_ind, col_ind]

    def unpack_x(x):
        c = x[:m]
        l = x[m:]
        L = np.zeros((m, m))
        L[row_ind, col_ind] = l
        return c, l, L

    def fun(x):
        c, l, L = unpack_x(x)
        Y = X - c
        Z = (L @ Y.T).T
        return np.hstack(
            [
                1 / n**0.5 * (np.sum(Z**2, axis=1) - 1),
                (alpha / len(l)) ** 0.5 * (expected_radius * l - I_flat),
            ]
        )

    def jac(x):
        c, l, L = unpack_x(x)
        Y = X - c
        Z = (L @ Y.T).T
        ZY = np.einsum("...i,...j->...ij", Z, Y)
        return np.block(
            [
                [-2 / n**0.5 * Z @ L, 2 / n**0.5 * ZY[:, row_ind, col_ind]],
                [
                    np.zeros((len(l), m)),
                    (alpha / len(l)) ** 0.5 * expected_radius * np.eye(len(l)),
                ],
            ]
        )

    x0 = np.hstack([c, L[row_ind, col_ind]])
    opt_result = scipy.optimize.least_squares(fun, x0, jac=jac, verbose=2)

    c, _, L = unpack_x(opt_result.x)
    L_inv = scipy.linalg.solve_triangular(L, np.eye(m), lower=True)
    P = L_inv @ L_inv.T

    return c, P, L, L_inv, opt_result
```

# Examples of fitting

Considering geometric nature of the problem it seems sufficient to demonstrate and assess fitting quality visually.

## Fitting ellipse on synthetic data

Here we simply generate points belonging to an ellipse with semi-axes of 4 and 2 tilted by 45 degrees and consider the algorithm behavior depending on which points are passed to the algorithm.
All runs were executed with `expected_radius = 3`.

[![ellipse_full](figs/ellipse_full.svg)](figs/ellipse_full.svg)
[![ellipse_full](figs/ellipse_2_sectors.svg)](figs/ellipse_2_sectors.svg)
[![ellipse_full](figs/ellipse_1_sector.svg)](figs/ellipse_1_sector.svg)

We can see that when the points fully define the shape and location of the ellipse the fit is good and regularization of 0.01 and 0.1 doesn't distort the shape significantly, nor it is necessary.
In case of the points occupying only one short arc, the problem is ill-conditioned and the fitted shape strongly depends on a regularization parameter, whereas "correct" fit seems impossible to find.

In general to get proper and meaningful results we must avoid cases depicted on the third plot.

## Fitting ellipsoid to magnetometer measurements

One of well-known applications of ellipsoid fitting is magnetometer calibration where the aim is to estimate its bias ("hard iron") and transform matrix ("soft iron").
The idea is that in a calibrated sensor the Earth magnetic field must occupy a sphere when the magnetometer rotates around.
To estimate the aforementioned intrinsic parameters we fit an ellipsoid to the recorded set of measured vectors, provided that the sensor underwent sufficient rotations.
The estimated matrix $L$ will in fact serve as required "soft iron" calibration matrix.

Note, that there are some intricacies of 2D vs 3D calibration and practical algorithms might rely on some other tricks and heuristics.
Here examples are shown for two simple cases where the sensor was rotated around all axes in full ranges.
No regularization were used as it seems unnecessary.

{{< plotly json="figs/calibration_1.json" >}}
{{< plotly json="figs/calibration_2.json" >}}
