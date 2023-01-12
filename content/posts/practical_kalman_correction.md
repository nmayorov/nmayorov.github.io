---
title: "Implementation of Kalman correction formulas"
date: 2023-01-12
katex: true
---

This is a short note on my view how to best implement Kalman correction formulas in a programing language.

# Basic formulas

A vector measurement is linearly related to a state vector with an additive noise:
$$
z = H x + v
$$
Where $z$ and $H$ are a known vector and matrix and the noise vector $v$ has a known covariance matrix $R$.
Thus the measurement model is defined by 3 elements -- $z, H, R$.

The measurement is processed to improve the state estimate and reduce its error variance according to formulas:
$$
x^+ = x^- + K (z - H x^-) \\\\
P^+ = (I - KH) P^- \\\\
\text{with } K = P^- H^T S^{-1} \text{ and } S = H P^- H^T + R
$$
Where $x^-, P^-$ and $x^+, P^+$ are state estimates and error covariance matrices before (aprior) and after (aposteriori) the measurement was processed respectively.
The gain vector $K$ optimally combines the apriori estimate and the measurement to achieve the smallest error variance of the $x^+$.

The formulas are valid as long as the innovation covariance matrix $S$ is positive definite, which makes them almost universally applicable in practical problems.

## Scalar measurement 

If $z$ contains only a single element we can consider it to be a scalar and write the measurement equation as
$$
z = h^T x + v
$$
Let $r$ be the variance of the scalar noise $v$.

The correction formulas are more simple and elegant in this case:
$$
x^+ = x^- + (z - h^T x^-) k \\\\
P^+ = P^- - s k k^T \\\\
\text{with } k = P^- h / s \space \text{ and } \space s = h^T P^- h + r
$$
Its free from the matrix inversion as we have the scalar innovation variance $s$.

The covariance formula has a nice symmetric rank-1 update form.
This operation is known as `SYR` BLAS subroutine.
For example in Eigen C++ library it is provided by the method `SelfAdjointView::rankUpdate`.
It is recommended to use a similar method of function in your implementation if possible.

# Independent measurements and sequential updates

Consider a measurement vector composed of two independent (uncorrelated) blocks:
$$
z = \begin{bmatrix} z_1 \\\\ z_2 \end{bmatrix} \\\\
H = \begin{bmatrix} H_1 \\\\ H_2 \end{bmatrix} \\\\
R = \begin{bmatrix} R_1 & 0 \\\\ 0 & R_2 \end{bmatrix}
$$
The matrix $R$ is block diagonal which represents the aforementioned independence.

From probabilistic and logical reasoning we can conclude that processing measurement $z$ is equivalent to processing measurements $z_1$ and $z_2$ sequentially in any order.
Surprisingly it's hard to prove that using formal matrix transformations of the correction formulas.
I haven't been able to find a simple approach to that.

A proof can be given using the equivalent <<inverse covariance>> form of the correction formulas:
$$
(P^+)^{-1} x^+ = (P^-)^{-1} x^- + H^T R^{-1} z \\\\
(P^+)^{-1} = (P^-)^{-1} + H^T R^{-1} H 
$$
They are more restrictive compared to the standard formulas as requiring positive definite $P$ and $R$ matrices.

Substituting the composed measurement matrices into it we get
$$
(P^+)^{-1} x^+ = (P^-)^{-1} x^- + H_1^T R_1^{-1} z_1 + H_2^T R_2^{-1} z_2 \\\\
(P^+)^{-1} = (P^-)^{-1} + H_1^T R_1^{-1} H_1 + H_2^T R_2^{-1} H_2
$$
We see that these formulas are additive in the measurement blocks and thus the sequential processing is correct.

If $R$ and $P$ are not positive definite we can consider an artificial positive diagonal regularization, do the sequential processing in the standard form and then set the regularization to zero (in a limit sense).
This is an informal proof that the sequential processing is correct as long as the standard formulas are applicable.

# Measurements with diagonal $R$

If the matrix $R$ is diagonal then all individual measurement components are independent.
Such measurement are processed as a sequence of scalar measurements $z_i, h_i, r_i$ with
$$
h_i = (H_{:, i})^T \\\\
r_i = R_{i, i}
$$
where $H_{:, i}$ -- $i$-th row vector of $H$.

# Measurements with non-diagonal $R$

We can apply the sequential measurement processing in this case too by introducing an augmented vector:
$$
x_a = \begin{bmatrix}
x \\\\
v
\end{bmatrix}
$$

The apriori state estimate and covariance matrix for the augmented vector are:
$$
x_a^- = \begin{bmatrix}
x^- \\\\
0
\end{bmatrix} \\\\
P_a^- = \begin{bmatrix}
P^- & 0 \\\\
0 & R
\end{bmatrix}
$$

The measurement in terms of the augmented vector:
$$
z = H x + v = H_a x_a + v_a \\\\
\text{with } H_a = \begin{bmatrix}
H & I
\end{bmatrix}
$$
The dummy noise vector $v_a$ has zero covariance matrix $R_a = 0$.

Then the augmented measurement model $z, H_a, R_a$ is processed for $x_a^-, P_a^-$  using the sequential approach to get $x_a^+, P_a^+$.
The head of $x_a^+$ and the top-left corner of $P_a^+$ contain the required aposteriori $x^+$ and $P^+$.

# Overall algorithm

The following algorithm can be used to process an arbitrary measurement model $z, H, R$:

1. If $R$ is diagonal: apply the sequential processing to the original model
2. If $R$ is not diagonal: build the augmented model and apply the sequential processing to it

The core algorithm in both cases is the scalar measurement processing described above, which is easy to implement and robust.
