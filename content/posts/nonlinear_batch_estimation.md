---
title: "Nonlinear batch estimation"
date: 2022-11-30
katex: true
---

In this note I present a development of nonlinear batch estimation algorithm.

# Model description

The state estimation problem in a nonlinear system is considered.
The formulations is analogous to the [linear case]({{<ref "/content/posts/rts_as_optimization.md#problem-formulation">}}), but with nonlinear transition and measurement equations.
We use uppercase letters to denote variables participating in the nonlinear model. 
Time transition and measurement equations for which are
$$
X_{k+1} = f_k(X_k, W_k) \\\\
Z_k = h_k(X_k) + V_k 
$$
All the other assumptions remain the same. 
The task is to estimate $X_k$ for epochs $k = 0, 1, \ldots, N$.
This will be done by solving an optimization problem.

# Optimization problem

Similarly to the linear case the cost function is defined as follows
$$
\begin{split}
E(X, W) =& \frac{1}{2} (X_0 - X_0^-)^T P_0^- (X_0 - X_0^-)  \\\\
        +& \frac{1}{2} \sum_{k = 0}^{N-1} (h_k(X_k) - Z_k)^T R_k^{-1} (h_k(X_k) - Z_k) \\\\
        +& \frac{1}{2} \sum_{k = 0}^{N-1} W_k^T Q_k^{-1} W_k
\end{split}
$$
It must be optimized taking into account the time transition equations:
$$
\min_{X, W} \space E(X, W) \text{ subject to} \space X_{k + 1} = f_k(X_k, W_k)
$$
It can be classified as a nonlinear least-squares problem with nonlinear equality constraints on the variables.

# Outline of the solution method

Consider a general optimization problem of the form:
$$
\min_x f(x) \text{ subject to } c(x) = 0
$$
One established method of solving it called sequential quadratic programming (SQP) goes like this:

1. At current estimate $x_i$ form a local quadratic approximation of the Lagrange function for the problem and a linear approximation of the constraints
2. Solve this quadratic programming task to obtain step $p_k$ to the next estimate
3. Form the next estimate $x_{i+1} = x_i + \alpha_i p_i$
4. Repeat until convergence

To form the mentioned local quadratic approximation, the gradient (gradient of $f$ might be used with the equivalent result) and the Hessian of the Lagrange function needs to be formed.
The method we apply differs in how the Hessian approximation is formed and might be called a Gauss-Newton approximation to the SQP method.

In a nonlinear least-squares problem the cost function has the form
$$
f(x) = \frac{1}{2} r(x)^T r(x)
$$
The Jacobian matrix of $r(x)$ is defined as
$$
J(x) = \begin{bmatrix} \nabla r_1(x) & \nabla r_2(x) & \ldots & \nabla r_n(x) \end{bmatrix}^T
$$
Then the gradient and Gauss-Newton approximation of the Hessian are:
$$
\nabla f(x) = J^T r(x) \\\\
\nabla^2 f(x) \approx J(x)^T J(x)
$$
Opposed to the SQP method we also neglect the second derivatives of the constraints, i. e. form a quadratic model for the cost function, not for the Lagrangian.

With these approximations the quadratic subproblem is
$$
\min_{p_i} (r_i + J_i p_i)^T (r_i + J_i p_i) \text{ subject to } c_i + A_i p_i = 0 \\\\
\text{with } r_i = r(x_i), J_i = J(x_i), c_i = c(x_i), A_i = \nabla c(x) \vert_{x_i}
$$
Logic of the algorithm iterations remain the same.

It must be said that the proposed method is not widely known. 
However my intuition and limited experience tells me that it is a legit approach similar to Gauss-Newton method of solving unconstrained nonlinear least squares.

For our specific problem the linear subproblem formulated above is solved by a Kalman smoother.

# Merit function and line search

To guarantee convergence of a nonlinear optimization method a step $p$ obtained from a linear subproblem must be applied with a properly selected multiplier as $\alpha p$.
In constrained optimization there are 2 objectives: minimize the cost function and satisfy the constraints.
One approach to account for both ot them is to consider a merit function ([1], section 11.2) with $\mu > 0$:
$$
\phi(x; \mu) = f(x) + \mu \lVert c(x) \rVert_1
$$
The parameter $\alpha$ is selected to achieve a sufficient decrease of $\phi(x; \mu)$ (the Armijo condition) for some $0 < \eta < 1$:
$$
\phi(x_i + \alpha_i p_i; \mu) \leq \phi(x_i; \mu) + \eta \alpha_i D(\phi(x_i; \mu); p_i)
$$
Where $D(g(x); l)$ is a directional derivative of $g(x)$ in direction $l$:
$$
D(g(x); l) \coloneqq \lim\_{\epsilon \rightarrow 0} \frac{g(x + \epsilon l) - g(x)}{\epsilon}
$$
The derivative of the introduced merit function in the direction $p_i$ (which satisfies linearized constraints!) can be shown to be ([1], section 18.3)
$$
D(\phi(x_i; \mu); p_i) = \nabla f_i^T p_i - \mu \lVert c_i \rVert_1
$$
In order to $p_i$ be the descent direction for the merit function the directional derivative must be negative.
This can be achieved by selecting large enough $\mu$.
Informally it means that if necessary the optimizer must take steps which don't reduce the cost function, but only constraint violation.
The specific strategy of selecting $\mu$ suggested in [1] (section 18.3) is explained next.

Using the quadratic model $q(p)$ of $f(x)$ and linear model of $c(x)$ near $x_i$ it is possible to estimate the change of the merit function for step $p_i$:
$$
\phi(x_i + p_i; \mu) - \phi(x_i) \approx q_i(p_i) + \mu \lVert c_i + A_i p_i \rVert_1 - q_i(0) - \mu \lVert c_i \rVert_1 = q_i(p_i) - q_i(0) - \mu \lVert c_i \rVert_1
$$
The aim is to have this change sufficiently negative
$$
q_i(p_i) - q_i(0) - \mu \lVert c_i \rVert_1 \leq -\rho \mu \lVert c_i \rVert_1
$$
for some $0 < \rho < 1$. From this the inequality for $\mu$ follows
$$
\mu \geq \frac{q_i(p_i) - q_i(0)}{(1 - \rho) \lVert c_i \rVert_1} = \frac{\nabla f_i^T p_i + (1/2) p_i^T H_i p_i}{(1 - \rho) \lVert c_i \rVert_1}
$$
Where $H_i$ is the Hessian of the Lagrange function or its approximation evaluated at $x_i$.

If $\mu$ satisfies the above inequality and $H_i$ is positive semidefinite we have for the directional derivative 
$$
D(\phi(x_i; \mu); p_i) \leq -\rho \mu \lVert c_i \rVert_1 < 0
$$
Meaning that $p_i$ is a descent direction for the merit function, which is a requirement for the line search procedure.
Note that the above inequality holds even when $\mu$ is computed without the Hessian part. 
However usage of the second order information results in better $\mu$ selection and larger steps.

## Summary of the algorithm for $\mu$ selection

If the current $\mu$ already satisfies the inequality it must not be changed.
This simple algorithm follows:

1. Start with $\mu_0 = 1$
2. On each iteration set $\mu_i = \max\left(\mu_{i - 1}, \dfrac{q_i(p_i) - q(0)}{(1 - \rho) \lVert c_i \rVert_1}  \right)$

## Summary of the line search algorithm

To find $\alpha_i$ which satisfies the sufficient decrease condition of the merit function the following algorithm is used:

1. Set $\alpha_i = 1$
3. Compute the directional derivatives $D_i \coloneqq D(\phi(x_i; \mu); p_i)$
2. While $\phi(x_i + \alpha_i p_i; \mu) > \phi(x_i; \mu) + \eta \alpha_i D_i$ set $\alpha_i \leftarrow \tau \alpha_i$ for some $0 < \tau < 1$

In the end the estimate is updated as 
$$
x_{i + 1} = x_i + \alpha_i p_i
$$

## Numerical values of the algorithm constants

The following reasonable values might be used:
$$
\eta = 0.5 \\\\
\rho = 0.5 \\\\
\tau = 0.5
$$

# Constructing the linearized subproblem

First introduce convenient notation. Let 

* $\hat{X}_k, \hat{W}_k$ be current estimates of $X_k$ and $W_k$
* $x_k, w_k$ be correction to the current estimates

The linearization relation is defined as:
$$
X_k = \hat{X}_k + x_k \\\\
W_k = \hat{W}_k + w_k
$$

The linearized version of the measurement residual:
$$
h(X_k) - Z_k = h(\hat{X}_k + x_k) - Z_k \approx h(\hat{X}_k) + H_k x_k - Z_k = H_k x_k - z_k \\\\
\text{with } H_k = \left.\frac{\partial h(X)}{\partial X}\right\vert\_{\hat{X}_k} \text{ and } z_k = Z_k - h(\hat{X}_k)
$$

The linearized version of the time propagation equation:
$$
X_{k + 1} = f_k(X_k, W_k) \\\\
\hat{X}_{k + 1} + x\_{k + 1} = f_k(\hat{X}_k + x_k, \hat{W}_k + w_k) \approx f_k(\hat{X}_k, \hat{W}_k) + F_k x_k + G_k w_k \\\\
x\_{k + 1} = F_k x_k + G_k w_k + u_k \\\\
\text{with } F_k = \left.\frac{\partial f(X, W)}{\partial X}\right\vert\_{\hat{X}_k, \hat{W}_k}
, G_k = \left.\frac{\partial f(X, W)}{\partial W}\right\vert\_{\hat{X}_k, \hat{W}_k}
, u_k = f_k(\hat{X}_k, \hat{W}_k) - \hat{X}\_{k + 1}
$$

## Summary of the linearized subproblem

The linearized cost function
$$
\begin{split}
J(x, w) =& \frac{1}{2} (x_0 - x_0^-)^T P_0^{-1} (x_0 - x_0^-) \\\\
        +& \frac{1}{2} \sum_{k = 0}^{N - 1} (H_k x_k - z_k)^T R_k^{-1} (H_k x_k - z_k) \\\\
        +& \frac{1}{2} \sum_{k = 0}^{N - 1} (w_k - w_k^-)^T Q_k^{-1} (w_k - w_k^-)
\end{split}
$$
is minimized subject to constraints:
$$
x_{k + 1} = F_k x_k + G_k w_k + u_k
$$

Where all matrices and vectors have the following values:
$$
x_0^- = X_0^- - \hat{X}_0 \\\\
z_k = Z_k - h(\hat{X}_k)\\\\
H_k = \left.\frac{\partial h(X)}{\partial X}\right\vert\_{\hat{X}_k} \\\\
w_k^- = -\hat{W}_k \\\\
F_k = \left.\frac{\partial f(X, W)}{\partial X}\right\vert\_{\hat{X}_k, \hat{W}_k} \\\\
G_k = \left.\frac{\partial f(X, W)}{\partial W}\right\vert\_{\hat{X}_k, \hat{W}_k} \\\\
u_k = f_k(\hat{X}_k, \hat{W}_k) - \hat{X}\_{k + 1}
$$

Solution $x_k$ and $w_k$ to it is found using the linear smoother algorithm summarized [here]({{<ref "/content/posts/rts_as_optimization.md#summary-of-the-algorithm">}}).
Then the estimates are updated as:
$$
\hat{X}_k \leftarrow \hat{X}_k + \alpha x_k \\\\
\hat{W}_k \leftarrow \hat{W}_k + \alpha w_k
$$
Where $0 < \alpha \leq 1$ is selected as described in the previous section.

# Initialization

Numerical optimization methods are sensitive to the initial guess -- they search a local minimum in its vicinity.
In our problem we can get a good initial guess on $X_k$ by estimating them with an Extended Kalman Filter.
For the noises $W_k$ a natural guess is just zero.

# Termination

The iterations terminate if two criteria are satisfied

1. The relative reduction of the cost functions is less than $t_f$: 
   $$E(\hat{X}^{i}, \hat{W}^{i}) < (1 + t_f) E(\hat{X}^{i - 1}, \hat{W}^{i - 1}),$$
   where superscript denotes iteration index
2. Relative violation of the time transition equations is less than $t_c$ for each epoch $k$: 
   $$\left|\hat{X}_{k + 1} - f_k(\hat{X}_k, \hat{W}_k) \right| < t_c \max\left(\left|\hat{X}\_{k + 1}\right|, 1 \right),$$
   where operations and inequalities are applied elementwise

Parameters $t_f$ and $t_c$ are passed to the optimization subroutine.

# Error covariance estimation

Error covariance estimates for $X_k$ and $W_k$ are taken from the linear smoother output from the last iteration.

# Conclusion

I've presented a nonlinear batch estimation algorithm which solves a nonlinear optimization problem with a help of linear Kalman smoother for iteration subproblems.
The algorithm has solid theoretical foundations and gives reliable and easy to interpret estimates when converged.

# References

1. J. Nocedal, S. J. Wright <<Numerical Optimization, 2nd edition>
