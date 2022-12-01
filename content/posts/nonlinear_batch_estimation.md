---
title: "Nonlinear batch estimation"
date: 2022-11-30
katex: true
draft: true
---

In this note I present a development of nonlinear batch estimation algorithm.

# Model description

The state estimation problem in a nonlinear system is considered.
The formulations is analogous to the [linear case]({{<ref "/content/posts/rts_as_optimization.md#problem-formulation">}}), but with nonlinear transition and measurement equations.
We use uppercase letters to denote variables participating in the nonlinear model. 
Time transition and measurement equations for which are
$$
Z_k = f_k(X_k) + V_k \\\\
X_{k+1} = h_k(X_k, W_k)
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
The method we will apply differs in how we form the Hessian approximation and might be called a Gauss-Newton approximation of the SQP method.

In a nonlinear least-squares problem the cost function has the form
$$
f(x) = \frac{1}{2} r(x)^T r(x)
$$
The Jacobian matrix of $r(x)$ is defined as
$$
J(x) = \begin{bmatrix} \nabla_x r_1(x) & \nabla_x r_2(x) & \ldots & \nabla_x r_n(x) \end{bmatrix}^T
$$
Then the gradient and Gauss-Newton approximation of the Jacobian are:
$$
\nabla_x f(x) = J^T r(x) \\\\
\nabla^2_{xx} f(x) \approx J(x)^T J(x)
$$
Opposed to SQP method we also don't include the Hessian of the constraints, i. e. form a quadratic model for the cost function, not for the Lagrangian.

With these approximations the quadratic subproblem is
$$
\min_{p_i} (r_i + J_i p_i)^T (r_i + J_i p_i) \text{ subject to } c_i + A_i p_i = 0 \\\\
\text{with } r_i = r(x_i), J_i = J(x_i), c_i = c(x_i), A_i = \nabla_x c(x) \vert_{x_i}
$$
Logic of the algorithm iterations remain the same.

It must be said that the proposed method is not widely known. 
However my intuition and limited experience tells me that it is a legit approach akin to Gauss-Newton method of solving unconstrained nonlinear least squares.

For our specific problem the linear subproblem formulated above is solved by a Kalman smoother.
The details will be presented in the next section.

# Constructing the linearized subproblem

First introduce convenient notation. Let 

* $\hat{X}_k, \hat{W}_k$ be current estimates of $X_k$ and $W_k$
* $x_k, w_k$ be correction to the current estimates

To conform to the estimation reasoning like in Extended Kalman Filter we define relation for the linearization as 
$$
\hat{X}_k = X_k + x_k \\\\
\hat{W}_k = W_k + w_k
$$
That is the current estimate equals the true value plus the correction or error to be estimated.
This form seems more natural for the estimation problems.

The linearized version of the measurement residual:
$$
Z_k - h(X_k) = Z_k - h(\hat{X}_k - x_k) \approx Z_k - h(\hat{X}_k) + H_k x_k = H_k x_k - z_k \\\\
\text{with } H_k = \left.\frac{\partial h(X)}{\partial X}\right\vert\_{\hat{X}_k} \text{ and } z_k = h(\hat{X}_k) - Z_k
$$

The linearized version of the time propagation equation:
$$
X_{k + 1} = f_k(X_k, W_k) \\\\
\hat{X}_{k + 1} - x\_{k + 1} = f_k(\hat{X}_k - x_k, \hat{W}_k - w_k) \approx f_k(\hat{X}_k, \hat{W}_k) - F_k x_k - G_k w_k \\\\
x\_{k + 1} = F_k x_k + G_k w_k + u_k \\\\
\text{with } F_k = \left.\frac{\partial f(X, W)}{\partial X}\right\vert\_{\hat{X}_k, \hat{W}_k}
, G_k = \left.\frac{\partial f(X, W)}{\partial W}\right\vert\_{\hat{X}_k, \hat{W}_k}
, u_k = \hat{X}\_{k + 1} - f_k(\hat{X}_k, \hat{W}_k)
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
x_0^- = \hat{X}_0 - X_0^- \\\\
z_k = h(\hat{X}_k) - Z_k \\\\
H_k = \left.\frac{\partial h(X)}{\partial X}\right\vert\_{\hat{X}_k} \\\\
w_k^- = \hat{W}_k \\\\
F_k = \left.\frac{\partial f(X, W)}{\partial X}\right\vert\_{\hat{X}_k, \hat{W}_k} \\\\
G_k = \left.\frac{\partial f(X, W)}{\partial W}\right\vert\_{\hat{X}_k, \hat{W}_k} \\\\
u_k = \hat{X}\_{k + 1} - f_k(\hat{X}_k, \hat{W}_k)
$$

Solution $x_k$ and $w_k$ to it is found using the linear smoother algorithm summarized [here]({{<ref "/content/posts/rts_as_optimization.md#summary-of-the-algorithm">}}).

# Step selection for the estimate update

After the solution $x_k$ and $w_k$ to the subproblem is obtained the estimates are updated as:
$$
\hat{X}_k \leftarrow \hat{X}_k - \alpha x_k \\\\
\hat{W}_k \leftarrow \hat{W}_k - \alpha w_k
$$
Where $0 < \alpha \leq 1$ is selected such to decrease the cost function and violation of the constraints.

The idea is to monitor a merit function which combines the cost function with constraint violation $l_1$-norm [1], section 11.2:
$$
\Phi(X, W; \mu) = E(X, W) + \mu \sum_{k = 0}^{N - 1} \left\lVert X_{k + 1} - f_k(X_k, W_k) \right\rVert_1
$$
A proper value for $\mu$ is generally unknown and must be adjusted during iterations (increasing only).
The whole theory is quite involved and requires experimentation.
For now I just say that we adjust $\mu$ and select $\alpha$ on each iteration such as to drive constraints violations to zero while improving the cost function (at a certain level of constraints violation).
The simplest strategy of selecting $\alpha = 1$ might also be effective because a good initial guess is available.

# Initialization and termination

Numerical optimization methods are sensitive to the initial guess -- they search a local minimum in its vicinity.
In our problem we can get a good initial guess on $X_k$ by estimating them with an Extended Kalman Filter.
For the noises $W_k$ a natural guess is just zero.

There are several possible criteria for stopping the iterations:

1. The cost function decrease on the last iteration is small compared to the cost function itself
2. The correction of $X_k$ and $W_k$ on the last iteration is small compared to their norm
3. The norm of the of the gradient of $E$ is small

It might be necessary to combine criteria 1, 2, 3 with a check on the constraints violation norm
Selecting a proper criteria requires experimentation as well.

# Conclusion

I've presented an outline of the nonlinear batch estimation algorithm which uses linear Kalman smoother on each iteration.
This is a preliminarily version and I will fill in missing details as I get practical experience with it.

# References

1. J. Nocedal, S. J. Wright <<Numerical Optimization, 2nd edition>
