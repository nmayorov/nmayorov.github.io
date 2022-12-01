---
title: "Derivation of Kalman smoother from an optimization perspective"
date: 2022-11-27
katex: true
---

In this post Kalman smoother formulas are derived as a solution to an optimization problem.

# Problem formulation

We consider an estimation problem for a discrete time linear stochastic system:
$$
\begin{gather*}
x_{k+1} = F_k x_k + G_k w_k \\\\
z_k = H_k x_k + v_k \\\\
\operatorname{E} x_0 = x_0^- \\\\
\operatorname{E} (x - x_0^-) (x - x_0^-)^T = P_0^- \succ 0 \\\\
\operatorname{E} w_k = 0 \\\\
\operatorname{E} w_k w_k^T = Q_k \succ 0 \\\\
\operatorname{E} v_k = 0 \\\\
\operatorname{E} v_k v_k^T = R_k \succ 0 \\\\
\operatorname{E} w_i w_j^T = 0 \text{ for } i \neq j \\\\
\operatorname{E} v_i v_j^T = 0 \text{ for } i \neq j \\\\
\operatorname{E} w_i v_j^T = 0 \\\\
\end{gather*}
$$

All symbols above are indexed by integer epoch and have the following meanings:

* $x_k$ -- state vector, shape $n \times 1$
* $F_k$ -- state transition matrix, shape $n \times n$
* $w_k$ -- process noise vector, shape $m \times 1$
* $G_k$ -- noise input matrix, shape $n \times m$
* $z_k$ -- measurement vector, shape $l \times 1$
* $H_k$ -- measurement matrix, shape $l \times n$
* $v_k$ -- measurement noise vector, shape $l \times 1$
* $x_0^-$ -- apriori state mean, shape $n \times 1$
* $P_0^-$ -- apriori state error covariance, shape $n \times n$
* $Q_k$ -- covariance of process noise,  shape $m \times m$
* $R_k$ -- covariance of measurement noise, shape $l \times l$

Additionally the following notation is used:

* $\operatorname{E} \ldots$ -- operation of taking mean
* $C \succ 0$ -- matrix $C$ is positive definite

The last 3 equations say that process and measurement noises are not time correlated and not correlated to each other.

The task is to find optimal estimates of $x_k$ for epochs $k = 0, 1, \ldots, N$ taking into account the system model and the observed measurements $z_0, z_1, \ldots, z_{N-1}$.

# Defining the optimization problem

We will determine optimal estimates of $x_k$ by minimizing a cost function.
The proposed cost function will be quadratic in $x_k$ and expressed as a sum of terms of 3 kinds:

1. A term for prior on $x_0$
2. Terms for measurements (associated with $z_k$)
3. Terms for the time propagation (relating $x_k$ and $x_{k+1}$)

## Case of square $G_k$

In the transition equation we can define the effective noise term as 
$$w_k^\prime = G_k w_k$$ 
$$\operatorname{E} w_k^\prime \left(w_k^\prime\right)^T = G_k Q_k G_k^T$$

If $G_k$ is square and invertible then $Q_k^\prime$ will be invertible too.
It means that each state is influenced by the noise and there is no deterministic part in the state propagation.
In this case the cost function can be defined naturally as follows:

$$
\begin{split}
J^\prime(x) =& \frac{1}{2} \left(x_0 - x_0^-\right)^T \left(P_0^-\right)^{-1} \left(x_0 - x_0^-\right) \\\\
            +& \frac{1}{2} \sum_{k = 0}^{N - 1} (H_k x_k - z_k)^T R_k^{-1} (H_k x_k - z_k) \\\\
            +& \frac{1}{2} \sum_{k = 0}^{N - 1} (F_k x_k - x_{k + 1})^T \left(Q_k^\prime\right)^{-1} (F_k x_k - x_{k + 1})
\end{split}
$$

Each term penalizes deviation of an actual measurement (or state) from its expectation weighed by inverses of corresponding covariance matrices.
It is also equal to the negative log probability of the joint distribution of all $x_k$ and $z_k$ assuming normal probability densities.

Minimization of $J^\prime(x)$ is a linear least-squares problem. 
The normal equation for it is a large positive definite tridiagonal system of $nN$ equations.
It can be solved efficiently by a two-pass algorithm in a time linear in $N$.
The resulting formulas can be interpreted as another form of a linear smoother (involving inverse covariance matrices).
Details can be found in [1].

This formulation is quite simple and has a straightforward solution.
Unfortunately the assumption of square $G_k$ and invertible $Q_k^\prime$ is limiting.
Problems where the number of states is larger than the number of noise sources are very common in practice.
There is no simple workaround for this issue and another approach needs to be developed.

## Case of rectangular $G_k$ with $m < n$

In this case the inverse of $Q_k^\prime$ can't be computed and thus the time propagation terms can't be formed directly.
To overcome this problem the noise vectors should be considered as unknown variables which need to be estimated along with $x$.
The following cost function is minimized:

$$
\begin{split}
J(x, w) =& \frac{1}{2} \left(x_0 - x_0^-\right)^T \left(P_0^-\right)^{-1} \left(x_0 - x_0^-\right) \\\\
        +& \frac{1}{2} \sum_{k = 0}^{N - 1} (H_k x_k - z_k)^T R_k^{-1} (H_k x_k - z_k) \\\\
        +& \frac{1}{2} \sum_{k = 0}^{N - 1} w_k^T Q_k^{-1} w_k
\end{split}
$$
Here $x$ and $w$ are not independent and the time propagation equations must be taken into account.
Doing this, we'll get the following constrained optimization problem:

$$
\min_{x, w} J(x, w) \text{ subject to} \space x_{k + 1} = F_k x_k + G_k w_k
$$

This is a linear least squares problem with linear equality constraints on the variables.
Such problem can be solved by a known linear algebra method.
The issue however is that this problem is large (in practically valuable cases) and sparse. 
Whether there are suitable general numerical solvers for this kind of sparse problems is a topic for another research.
In what follows we will derive an efficient solver specifically for this problem. 

# Lagrange function and its stationary points

In order to solve an equality constrained optimization problem one should form the Lagrange function as a sum of an original cost function and constraint equations multiplied by Lagrange multipliers.
Let's denote the multiplier for the time propagation equation from epoch $k$ to epoch $k + 1$ as $\lambda_{k+1}$.
The Lagrange function is then
$$
\begin{split}
L(x, w, \lambda) =& \frac{1}{2} \left(x_0 - x_0^-\right)^T \left(P_0^-\right)^{-1} \left(x_0 - x_0^-\right) \\\\
                 +& \frac{1}{2} \sum_{k = 0}^{N - 1} (H_k x_k - z_k)^T R_k^{-1} (H_k x_k - z_k) \\\\
                 +& \frac{1}{2} \sum_{k = 0}^{N - 1} w_k^T Q_k^{-1} w_k \\\\
                 +& \sum_{k = 0}^{N - 1} \lambda_{k + 1}^T (x_{k + 1} - F_k x_k - G_k w_k)
\end{split}
$$

Candidates for the minimum are searched as stationary points of the Lagrange function:
$$\frac{\partial L}{\partial x_k} = 0; \frac{\partial L}{\partial w_k} = 0; \frac{\partial L}{\partial \lambda_k} = 0$$
The partial derivative with respect to $x_k$ requires a separate treatment depending on $k$. 

For $k = 0$:
$$
H_0^T R_0^{-1} \left(H_0 x_0 - z_0\right) + \left(P_0^-\right)^{-1} \left(x_0 - x_0^-\right) - F_0^T \lambda_1 = 0
$$

For $0 < k < N$:
$$H_k^T R_k^{-1} \left(H_k x_k - z_k\right) + \lambda_k - F_k^T \lambda_{k + 1} = 0$$

For $k = N$:
$$\lambda_N = 0$$

To conform the first equation with the second one we introduce $$\lambda_0 \coloneqq \left(P_0^-\right)^{-1} \left(x_0 - x_0^-\right)$$

Then computing partial derivatives with respect to $w_k$ and $\lambda_k$ and joining all the equations together we get the following system:
$$
\begin{gather*}
x_{k + 1} = F_k x_k + G_k w_k \\\\
\lambda_k = F_k^T \lambda_{k + 1} + H_k^T R_k^{-1} \left(z_k - H_k x_k\right) \\\\
w_k = Q_k G_k^T \lambda_{k+1} \\\\
\lambda_0 = \left(P_0^-\right)^{-1} \left(x_0 - x_0^-\right) \\\\
\lambda_N = 0
\end{gather*}
$$

Let's analyze what is this system and how it can be solved.
First we can eliminate $w_k$ as it's directly expressed from $\lambda_{k+1}$.
Then we see that it's a system of difference equations for epochs $k = 0, 1, \ldots, N$ for vectors $x_k$ and $\lambda_k$ (both with $n$ elements).
To solve it, $2n$ initial or boundary conditions are required.
We have $n$ conditions at $k = 0$ (involving $\lambda_0$ and $x_0$) and $n$ conditions at $k = N$ (as $\lambda_N = 0$).
So this is a boundary value problem which potentially has a solution (if boundary conditions are consistent) and if it's the case we can find it.

# Kalman filter notation and formulas

To derive the solution we will use Kalman filter notation and formulas.
So here is a recap for these.
The following notation is used:

* $x_k^-, P_k^-$ -- a priori estimate and covariance at $k$, i.e. before processing measurement $z_k$
* $x_k^+, P_k^+$ -- a posteriori estimate and covariance at $k$, i.e. after processing measurement $z_k$

The correction equations (processing of $z_k$) are
$$
x_k^+ = x_k^- + K_k (z_k - H_k x_k^-) \\\\
P_k^+ = (I - K_k H_k) P_k^- \\\\
\text{ with } K_k = P_k^- H_k^T (H_k P_k H_k^T + R_k)^{-1}
$$

They can be expressed alternatively as
$$
(P_k^+)^{-1} = (P_k^-)^{-1} + H_k^T R_k^{-1} H_k \\\\
(P_k^+)^{-1} x_k^+ = (P_k^-)^{-1} x_k^- + H_k^T R_k^{-1} z_k
$$

And the prediction equations are
$$
x_{k+1}^- = F_k x_k^+ \\\\
P_{k+1}^- = F_k P_k^+ F_k^T + G_k Q_k G_k^T
$$


# Solving the difference equations

First rewrite the equations with $w_k$ eliminated:
$$
x_{k + 1} = F_k x_k + G_k Q_k G_k^T \lambda_{k+1} \\\\
\lambda_k = F_k^T \lambda_{k + 1} + H_k^T R_k^{-1} \left(z_k - H_k x_k\right) \\\\
\lambda_0 = \left(P_0^-\right)^{-1} \left(x_0 - x_0^-\right) \\\\
\lambda_N = 0
$$

Looking at the expression for $\lambda_0$ we make an assumption that
$$
x_k = x_k^- + P_k^- \lambda_k
$$
This assumption holds for $k = 0$ as the boundary condition.
Let's prove that if it's true for $k$ then it's true for $k + 1$ as well.
By induction principle then it will be true for any $k$.

Substituting $\lambda_k$ from the second equation into the assumption equation we get
$$
\begin{gather*}
x_k = x_k^- + P_k^- F_k^T \lambda_{k + 1} + P_k^- H_k^T R_k^{-1} z_k - P_k^- H_k^T R_k^{-1} H_k x_k \\\\
(I + P_k^- H_k^T R_k^{-1} H_k) x_k = x_k^- + P_k^- F_k^T \lambda_{k + 1} + P_k^- H_k^T R_k^{-1} z_k \\\\
P_k^-((P_k^-)^{-1} + H_k^T R_k^{-1} H_k) x_k = x_k^- + P_k^- F_k^T \lambda_{k + 1} + P_k^- H_k^T R_k^{-1} z_k \\\\
P_k^- (P_k^+)^{-1} x_k = x_k^- + P_k^- F_k^T \lambda_{k + 1} + P_k^- H_k^T R_k^{-1} z_k \\\\
x_k = P_k^+ \left( (P_k^-)^{-1} x_k^- + H_k^T R_k^{-1} z_k \right)  + P_k^+ F_k^T \lambda_{k + 1} \\\\
x_k = x_k^+ + P_k^+ F_k^T \lambda_{k + 1}
\end{gather*}
$$

Then substituting the new expression for $x_k$ into the time propagation equations we get
$$
\begin{gather*}
x_{k + 1} = F_k x_k^+ + (F_k^T P_k^+ F_k + G_k Q_k G_k^T) \lambda_{k + 1} \\\\
x_{k + 1} = x_{k + 1}^- + P_{k + 1}^- \lambda_{k + 1}
\end{gather*}
$$
Starting from the identity for $k$ we derived the same identity for $k + 1$ which proves the step of induction and with that the identity for any $k$.

Then consider the two derived identities together:
$$
\begin{gather*}
x_k = x_k^+ + P_k^+ F_k^T \lambda_{k + 1} \\\\
x_{k + 1} = x_{k + 1}^- + P_{k + 1}^- \lambda_{k + 1}
\end{gather*}
$$
Expressing $\lambda_{k + 1}$ from the second equation and substituting it into the first gets us
$$
\begin{gather*}
x_k^s = x_k^+ + C_k \left(x_{k + 1}^s - x_{k + 1}^-\right) \text{with} \space C_k = P_k^+ F_k^T \left(P_{k + 1}^-\right)^{-1}
\end{gather*}
$$
Here the superscript <<s>> denotes an optimal <<smoothed>> estimate.
This is a formula for the RTS smoother backward recursion [2].

From the boundary condition $\lambda_N = 0$ we get the starting point for this recursion $$x_{N}^s = x_{N}^-$$
It means that the optimal estimate at $k = N$ is equal to the filter estimate. 

# Error covariance computation

The errors of the estimates are defined as follows
$$
\Delta x_k^- = x_k^- - x_k \\\\
\Delta x_k^+ = x_k^+ - x_k \\\\
\Delta x_k^s = x_k^s - x_k \\\\
$$
with $x_k$ being the true state. The error covariances:
$$
\begin{gather*}
P_k^- = \operatorname{E} \Delta x_k^- \left(\Delta x_k^-\right)^T \\\\
P_k^+ = \operatorname{E} \Delta x_k^+ \left(\Delta x_k^+\right)^T \\\\
P_k^s = \operatorname{E} \Delta x_k^s \left(\Delta x_k^s\right)^T 
\end{gather*}
$$
The filter covariances $P_k^-$ and $P_k^+$ are known and the task is to determine $P_k^s$.

Turns out that derivation of the error covariance straight from the RTS recursion formula is difficult.
In fact I haven't seen such derivation anywhere. 
Instead authors rely on some indirect properties to take a shortcut to the result.
I will resort to that too with a solid justification.

## Kalman filter errors correlation property

Consider a state vector estimate before ($x^-$) and after ($x^+$) a vector $z$ is processed in a Kalman filter:
$$
x^+ = x^- + K (z - H x^-)
$$
For the errors we have
$$
\Delta x^+ = \Delta x^- + K (v - H \Delta x^-)
$$
The cross covariance is then
$$
\operatorname E \Delta x^+ (\Delta x^-)^T = (I - K H) P^- = P^+
$$
Using this property we get
$$
\operatorname{E}(\Delta x^+ - \Delta x^-)(\Delta x^+ - \Delta x^-)^T = P^- - P^+
$$
To derive the smoother covariance formula we need this property for pairs $\Delta x^s, \Delta x^+$ and $\Delta x^s, \Delta x^-$.
It holds for them because $x^s_k$ is a corrected version (in an optimal Kalman sense) of $x^-_k$ or $x^+_k$ using all measurements $z_i$ for $i \geq k$ or $i > k$ respectively.
This can be rigorously shown by forming an augmented state vector and considering the smoothing problem as a filtering problem for this augmented vector.

## Deriving the covariance recursion formula

Rewrite the equation for $x_k^s$ as
$$
x_k^s - x_k^+= C_k \left(x_{k + 1}^s - x_{k + 1}^-\right)
$$
From it the equivalent equation for the errors follows:
$$
\Delta x_k^s - \Delta x_k^+= C_k \left(\Delta x_{k + 1}^s - \Delta x_{k + 1}^-\right)
$$
Taking the covariance of both sides using the above stated properties we get:
$$
P_k^+ - P_k^s = C_k \left( P_{k+1}^- - P_{k+1}^s\right) C_k^T
$$
Or finally:
$$
P_k^s = P_k^+ + C_k \left(P_{k+1}^s - P_{k+1}^-\right) C_k^T
$$

# Estimation of noise vectors

Noise estimates are computed directly from $\lambda_{k+1}$ as
$$
w_k^s = Q_k G_k^T \lambda_{k+1}^s = Q_k G_k^T \left( P_{k+1}^- \right)^{-1}  \left(x_{k+1}^s - x_{k+1}^-\right)
$$
We are also interested in its error:
$$
\Delta w_k^s = w_k^s - w_k
$$
The covariance of $\Delta w_k^s$ can be computed using the same approach as for the state estimate error covariance.
The derivation is omitted for the sake of brevity and the final formula looks like this:

$$
Q_k^s \coloneqq \operatorname{E} \Delta w_k^s \left(\Delta w_k^s\right)^T = Q_k + B_k \left(P^s_{k + 1} - P^-_{k + 1} \right) B_k^T \\\\
\text{with } B_k = Q_k G_k^T \left(P\_{k+1}^-\right)^{-1}
$$

# Minor generalization of the model

When applying the linear smoother just derived to solve subproblems arising in nonlinear estimation a minor generalization of the model is required:

1. Consider a deterministic control signal $u_k$ in the time propagation equations
2. Consider that the process noise vectors have known nonzero mean values $w_k^-$

The Kalman prediction equation for the state then becomes
$$
x_{k + 1}^- = F_k x_k^+ + G_k w_k^- + u_k
$$

And the noise estimates are computed as
$$
w_k^s = w_k^- + Q_k G_k^T \left( P_{k+1}^- \right)^{-1}  \left(x_{k+1}^s - x_{k+1}^-\right)
$$

And everything else stays the same.

# Summary of the algorithm

To compute the state estimate and its covariance as the solution of the equality constrained linear least-squares problem the following algorithm must be used.

## Forward pass

The initial filter state and covariance estimates are given as $x_0^-$ and $P_0^-$.
Alternate correction and prediction steps for $k = 0, 1, \ldots, N - 1$.

### Correction step

$$
x_k^+ = x_k^- + K_k (z_k - H_k x_k^-) \\\\
P_k^+ = (I - K_k H_k) P_k^- \\\\
\text{with} \space K_k = P_k^- H_k^T \left(H_k P_k^- H_k^T + R_k \right)^{-1}
$$

### Prediction step

$$x_{k + 1}^- = F_k x_k^+ + G_k w_k^- + u_k$$
$$P_{k + 1}^- = F_k P_k^+ F_k^T + G_k Q_k G_k^T$$

All the variables $x_k^-, x_k^+, P_k^-, P_k^+$ are saved for the backward pass.

## Backward pass

Assign the smoothed estimate and covariance for epoch $N$ to the forward pass variables:
$$
x_N^s = x_N^- \\\\
P_N^s = P_N^-
$$

For $k = N - 1, N - 2, \ldots, 0$ apply the following formulas:
$$x_k^s = x_k^+ + C_k (x_{k+1}^s - x_{k+1}^-)$$
$$w_k^s = w_k^- + B_k (x_{k+1}^s - x_{k+1}^-)$$
$$P_k^s = P_k^+ + C_k (P_{k+1}^s - P_{k+1}^-) C_k^T$$
$$Q_k^s = Q_k + B_k (P_{k+1}^s - P_{k+1}^-) B_k^T$$
$$\text{with } C_k = P_k^+ F_k^T \left(P_{k + 1}^-\right)^{-1} \text{ and } B_k = Q_k G_k^T \left( P_{k+1}^- \right)^{-1}$$

In the end, the variables $x_k^s, P_k^s$ and $w_k^s, Q_k^s$ contain optimal estimates and error covariances of state and noise vectors respectively for $k = 0, 1, \ldots, N$.

# Conclusion

A step-by-step derivation of Kalman smoother (RTS form [2]) as a solution to an optimization problem was presented.
Expressing this state estimation algorithm as a solution to the optimization problem opens up opportunities to develop other advanced estimation methods (nonlinear, robust, etc.) based on an optimization approach.

# References

1. [Optimization viewpoint on Kalman smoothing, with applications to robust and sparse estimation](https://arxiv.org/abs/1303.1993)
2. H. E. Rauch, F. Tung, C. T. Striebel <<Maximum Likelihood Estimates of Linear Dynamic Systems>>
