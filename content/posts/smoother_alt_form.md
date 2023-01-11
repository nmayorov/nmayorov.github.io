---
title: "Alternative form of Kalman smoother"
date: 2023-01-05
katex: true
---

In [this post]({{<ref "/content/posts/rts_as_optimization.md">}}) (recommended to read beforehand) I've shown a derivation of the Kalman smoother as a solution to an optimization problem.
The resulting formulas are surprisingly elegant, however their applicability depends on the assumption that apriori filter covariance matrices $P_k^-$ are positive definite and invertible.
This assumption might be limiting in practical problems and thus another form of the Kalman smoother is derived here.

# Motivating example

A singular covariance matrix may arise in the following practical scenario.
Imagine that we want to process measurements which relate states at the current and previous epochs -- $x_k$ and $x_{k-1}$.
A possible example might be processing of distance increments from an odometeter in a navigation algorithm.
Such measurements are not directly supported by Kalman filter or smoother algorithms.
However we can cast them into the Kalman framework by considering an augmented vector
$$
x^a_{k + 1} = \begin{bmatrix}
x_{k + 1} \\\\
x_k
\end{bmatrix}
$$
The approach is also known as <<stochastic cloning>>.

Now consider a time transition from $x_k$ to $x^a_{k + 1}$:
$$
x^a_{k + 1} = F^a_k x_k + G^a_k w_k \\\\
\text{with } F^a_k = \begin{bmatrix}
F_k \\\\\
I
\end{bmatrix}, G_k^a = \begin{bmatrix}
G_k \\\\
0
\end{bmatrix}
$$
For the covariance we then get
$$
P^a_{k + 1} = F_k^a P_k (F_k^a)^T + G^a_k Q_k (G^a_k)^T = \begin{bmatrix}
F_k P_k F_k^T + G_k Q_k G_k^T & F_k P_k \\\\
P_k F_k^T & P_k
\end{bmatrix}
$$
If $G_k Q_k G_k^T$ is singular, then the covariance $P^a_{k + 1}$ will be singular too.
For example if there is no noise we get
$$
P^a_{k + 1} = \begin{bmatrix}
F_k P_k F_k^T & F_k P_k \\\\
P_k F_k^T & P_k
\end{bmatrix}
$$
Multiplying the second block row by $F_k$ we get identical block rows, which shows that $P^a_{k + 1}$ is singular.

Also even if we don't use stochastic cloning, the covariance matrix may gradually approach singularity from the perspective of floating point numbers as measurements are processed.
In this case we want to avoid computing its inverse, which is another solid motivation to develop an alternative form where inverses are not needed.

# Problem formulation and the boundary value problem

Let's recap the optimization problem.

We want to make a small generalization that the initial apriori covariance $P_0^-$ might be singular.
In this case we can't add a term $1/2 (x_0 - x_0^-)^T (P_0^-)^{-1} (x_0 - x_0^-)$ to the cost function.
But any covariance matrix can be represented as 
$$
P_0 = G_{-1} Q_{-1} G_{-1}^T
$$
where $Q_{-1}$ is positive definite (subscripts -1 is purely a notational choice).
Keeping in mind such factorization we express the initial state $x_0$ as 
$$
x_0 = x_0^- + G_{-1} w_{-1}
$$
where noise $w_{-1}$ is assumed to have covariance $Q_{-1}$.

The state estimates are searched by considering the cost function:
$$
\begin{split}
J(x, w) =& \frac{1}{2} w_{-1}^T Q_{-1}^{-1} w_{-1} \\\\
        +& \frac{1}{2} \sum_{k = 0}^{N - 1} (H_k x_k - z_k)^T R_k^{-1} (H_k x_k - z_k) \\\\
        +& \frac{1}{2} \sum_{k = 0}^{N - 1} w_k^T Q_k^{-1} w_k
\end{split}
$$
Which is minimized subject to the state transition equations:
$$
\min_{x, w} J(x, w) \text{ subject to } \\\\
x_0 = x_0^- + G_{-1} w_{-1} \\\\ 
x_{k + 1} = F_k x_k + G_k w_k
$$

The problem is solved by finding the stationary point of the Lagrange function.
It gives the following discrete-time boundary value problem:
$$
x_{k + 1} = F_k x_k + G_k Q_k G_k^T \lambda_{k+1} \\\\
\lambda_k = F_k^T \lambda_{k + 1} + H_k^T R_k^{-1} \left(z_k - H_k x_k\right) \\\\
P_0^- \lambda_0 = x_0 - x_0^- \\\\
\lambda_N = 0
$$
Where $\lambda_k$ are introduced Lagrange multipliers. From them the noise vectors are computed as
$$
w_k = Q_k G_k^T \lambda_{k+1} \\\\
$$

Note that our approach to consider the singular initial covariance came down to writing the boundary condition $\lambda_0 = \left(P_0^-\right)^{-1} \left(x_0 - x_0^-\right)$ as $P_0^- \lambda_0 = x_0 - x_0^-$.
The rigorous derivation gave the only plausible result.
Also note that at this point there is no notion of covariance matrices (besides $P_0$), let alone their positive definitiveness.

# Some matrix identities related to Kalman filtering

As a starting point we take the standard Kalman correction formula with the optimal gain matrix $K$:
$$
K = P^- H^T (H P^- H^T + R)^{-1} \\\\
P^+ = (I - K H) P^-
$$
The formula is applicable for any covariance matrix $P^-$ as long as $R$ is positive definite, which we already assume in our problem.

The last equation can be written as
$$
P^- - P^+ = K H P^- = P^- H^T K^T
$$

Let's combine it with the gain equation as follows:
$$
K (H P^- H^T + R) = P^- H^T \\\\
K H P^- H^T = P^- H^T - P^+ H^T
$$
Subtracting the equations we get the following alternative expression for $K$:
$$
K = P^+ H^T R^{-1}
$$

Using this alternative expression we can also write
$$
P^- - P^+ = P^+ H^T R^{-1} H P^- = P^- H^T R^{-1} H P^+
$$

Now let's prove the identity:
$$
(I + H^T R^{-1} H P^-)^{-1} = I - H^T R^{-1} H P^+ = I - H^T K^T
$$
To do that we show that the multiplications on both sides give the identity matrix:
$$
\begin{split}
(I + H^T R^{-1} H P^-) (I - H^T R^{-1} H P^+) = I + H^T R^{-1} H (P^- - P^+ - P^- H^T R^{-1} H P^+) = I \\\\
(I - H^T R^{-1} H P^+) (I + H^T R^{-1} H P^-) = I + H^T R^{-1} H (P^- - P^+ - P^+ H^T R^{-1} H P^-) = I
\end{split}
$$

By taking its transpose we get
$$
(I + P^- H^T R^{-1} H)^{-1} = I - P^+ H^T R^{-1} H = I - K H
$$

Similarly another formula can be proven
$$
(R + H P^- H^T)^{-1} = R^{-1} - R^{-1} H P^+ H^T R^{-1}
$$

Equipped with the derived identities we can carry on with solving the boundary-value problem.

# Equivalent problem

First we prove the equivalence of the problems:
$$
x_{k + 1} = F_k x_k + G_k Q_k G_k^T \lambda_{k+1} \\\\
\lambda_k = F_k^T \lambda_{k + 1} + H_k^T R_k^{-1} \left(z_k - H_k x_k\right) \\\\
P_0^- \lambda_0 = x_0 - x_0^- \\\\
\lambda_N = 0 \\\\
\downdownarrows \upuparrows \\\\
x_k = x_k^- + P_k^- \lambda_k \\\\
\lambda_k = F_k^T \lambda_{k + 1} + H_k^T R_k^{-1} \left(z_k - H_k x_k\right) \\\\
\lambda_N = 0
$$
Meaning that any solution to the first problem is the solution to the second problem and vice versa.
Or that the set of solutions are identical.

For convenience let's call $x_k = x_k^- + P_k^- \lambda_k$ the <<magic equation>>.

## Proof from up to down

We need to prove that the magic equation follows from the difference equations.

We do it by induction.
For $k = 0$ it's true as the boundary condition.
Now assume that it's true for some $k$ and substitute the expression for $\lambda_k$ into it:
$$
x_k = x_k^- + P_k^- F_k^T \lambda_{k + 1} + P_k^- H_k^T R_k^{-1} z_k - P_k^- H_k^T R_k^{-1} H_k x_k \\\\
(I + P_k^- H_k^T R_k^{-1} H_k) x_k = x_k^- + P_k^- F_k^T \lambda_{k + 1} + P_k^- H_k^T R_k^{-1} z_k \\\\
$$
Now using the inverse formula for the matrix multiplier on the left we get
$$
x_k = (I - K_k H_k) x_k^- + P_k^+ F_k^T \lambda_{k + 1} + P_k^+ H_k^T R_k^{-1} z_k
$$
Here we spot the expression for $K_k$ before $z_k$ and rewrite the equation as
$$
x_k = x_k^- + K_k (z_k - H_k x_k^-) + P^+_k F_k^T \lambda\_{k+1} = x_k^+ + P_k^+ F_k^T \lambda\_{k+1}
$$
And finally substitute this expression into recursion for $x$:
$$
x\_{k + 1} = F_k x_k^+ + (F_k^T P_k^+ F_k + G_k Q_k G_k^T) \lambda\_{k + 1} \\\\
x\_{k + 1} = x\_{k + 1}^- + P\_{k + 1}^- \lambda\_{k + 1}
$$
Which proves the step of the induction and with that the equation for all $k$.

## Proof from down to up

We need to prove that the difference equation for $x$ follows from the magic equation and the difference equation for $\lambda$.

Substitute the expression for $\lambda_k$ into the magic equation to get (derived above already)
$$
x_k = x_k^+ + P_k^+ F_k^T \lambda\_{k+1} \\\\
x_k^+ = x_k - P_k^+ F_k^T \lambda\_{k+1}
$$
Now substitute $x_k^+$ into the magic equation for $k + 1$:
$$
\begin{split}
x_{k + 1} = x_{k + 1}^- + P_{k + 1}^- \lambda\_{k + 1} = F_k x_k^+ + P^-_{k + 1} \lambda\_{k + 1} 
= F_k x_k + (P\_{k + 1}^- - F_k P_k^+ F_k^T) \lambda\_{k + 1} = \\\\  =  F_k x_k + G_k Q_k G_k^T \lambda\_{k + 1}
\end{split}
$$

# Solution to the problem

With the proof from the previous section we can solve the more simple equivalent problem
$$
x_k = x_k^- + P_k^- \lambda_k  \\\\
\lambda_k = F_k^T \lambda_{k + 1} + H_k^T R_k^{-1} \left(z_k - H_k x_k\right) \\\\
\lambda_N = 0
$$
and be sure to get the solution of the original problem.

## RTS form

To derive the RTS form we consider the already derived corollary equation along with the magic equation:
$$
x_k = x_k^- + P\_k^- \lambda_k \\\\
x_k = x_k^+ + P_k^+ F_k^T \lambda\_{k+1}
$$
If all $P_k^-$ are invertible the difference equation for $\lambda$ follows from them:
$$
x_k^- + P_k^- \lambda_k = x_k^+ + P_k^+ F_k^T \lambda\_{k+1} \\\\
\lambda_k = (P_k^-)^{-1} P_k^+ F_k^T \lambda\_{k+1} + (P_k^-)^{-1} (x_k^+ - x_k^-) \\\\
\lambda_k = F_k^T \lambda_{k + 1} + H_k^T R_k^{-1} \left(z_k - H_k x_k\right) 
$$
It means the two equations are equivalent to the original problem.

Rewriting them slightly different
$$
x_{k + 1} = x_{k + 1}^- + P\_k^- \lambda_{k + 1} \\\\
x_k = x_k^+ + P_k^+ F_k^T \lambda\_{k+1}
$$
we can eliminate $\lambda_{k + 1}$ and come up with the unique solution:
$$
x_N^s = x_N^- \\\\
x_k^s = x_k^+ + P_k^+ F_k^T \left(P_{k + 1}^-\right)^{-1} \left(x_{k + 1}^s - x_{k + 1}^-\right)
$$
and
$$
\lambda_k^s = (P_k^-)^{-1} (x_k^s - x_k^-) \\\\
w_k^s = Q_k G_k^T \lambda_{k+1}^s
$$

## Form with explicit recursion for $\lambda$

If $P_k^-$ is singular then the two equations from the RTS derivation are not equivalent to the original problem.
Meaning that there is no unique solution to them and an arbitrary solution is not guaranteed to be the solution of the original problem.

We go back to the original equation for $\lambda_k$ and substitute $x_k = x_k^- + P_k^- \lambda_k$ into it:
$$
\lambda_k = F_k^T \lambda_{k+1} + H_k^T R_k^{-1} (z_ k - H_k x_k^- - H_k P_k^- \lambda_k) \\\\
(I + H_k^T R_k^{-1} H_k P_k^-) \lambda_k = F_k^T \lambda_{k+1} + H_k^T R_k^{-1} (z_k - H_k x_k^-)
$$
Taking the inverse of the matrix multiplier on the left we get
$$
\lambda_k = (I - K_k^T H_k^T) F_k^T \lambda_{k + 1} + (I - H_k^T R_k^{-1} H_k P^+_k) H_k^T R_k^{-1} (z_k - H_k x_k^-) \\\\
\lambda_k = (I - K_k^T H_k^T) F_k^T \lambda\_{k + 1} + H_k^T (H_k P^-_k H_k^T + R_k)^{-1} (z_k - H_k x_k^-) \\\\
$$
Introduce some notation:

* $\Psi_k = F_k (I - K_k H_k)$ -- can be though of as a total matrix for the state propagation, combing filter correction and prediction steps
* $S_k = H_k P^-_k H_k^T + R_k$ -- covariance matrix of the innovation vector $z_k - H_k x_k^-$

Using it the recursion can be written as
$$
\lambda_N^s = 0 \\\\
\lambda_k^s = \Psi_k^T \lambda^s\_{k + 1} + H_k^T S_k^{-1} (z_k - H_k x_k^-)
$$

From $\lambda^s_k$ we can compute optimal state and noise estimates as 
$$
x_k^s = x_k^- + P_k^- \lambda^s_k = x_k^+ + P_k^+ F_k^T \lambda_{k + 1}^s  \\\\
w_k^s = Q_k G_k^T \lambda^s_{k + 1}
$$ 

## Uniqueness of the solution

We see that the solution is unique, because $\lambda_k^s$ obeys a deterministic recursion process starting from $\lambda_N^s = 0$.
And $x_k^s$ and $w_k^s$ are uniquely determined from $\lambda_k^s$.
This fact doesn't depend on whether $P_k^-$ is singular or not.

One thing that may seem a bit confusing is that when $P_k^-$ is singular the equation $x_k = x_k^- + P_k^- \lambda_k$ is satisfied for an infinite set of $\lambda_k$.
However $\lambda_k$ determine the noise vectors $w_k$ which purpose is to make state transition equations consistent.
Thus arbitrary adjustment of individual $\lambda_k$ is not possible, they are determined by the recursion equation uniquely.

# Covariance computation

To compute covariance of the estimate errors we recursively update covariance of $\lambda_k^s$ as:
$$
\Lambda^s_k = \Psi^T_k \Lambda^s_{k + 1} \Psi_k + H_k^T S_k^{-1} H_k \\\\
\Lambda^s_N = 0
$$
Then using the same reasoning as for the RTS form we get the following expressions for error covariances:
$$
P_k^s = P_k^- - P_k^- \Lambda^s_k P_k^- = P_k^+ - P_k^+ F_k^T \Lambda^s_{k + 1} F_k P_k^+ \\\\
Q_k^s = Q_k - Q_k G_k^T \Lambda^s_{k + 1} G_k Q_k
$$

# Algorithmic aspects: multiple or missing measurements
In practical algorithms its better to assume that at each epoch there is an arbitrary number of independent measurements (possibly zero).
In this regard its better to split $\lambda_k^s$ and $\Lambda_k^s$ computation into (backward) <<prediction>> and <<correction>> steps similar to Kalman filter.

The prediction step:
$$
\lambda_{k}^- = F_k^T \lambda_{k + 1}^+ \\\\
\Lambda_{k}^- = F_k^T \Lambda_{k + 1}^+ F_k
$$
The correction step:
$$
\lambda_{k}^+ = (I - K_k H_k)^T \lambda_{k}^- + H_k^T S_k^{-1} (z_k - H_k x_k^-) \\\\
\Lambda_{k}^+ = (I - K_k H_k)^T \Lambda_{k}^- (I - K_k H_k) + H_k^T S_k^{-1} H_k
$$
It is repeated for each measurement, applying the formulas to the current estimate. 
The gain vectors, innovations and their covariances are saved from the filter pass.
The crucial point is that the formulas are applied in the *reverse* to the order in which the measurements were processed in the filter.
In the end the result is the <<smoothed>> estimate $\lambda_k^s$ which is used to compute $x_k^s$ and $w_k^s$.

The correctness of such approach can be understood by considering intermediate states between the measurements.

# Conclusion

An alternative form of the Kalman smoother which doesn't assume or require positive definite covariance matrices was derived.
This form is somewhat more involved to implement.
However for practical algorithms the main requirements are robustness and universality and this alternative algorithm seems to fullfil them better.
