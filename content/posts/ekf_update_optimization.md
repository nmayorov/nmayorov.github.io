---
title: "Extended Kalman Filter update step as an optimization problem"
date: 2023-05-04
katex: true
---

The Extended Kalman Filter is classically built as an extension of the linear Kalman filter using system and measurement models linearization.
In this regard the update step in EKF is naturally done in a single step.
The iterated EKF aims to improve the linearization point by doing several update iterations and recomputing measurement Jacobian each time.
It's known to be connected with nonlinear optimization.

In this note I want to derive possible EKF update strategies starting from the optimization viewpoint.

# The update step in the linear Kalman filter

In the update step a prior distribution on the state vector $x \sim \mathcal{N}(x^-, P^-)$ is optimally combined with a measurement $z = H x + v$ with $v \sim \mathcal{N}(0, R)$.
To do so the following optimization problem is formulated:
$$
\min_x J(x) = \frac{1}{2} (x - x^-)^T (P^-)^{-1} (x - x^-) + \frac{1}{2} (H x - z)^T R^{-1} (H x - z)
$$
This is a linear least-squares problem -- the quadratic function has a global minimum where its gradient equals zero:
$$
(P^-)^{-1} (x - x^-) + H^T R^{-1} (H x - z) = 0
$$
This equation can be transformed as follows (with $x$ denoted as $x^+$):
$$
\left((P^-)^{-1} + H^T R^{-1} H \right) x^+ = (P^-)^{-1} x^- + H^T R^{-1} z
$$
And this is an alternative form of the Kalman correction formula equivalent to
$$
x^+ = x^- + K (z - H x^-) \text{ with } K = P^- H^T (H P^- H^T + R)^{-1}
$$

# The update step in the Extended Kalman Filter 

In a similar fashion, a prior distribution $X \sim \mathcal{N}(X^-, P^-)$ is optimally combined with a nonlinear measurement $Z = h(X) + v$ with $v \sim \mathcal{N}(0, R)$, by forming an optimization problem:
$$
\min_X E(X) = \frac{1}{2} (X - X^-)^T (P^-)^{-1} (X - X^-) + \frac{1}{2} (h(X) - Z)^T R^{-1} (h(X) - Z)
$$
This is a nonlinear least-squares problem which can be solved iteratively by forming a sequence of linear subproblems using linearization around the current estimate.

Let $X_i$ be an estimate after $i$ iterations with $X_0 \coloneqq X^-$. 
An update $x_i$ is sought by substituting $X = X_i + x$ into $E(X)$ and keeping terms linear in $x_i$.
It results in the following linear least-squares problem:
$$
\min_x J(x) = \frac{1}{2} (x - x_i^-)^T (P^-)^{-1} (x - x_i^-) + \frac{1}{2} (H_i x - z_i)^T R^{-1} (H_i x - z)
$$
Where the following variables were introduced:
$$
x^- = X^- - X_i \\\\
z_i = Z - h(X_i) \\\\
H_i = \left.\dfrac{\partial h(X)}{\partial X}\right\vert\_{X_i}
$$
The linear subproblem has the form as in the linear Kalman filter update step for which the solution is known:
$$
x_i = x_i^- + K_i (z_i - H_i x_i^-) \text{ with } {K_i = P^- H_i^T (H_i P^- H_i^T + R)^{-1}}
$$

The estimate is then updated as
$$
X_{i + 1} = X_i + \alpha_i x_i
$$
Where $0 < \alpha_i \leq 1$ is selected to ensure the convergence using a procedure known as a line search.
The simplest strategy is to set $\alpha_i = 1$ which works well when $X_i$ is close to the minimizer.

In the iterated EKF the full corrections are applied with $\alpha_i = 1$:
$$
X_{i + 1} = X_i + x_i = X_i + X^- - X_i + K_i (Z - h(X_i) - H_i (X^- - X_i)) \approx X^- + K_i (Z - h(X^-))
$$
The last approximate equality follows from the first-order Taylor expansion of $h(X)$.
The whole procedure is already based on this and thus such substitution is justifiable.
If implemented using the last equation, the iterated EKF update comes down to updating the point for Jacobian computation $H_i$.
It is generally known that a correct Jacobian value is the most important for consistent EKF update and so this scheme makes practical and conceptual sense.

The standard EKF update is essentially a single step of the optimization procedure taken with $\alpha_0 = 1$.
It is generally known from optimization theory that such step in <<good>> conditions can yield estimate very close to the optimum.
This can serve as a conceptual explanation why the standard EKF update often works well in practice.

# Discussion

There are 3 choices when implemented the EKF update:

1. Standard EKF update -- the most simple and computationally efficient
2. Iterated EKF update using one of two slightly different numerical schemes -- quite simple but at least doubles required computations (because at least 2 iterations are required to check for the convergence)
3. Proper optimization with the line search (or other means to ensure convergence) -- the most rigorous, but requires understanding of optimization algorithms and some tuning. 
   At least doubles required computations

In approaches 2 and 3 a care of proper convergence control must be taken to avoid infinite or excessively large amount of iterations.

I believe that the right approach depends on the application (the estimation problem at hand) and should be selected empirically.
I personally have the most experience with the approach 1 and know that it works well for many practical problems with mild measurement nonlinearities.
After all, it was successfully used in many practical systems as the basic EKF is a <<go-to>> algorithm for the estimation problems.

More advanced approaches 2 and 3 are worth considering, especially if the computation cost is not a concern.
They can't degrade estimation performance, but might alleviate EKF deficiencies in difficult cases.  

# Further leveraging the optimization approach

The optimization approach is powerful because it naturally allows for modification of the cost function.
For example, imagine that the measurement noise instead of normal distribution has Laplace distribution.
Empirically it means that it is more heavy-tailed with more likely outliers.
This situation can be handled by considering a proper loss function (can be thought as negative logarithm of probability):
$$
E(X) = \frac{1}{2} (X - X^-)^T (P^-)^{-1} (X - X^-) + \sqrt{2} \sum_{j = 0}^{N - 1} \frac{|h_j(X) - Z_j|}{\sigma_j}
$$
Here $j$ denotes vector index and $\sigma_j$ is a standard deviation of independent noise components.
This optimization problem is considerably more difficult, but it is well defined and can be solved by a proper optimization algorithm.

Other cost functions might be proposed to account for required empirical or theoretical properties of the measurements.

# Conclusion

The update step of the Extended Kalman Filter was considered as an optimization problem.
Iterated and standard EKF update schemes are shown to be simplified approaches to solve this problem.
The optimization viewpoint is advantageous as it allows for modification of the cost function to account for different measurement model properties.
