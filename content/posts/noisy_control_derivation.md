---
title: 'Rigorous consideration of noisy control input in Kalman filter state propagation'
date: 2024-12-20
katex: true
---

In a [previous post]({{<ref "/posts/time_propagation_noise">}}) I've considered an intuitive engineering approach of handling noisy control input signals in Extended Kalman Filters.
Here I want to develop a more rigorous view to this approach.

# Problem formulation

Let's consider a linear state transition equation for $x_k$ of the form:
$$
x_{k + 1} = F_k x_k + B_k u_k,
$$
where $u_k$ is the deterministic sequence of control vectors.
Initial distribution of $x$ is given as normal:
$$
x_0 \sim \mathcal{N}(\hat{x}_0, P_0)
$$

Now imagine that only noisy measurements $\tilde{u}_k$ of control signals are available (think IMU readings in INS).
To be specific let's assume
$$
\tilde{u}_k = u_k + w_k, \text{with } w_k \sim \mathcal{N}(0, Q_k)
$$

Quite naturally the prediction step of a Kalman filter for this formulation has the form:
$$
\hat{x}_{k + 1} = F_k \hat{x}_k + B_k \tilde{u}_k \\\\
P\_{k + 1} = F_k P_k F_k^T + B_k Q_k B_k^T
$$
Where $\hat{x}_k$ and $P_k$ are conditional mean and covariance, the mean being the optimal estimate by all criteria.

It makes sense intuitively and was empirically explained for an error state of a nonlinear problem in my post linked above.
However, how do we exactly derive these equations in a rigorous way similarly to [1]?
How do we interpret the result?
Are these the same conditional mean and covariance arising from Bayesian approach?

# Conceptual solution

The easiest satisfactory solution I've found is to consider a control input vector as a part of an augmented state vector.
Then we can proceed with the following reasoning:

1. The availability of a measured control $u_k$ is nothing but a *measurement*, which we already know how to handle
2. The fact that we use these measurements means that our estimated mean and covariance are *conditioned* on these noisy control measurements
3. The fact that we apply the same steps -- measurement and prediction, allows for exactly the same interpretation of linear Kalman filter as derived in [1]

Before moving on to the detailed solution let's consider a supplementary problem.

# State initialization from complete uncertainty by a measurement

Consider a situation where we start from a complete uncertainty for the state vector $x$ and receive a measurement of the form:
$$
z = x + v, v \sim \mathcal{N}(0, R)
$$

By the Bayes rule we can compute a posterior distribution:
$$
p(x | z) \propto p(z | x) p (x) \propto p(z | x) = \mathcal{N}(z | x, R) = \mathcal{N}(x | z, R)
$$
Here we've used that $p(x)$ is constant (i.e. complete uncertainty) and that the normal distribution is symmetric with respect to the argument and the mean.
The result holds for other noise distributions as long as its PDF depends only on the absolute deviation from the mean.

This result means that if, for example, a GPS receiver delivers a rover position $\tilde{p}$ with uncertainty covariance $R$ then it is correct to describe the knowledge of our position as a normal distribution $\mathcal{N}(\tilde{p}, R)$.
Note that the measurement process starts with the *true* position $p$ and adds the noise to it.
I believe this result is more subtle than people realize and not trivially <<intuitive>>.
It can be correctly understood only in a Bayesian paradigm as a posterior conditional distribution.

In practice, we skip this derivation and directly initialize state mean and covariance (or their parts) from such measurements (keeping in mind its exact meaning).

# Detailed technical solution

Now let's formally derive how we should update conditional mean and covariance from epoch $k$ to $k + 1$.
First we formally introduce an augmented state vector including the control signal: 
$$
x^a_k = \begin{bmatrix} x_k \\\\ u_k \end{bmatrix}
$$

At this point we can assume that $u_k$ is completely uncertain. 
As described in the previous section we initialize it by processing the measurement $\tilde{u}_k$ to get the following conditional mean and covariance:
$$
\hat{x}^a_k = \begin{bmatrix}
\hat{x}_k \\\\
\tilde{u}_k
\end{bmatrix} \\\\
P^a_k = \begin{bmatrix}
P_k & 0 \\\\
0 & Q_k
\end{bmatrix}
$$

The transition equation can be written for the augmented vector in the form:
$$
x^a_{k + 1} = F^a_k x^a_k, \text{where } 
F^a_k = \begin{bmatrix}
F_k & B_k \\\\
0 & 0
\end{bmatrix}
$$
The prediction of the augmented condition mean and covariance has the simple form:
$$
\hat{x}^a_{k + 1} = F^a_k \hat{x}^a_k \\\\
P^a_{k + 1} = F^a_k P^a_k (F^a_k)^T
$$
Substituting all the augmented forms and keeping the result only for $\hat{x}$ and $P$ we get:
$$
\hat{x}_{k + 1} = F_k \hat{x}_k + B_k \tilde{u}_k \\\\
P\_{k + 1} = F_k P_k F_k^T + B_k Q_k B_k^T
$$
We've arrived to the expected result.

For nonlinear estimation algorithms we use this result as a building block and apply usual linearization techniques etc. with the same justification.

# Conclusion

Here I've presented a somewhat tedious and formal deviation which showed that the noisy control inputs often arising in estimation problems can be treated in a usual intuitive manner while staying in the rigorous Bayesian framework.

# References

1. P. S. Maybeck <<Stochastic Models Estimation and Control, Vol. 1>>