---
title: 'Linear Gaussian mixture filter'
date: 2025-09-09
katex: true
---

Gaussian mixture state representation seems like a powerful approach, which can work as a useful generalization of unimodal Gaussian representation.
However, its coverage in literature is scarce and unsatisfactory, not to say confusing.
It makes sense to start building understanding from the a, but fundamental problem --- estimation in linear systems.
Here a linear filter which estimates parameters of Gaussian mixture probability density is rigorously developed.
Conceptually it will be very similar to classical Kalman filter.

# Notation and definitions

First define a multivariate normal probability density function with parameters of mean $\mu$ and covariance matrix $P$:
$$
\mathcal{N}(x; m, P) = \frac{1}{(2 \pi)^{n/2} |P|^{1/2}} \exp \left[ -\frac{1}{2} (x - m)^T P^{-1} (x - m) \right]
$$
The Gaussian mixture is defined as a weighted sum of such functions:
$$
\mathcal{N}(x; w, m, P) = \sum_{i = 1}^N w^{(i)} N(x; m^{(i)}, P^{(i)}), \\\\
\text{with } w^{(i)} > 0 \text{ and } \sum_{i = 1}^N w^{(i)} = 1
$$
Here $w, m, P$ denote a set of $N$ parameters.

A sampling from such distribution is done in the following way:

1. Randomly select from which component to generate a sample with probabilities equal to weights $w^{(i)}$
2. Generate a sample using a normal distribution of the selected component

The Gaussian mixture should not be confused in any way with sum of normal variables.

# Characteristic functions

A characteristic function of a random variable $x$ with probability density $f_x(\xi)$ is defined as the following expectation ($j$ is a imaginary unit):
$$
\overline{f}_x(\alpha) = \mathrm{\mathbb{E}}_x e^{j \alpha^T x} = \int f_x(\xi) e^{j \alpha^T \xi} d \xi
$$

Because this transform is similar to Fourier, there exists a one-to-one mapping between probability density and characteristic functions.
Let's write out characteristic functions of normal and Gaussian mixture distributions:

1. For normal distribution:
   $$
   \overline{\mathcal{N}}(\alpha; m, P) = \exp \left(j \alpha^T m - \frac{1}{2} \alpha^T P \alpha \right)
   $$
   This can be shown directly from the definition by computing the integral.
2. For Gaussian mixture distribution:
   $$
   \overline{\mathcal{M}}(\alpha; w, m, P) = \sum_{i = 1}^N w^{(i)} \exp \left(j \alpha^T m^{(i)} - \frac{1}{2} \alpha^T P^{(i)} \alpha \right)
   $$
   It follows from linearity of integration and the form of the density as a linear combination of individual normal densities.

We will also use two basic properties of a characteristic function:

1. Characteristic function of sum of two independent random variables $z = x + y$:
   $$
   \overline{f}_z(\alpha) = \overline{f}_x(\alpha) \overline{f}_y(\alpha)
   $$
   It follows from the factorization property of exponent and ability to split the integral into independent parts.
2. Characteristic function of linearly transformed random variable $y = A x$:
   $$
   \overline{f}_y(\alpha) = \overline{f}_x(A^T \alpha)
   $$
   It follows directly from definition.

# Operations on Gaussian mixture variables

First consider how a Gaussian mixture variable transforms when adding an independent normal variable to it:
$$
x \sim \mathcal{M}(w_x, m_x, P_x) \\\\
y \sim \mathcal{N}(m_y, P_y) \\\\
z = x + y
$$
By using the sum property of characteristic functions we get:
$$
\begin{split}
\overline{f}_z(\alpha) = \overline{\mathcal{M}}(\alpha; w_x, m_x, P_x) \overline{\mathcal{N}}(\alpha; m_y, P_y) = \sum\_{i = 1}^N w^{(i)} \exp \left(j \alpha^T m_x^{(i)} - \frac{1}{2} \alpha^T P_x^{(i)} \alpha \right) \exp \left(j \alpha^T m_y - \frac{1}{2} \alpha^T P_y \alpha \right) = \\\\ = \sum\_{i = 1}^N w^{(i)} \exp \left(j \alpha^T (m_x^{(i)} + m_y) - \frac{1}{2} \alpha^T (P_x^{(i)} + P_y) \alpha \right) = \overline{\mathcal{M}}(\alpha; w_x, m_x + m_y, P_x + P_y)
\end{split}
$$
We can identify that $z$ has characteristic function of a Gaussian mixture and thus (the notation assumes that each element of the parameter set is modified):
$$
z \sim M(w_x, m_x + m_y, P_x + P_y)
$$
In other words, the mixture weights don't change and parameters of individual components are changed by the usual rules for mean and variance of a sum.

Now consider how a Gaussian mixture variable changes under a linear transform:
$$
x \sim M(w_x, m_x, P_x) \\\\
y = A x
$$
By using the linear transform rule of characteristic function we get (with the same simplified notation):
$$
\overline{f}_y(\alpha) = \overline{\mathcal{M}}(A^T \alpha; w_x, m_x, P_x) = \overline{\mathcal{M}}(\alpha; w_x, A m_x, A P_x A^T)
$$
We see, that $y$ is a Gaussian mixture variable with the following parameters:
$$
y \sim M(w_x, A m_x, A P_x A^T)
$$
In other words, the mixture weights don't change and parameters of individual components are transformed in the usual way.

# Recap of the Bayesian derivation of Kalman update step

Let's recap how Kalman update formulas are derived using the Bayes rule.
Starting from state and measurement models:
$$
x \sim \mathcal{N}(m, P) \\\\
z = H x + v, v \sim \mathcal{N}(0, R)
$$
The Kalman update step infers the probability density of $x$ conditioned on measurement $z$ using the Bayes rule:
$$
f\_{x | z}(\xi | \zeta) = \frac{f\_{z | x} (\zeta | \xi) f_x(\xi)}{f_z(\zeta)}
$$
Let's write individual probability densities:
$$
f_x(\xi) = \mathcal{N}(\xi; m, P) \\\\
f_{z | x}(\zeta | \xi) = \mathcal{N}(\zeta; H \xi, R) \\\\
f_z(\zeta) = \mathcal{N}(\zeta; H m, H P H^T + R)
$$
The key result is that after substitution it can be reduced to another normal density of the form:
$$
f_{x | z}(\xi | \zeta) = \mathcal{N}(\xi; m^+, P^+)
$$
With the updated mean and covariance given by the famous Kalman formulas:
$$
m^+ = m + K (\zeta - H m) \\\\
P^+ = (I - K H) P \\\\
K = P H^T (H P H^T + R)^{-1}
$$
Here $\zeta$ denotes a particular observed realization of $z$.

Later we will use this result in the following form:
$$
\mathcal{N}(\zeta; H \xi, R) \mathcal{N}(\xi; m, P) = \mathcal{N}(\zeta; H m, H P H^T + R) \mathcal{N}(\xi; m^+, P^+)
$$

# Derivation of linear Gaussian mixture filter

Now we have everything setup to derive a linear Gaussian mixture filter.
The filter will estimate probability density function of a vector $x$ given the following transition and measurement equations:
$$
x_{k + 1} = F_k x_k + w_k, w_k \sim \mathcal{N}(0, Q_k) \\\\
z_k = H_k x_k + v_k, v_k \sim \mathcal{N}(0, R_k)
$$
But as opposed to Kalman filter the initial distribution is a Gaussian mixture:
$$
x_0 \sim \mathcal{M}(w_0, m_0, P_0)
$$
As will be seen by an induction argument, the conditional distribution of $x$ will remain a Gaussian mixture. 
And thus the computational state of the filter is the set of weights, mean vectors and covariance matrices --- $w_k, x_k, P_k$ (here $k$ denotes the time index).

Now let's derive formulas for prediction and update phases of the filter.

## Prediction step

The prediction step propagates distribution parameters from epoch $k$ to $k + 1$:
$$
x_k \sim \mathcal{M}(w_k, m_k, P_k) \\\\
x_{k + 1} = F_k x_k + w_k
$$
Considering that $w_k$ is independent of $x_k$ and applying derived formulas for Gaussian mixture variable transformations we get:
$$
w_{k + 1} = w_k \\\\
m_{k + 1}^{(i)} = F_k m_k^{(i)} \\\\
P_{k + 1}^{(i)} = F_k P_k^{(i)} F_k^T + Q_k
$$
That is the weights of components don't change and parameters of individual components are transformed as in Kalman filter.

## Update step

The update step incorporates measurement information by computing the conditional probability density according to the measurement model:
$$
x_k \sim \mathcal{M}(w_k, m_k, P_k) \\\\
z_k = H_k x_k + v_k, v_k \sim \mathcal{N}(0, R_k)
$$
Analogous to Kalman filter we apply the Bayes rule:
$$
f\_{x | z}(\xi | \zeta) = \frac{f\_{z | x} (\zeta | \xi) f_x(\xi)}{f_z(\zeta)}
$$
Where the individual probability densities are (time index $k$ is omitted):
$$
f_x(\xi) = \mathcal{M}(\xi; w, m, P) \\\\
f_{z | x}(\zeta | \xi) = \mathcal{N}(\zeta; H \xi, R) \\\\
f_z(\zeta) = \mathcal{M}(\zeta; w, H m, H P H^T + R)
$$
By substituting the density functions into the Bayes rule we get:
$$
f_{x | z}(\xi | \zeta) = \frac{\sum_{i = 1}^N w^{(i)} \mathcal{N}(\zeta; H \xi, R) \mathcal{N}(\xi; m^{(i)}, P^{(i)})}{\sum_{i = 1}^N w^{(i)} \mathcal{N}(\zeta; H m^{(i)}, H P^{(i)} H^T + R)}
$$
Now in the numerator use the probability density relation from the Kalman filter update step:
$$
f_{x | z}(\xi | \zeta) = \frac{\sum_{i = 1}^N w^{(i)} \mathcal{N}\left(\zeta; H m^{(i)}, H P^{(i)} H^T + R \right) \mathcal{N}\left(\xi; m^{(i)+}, P^{(i)+}\right) }{\sum_{i = 1}^N w^{(i)} \mathcal{N}(\zeta; H m^{(i)}, H P^{(i)} H^T + R)}
$$
Now the variable $\xi$ appears only in one multiplier and the whole probability density can be identified as a Gaussian mixture with the parameters:
$$
w^{(i)+} = \frac{w^{(i)} p^{(i)}}{\sum_{i = 1}^N w^{(i)} p^{(i)}} \\\\ [3pt]
m^{(i)+} = m^{(i)} + K^{(i)} (\zeta - H m^{(i)}) \\\\
P^{(i)+} = (I - K^{(i)} H) P^{(i)}
$$
with
$$
K^{(i)} = P^{(i)} H^T (H P^{(i)} H^T + R)^{-1} \\\\
p^{(i)} = \mathcal{N}\left(\zeta; H m^{(i)}, H P^{(i)} H^T + R \right)
$$
The weights adjustment factors $p^{(i)}$ are essentially likelihoods to observe a measurement $\zeta$ assuming $x$ comes from a particular component.

## Filter summary

The filter operates in a standard sequence:

1. The state parameters $w_0, x_0, P_0$ are initialized from a prior knowledge 
2. Then the parameters are updated sequentially for each time epoch $k$:
    1. Measurements at epoch $k$ are processed as described in the corresponding section to obtain the updated parameters $w_k^+, m_k^+, P_k^+$
    2. The parameters are propagated from epoch $k$ to $k + 1$ as described in the corresponding section to obtain the new parameters $w_{k +1}, m_{k + 1}, P_{k + 1}$

The filter seems to be well and rigorously defined and will correctly propagate the Gaussian mixture distribution according to the given linear model.

# Conclusion

It was shown that for a linear model we can rigorously build an estimation filter for a Gaussian mixture distribution, which is very similar to Kalman filter.
It doesn't seem to have any problematic or questionable aspects.
However, its generalization to nonlinear systems looks to be more difficult and subtle.
This is something I want to figure out next.
