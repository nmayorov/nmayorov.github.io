---
title: "Strapdown IMU synthesis"
date: 2024-03-03
---

Synthesis of strapdown Inertial Measurement Unit (IMU) readings is required for simulation of inertial navigation systems and related algorithms development.
Here I want to describe the algorithm of IMU synthesis implemented in [pyins](github.com/nmayorov/pyins), which is probably one of the most useful functions in the library.

# Problem formulation

The task is to compute strapdown IMU readings from a time series of position and attitude. 
The following points define the algorithm requirements:

- Data is provided at known time points $t_k$, not necessary equispaced
- Position is given as a time series of latitude, longitude and altitude -- $\varphi(t_k), \lambda(t_k), h(t_k)$
- Attitude is given as a time series of roll, pitch and heading angles -- $\gamma(t_k), \theta(t_k), \psi(t_k)$
- The algorithm must account for Earth rotation and ellipticity and use realistic gravity model
- The algorithm must be able to compute rate and increment (integral) readings

# Preliminarily considerations

We want to generate IMU from sampled trajectory points instead of some predefined continuous functions to have more flexibility and convenience.
Figuring out and implementing appropriate continuous functions for position and attitude is difficult task on its own, which is better to be avoided.

With sampled trajectory points it is obvious that formally the derivatives and thus IMU readings can be arbitrary, i.e. the function might exhibit erratic behavior still passing through sampled points.
The key and necessary assumption here is that the underlying function is *bandlimited* and sampled with sufficient frequency (formally more frequent than Nyquist frequency if time points are equispaced).
Frequency domain approaches might in fact be the most appropriate and insightful, but the fact that rotations belong to a non-euclidean space (SO(3) manifold) complicates things significantly.

Instead, the proposed algorithm relies on cubic spline interpolation, which provides twice differentiable functions.
In my understanding it works as some kind of approximation of Shannon's optimal interpolation (sinc interpolation).

Another known approach to generate IMU is to invert a high accuracy strapdown integration algorithm.
I find this approach conceptually not satisfactory as quite complicated, coupled with the integration scheme and not clear in terms of properties of the resulting IMU readings.

Now leaving aside somewhat vague conceptual speculations let's move to the actual algorithm.

# Algorithm derivation

Coordinate frames, navigation definitions and equations are explained [here]({{<ref "/posts/ins_error_equations">}}).

## Definition of IMU readings

First let's define precisely what IMU measures.

Accelerometers measure "specific force", i.e. total force excluding gravity, computed for a unit mass as acceleration:
$$
f^b = C^b_i (\ddot{r}^i - g_i^i)
$$
Where $r^i$ is the body position relative to ECI frame, twice differentiated to get the true inertial acceleration.
And $g_i^i$ is acceleration due to Earth mass attraction, not including the centrifugal force.

Gyroscopes measure angular rate around its axes relative to the inertial frame -- $\omega_{ib}^b$.

The rate IMU readings are then $f^b(t_k)$ and $\omega_{ib}^b(t_k)$.
The increment IMU readings are defined as integrals:
$$
\upsilon(t_k) = \int_{t_{k - 1}}^{t_k} f^b(t) d t \\\\
\alpha(t_k) = \int_{t_{k - 1}}^{t_k} \omega_{ib}^b(t) d t
$$

## Building interpolation splines

First we compute "inertial longitude" to account for Earth rotation:
$$
\lambda_i(t_k) = \lambda(t_k) + \omega_e t_k
$$
And then compute position and attitude relative to ECI:
$$
r^i(t_k) = r^i(\varphi(t_k), \lambda_i(t_k), h(t_k)) \\\\
C^i_b(t_k) = C^i_n(\varphi(t_k), \lambda_i(t_k)) C^n_b(\gamma(t_k), \theta(t_k), \psi(t_k))
$$
All the involved variables are computed as known functions of the inputs.

For $r^i$ and $C^i_n$ we then build cubic interpolation splines, which in `scipy` are implemented as `CubicSpline` and `RotationSpline` respectively.
Specific boundary conditions for `CubicSpline` don't seem to be very important, default conditions are used in `pyins`.
Note that `RotationSpline` interpolates rotations with continuous angular rate and acceleration and does a lot of heavy lifting inside.

Let's formally denote these splines as $S_p$ and $S_a$ respectively.

## Computing rate readings

Rate readings are trivially computed from the splines:
$$
f^b(t_k) = C^b_i(t_k) (\ddot{S}\_p(t_k) - g_i^i(\varphi(t_k), \lambda_i(t_k), h(t_k))) \\\\
\omega_{ib}^b(t_k) = \dot{S}_a(t_k)
$$
Here $g_i^i$ is computed as known function of latitude, inertial longitude and altitude and the derivative of $S_a$ is assumed to give the angular rate.

## Computing incremental readings

Computation of incremental readings is more difficult because of two facts:

- Coefficients of `RotationSpline` define a rotation vector, but the angular velocity is a nonlinear function of the vector and its derivative
- Computation of the specific force in the body frame couples rotation and position splines, moreover the rotation transform is a nonlinear function of the rotation vector

Both of these difficulties can be resolved analytically with the assumption that the rotation vector on each interval is small.
Note that if this is not true, then IMU readings are almost useless for navigation anyway.

We use a second-order (in the norm of the rotation vector $\theta_k$) approximations:
$$
C^i_{bk}(\tau) \approx C^i_b(t_{k - 1}) \left(I + (\theta_k(\tau) \times) + \frac{1}{2} (\theta_k(\tau) \times)^2 \right) \\\\
\omega_{ibk}^b(\tau) \approx \left(I - \dfrac{1}{2} (\theta_k(\tau) \times) + \frac{1}{6} (\theta_k(\tau) \times)^2 \right) \dot{\theta}_k(\tau) \\\\
\text{ with } \tau = t - t\_{k - 1} \text{ and } t\_{k - 1} \leq t < t_k \\\\ 
$$

### Gyro readings

The rotation vector at $k$-th interval is expressed as a cubic function of the time offset (with coefficients available in `RotationSpline`):
$$
\theta_k(\tau) = a_k \tau + b_k \tau^2 + c_k \tau^3
$$
Substituting this into the approximate formula for the angular rate we get a 7-th order power series in $\tau$.
The coefficients in this power series contain the spline coefficients, their single and double cross-products.
The cross-products can be attributed to so-called rotation "coning" effects.
It then integrated from 0 to $t_k - t_{k - 1}$ to get $\alpha(t_k)$.

### Accelerometer readings

A linear model for the specific force in the inertial frame is built by twice differentiating $S_p$ and using piecewise-linear interpolation for $g_i^i$ (it changes slowly with coordinates):
$$
f_k^i(\tau) = d^i_k + e^i_k \tau
$$
Specific force in the body frame is then:
$$
\begin{split}
f_k^b(\tau) = C^b_{ik}(\tau) f_k^i(\tau) \approx \left(I - (\theta_k(\tau) \times) + \frac{1}{2} (\theta_k(\tau) \times)^2 \right) C^b_i(t_{k - 1}) (d^i_k + e^i_k \tau) = \\\\ 
= \left(I - (\theta_k(\tau) \times) + \frac{1}{2} (\theta_k(\tau) \times)^2 \right) (d_k + e_k \tau)
\end{split}
$$
Where the projected coefficients were defined:
$$
d_k \coloneqq C^b_i(t_{k - 1}) d^i_k \\\\
e_k \coloneqq C^b_i(t_{k - 1}) e^i_k
$$
Substituting $\theta_k(\tau)$ we get another 7-th order power series in $\tau$ (with different cross-products attributed to "sculling" effects).
It then integrated from 0 to $t_k - t_{k - 1}$ to get $\upsilon(t_k)$.

The coefficients of the aforementioned series can be found in [pyins.sim](https://github.com/nmayorov/pyins/blob/master/pyins/sim.py) source file.

# Practical considerations and discussion

From the experience, the implemented method does require more or less smooth input trajectory.
Typically, the required level of smoothing is determined empirically by checking if generated readings have adequate magnitude and whether they give an accurate trajectory after the integration.
The [example](https://github.com/nmayorov/pyins/blob/master/examples/ideal_imu.ipynb) illustrates this point.

How to characterize or decrease this "noise amplification", what is the required ratio between sampling and motion frequency -- it requires investigation, perhaps digging into frequency domain properties of spline interpolation.
