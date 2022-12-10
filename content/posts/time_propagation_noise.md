---
title: 'Noisy input in state time propagation: "control" vs. "navigation"'
date: 2022-12-10
katex: true
---

In this note I want to show that there are two possible scenarios how noisy input can influence state time propagation in estimation problems.
This is not well articulated in classical literature and may create confusion when applying estimation algorithms to practical systems.

I will illustrate the concept on a simple example.
Consider a robot which can move in the direction of its longitudinal axis and rotate.
Its state is described by 4 variables (2D case):

* $x, y$ -- euclidean coordinates in a selected global frame of reference
* $v$ -- velocity along the longitudinal axis
* $\psi$ -- heading angle relative to the global frame of reference

The state variables follow the differential equations:
$$
\dot{x} = v \cos \psi \\\\
\dot{y} = v \sin \psi \\\\
\dot{v} = a \\\\
\dot{\psi} = \omega
$$
Where $a$ and $\omega$ are acceleration and angular rate respectively.

Now consider two possible scenarios.

# Control scenario

In this scenario required control signals $a(t)$ and $\omega(t)$ are determined before the mission, but they are executed imperfectly by the robot as 
$$
\tilde{a} = a + w_a \\\\
\tilde{\omega} = \omega + w_\omega
$$
where $w_a$ and $w_\omega$ are random disturbances of acceleration and angular rate respectively.

The true state variables are stochastic and obey the equations with disturbed acceleration and angular rate:
$$
\dot{x} = v \cos \psi \\\\
\dot{y} = v \sin \psi \\\\
\dot{v} = \tilde{a}\\\\
\dot{\psi} = \tilde{\omega}
$$

An estimation algorithm (like Extended Kalman Filter) executes the time propagation equations with nominal acceleration ang angular rate:
$$
\dot{\hat{x}} = \hat{v} \cos \hat{\psi} \\\\
\dot{\hat{y}} = \hat{v} \sin \hat{\psi} \\\\
\dot{\hat{v}} = a \\\\
\dot{\hat{\psi}} = \omega
$$
where $\hat{\cdot}$ denotes estimate of the corresponding variable.

Here we have that the actual motion of the robot is stochastic and the estimation equations are deterministic.
This is how estimation problems and algorithms are usually formulated in literature.

# Navigation scenario

In this scenario the control signals to the robot are not set (it can be controlled by an operator, for example), but acceleration and angular rate are measured by sensors installed on the robot as
$$
\tilde{a} = a + w_a \\\\
\tilde{\omega} = \omega + w_\omega
$$
where $w_a$ and $w_\omega$ are random measurement noise of acceleration and angular rate respectively.

The true state variables are unknown but deterministic and obey the equations with true acceleration and angular rate:
$$
\dot{x} = v \cos \psi \\\\
\dot{y} = v \sin \psi \\\\
\dot{v} = a\\\\
\dot{\psi} = \omega
$$

An estimation algorithm executes the time propagation equations with measured acceleration and angular rate:
$$
\dot{\hat{x}} = \hat{v} \cos \hat{\psi} \\\\
\dot{\hat{y}} = \hat{v} \sin \hat{\psi} \\\\
\dot{\hat{v}} = \tilde{a} \\\\
\dot{\hat{\psi}} = \tilde{\omega}
$$

Here we have that the actual motion of the robot is deterministic and the estimation equations are stochastic.
This situation occurs when the noisy measurements are used in the estimate time propagation equations. 
An important example is an Inertial Navigation System (INS), where noisy measurements of gyros and accelerometers are used to propagate the system state estimate (i. e. position, velocity and attitude).
Our example in this section can be though of as a simplified INS.

# General treatment and consequences on estimation algorithms design

Let's write the state and estimate propagation equations in a general form for a state vector $x$ with a noise vector $w$.

For the control scenario:
$$
\dot{x} = f(x, w) \\\\
\dot{\hat{x}} = f(\hat{x}, 0)
$$
For the navigation scenario:
$$
\dot{x} = f(x, 0) \\\\
\dot{\hat{x}} = f(\hat{x}, w)
$$

Consider the error of the estimate:
$$
\Delta x = \hat{x} - x
$$
In both scenarios an estimation algorithm is sought to be unbiased with the minimal error variance:
$$
\operatorname{E} \Delta x = 0 \text{ and } \operatorname{E} (\Delta x)^T \Delta x \text{ is minimal}
$$
A linearized differential equation for the error can be obtained by subtracting the state propagation equations from the estimate propagation equations and linearizing:
$$
\Delta \dot{x} = F(x) \Delta x \pm G(x) w \\\\
\text{with } F(x) = \frac{\partial f(x)}{\partial x} \text{ and } G(x) = \frac{\partial f(x)}{\partial w}
$$
Plus and minus are used for navigation and control scenarios respectively. 
Considering that $w$ is usually modelled as a zero mean noise, this difference is not substantial.

With the assumption that $w$ is a zero mean white noise, the equation for the error covariance time propagation is the same in both scenarios:
$$
\dot{P} = F P F^T + G Q G^T
$$
where $Q$ is power spectral density matrix of the noise $w$.

We can see that there is no difference in estimation algorithms design besides whether noisy input is used in the estimate time propagation (<<navigation scenario>>) or not (<<control scenario>>).
Mixed scenarios are possible too -- with a conceptual understanding of the subject it's not difficult to properly design and reason about estimation algorithms in such cases too.

# Conclusion

I've highlighted that in estimation problems a noisy input can influence the true state time propagation or the state estimate time propagation.
It doesn't change an estimation algorithm's objective or implementation logic, but it is useful to keep that in mind for clear intuitive engineering understanding of estimation theory and algorithms.
Without understanding this point one can get a bit stuck when trying to design an estimation algorithm for a system like INS, where noisy measurements are used for the state estimate propagation.
