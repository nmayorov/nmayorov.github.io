---
title: "Heading observability in Inertial Navigation Systems"
date: 2023-04-02
katex: true
---

It's generally known that the heading angle in Inertial Navigation Systems (INS) is observable when the system undergoes acceleration (including turns at nonzero velocity).
However if we look at the ["modified error model"]({{<ref "/content/posts/ins_error_equations.md#modified-phi-vector-error-model">}}), 
it may seem that merely significantly nonzero constant velocity is sufficient.
Here I want to figure out this subject using formal observability criteria and then give it intuitive physical explanation.

# INS equations in 2 dimensions

To investigate heading observability, a simplified version of INS for a 2-dimensional motion is fully sufficient.
The state is described relative to the world frame $n$ by 5 variables:

- Position vector $r^n$ with 2 components
- Velocity vector $v^n$ with 2 components
- Heading angle $\psi$

The state variables obey the following differential equations:
$$
\dot{r}^n = v^n \\\\
\dot{v}^n = C^n_b f^b \\\\
\dot{\psi} = \omega
$$
Here $f^b$ and $\omega$ are acceleration and (vertical) angular rate measured by IMU.
And $C^n_b = R(\psi)$ is a 2-dimensional rotation matrix computed from the heading angle:
$$
R(\alpha) \coloneqq \begin{bmatrix}
\cos \alpha & -\sin \alpha \\\\
\sin \alpha & \cos \alpha
\end{bmatrix}
$$

The estimates state variables we denote as $\hat{r}, \hat{V}, \hat{\psi}$ 

# Error equations

To use a linear theory we need to derived a linear system for the INS errors.

## Natural error parametrization

Natural error parameters:

- Position error: $\Delta r = \hat{r} - r^n$
- Velocity error: $\Delta V = \hat{V} - V^n$
- Heading error: $\Delta \psi = \hat{\psi} - \psi$

The error equations:
$$
\Delta \dot{r} = \Delta v \\\\
\Delta \dot{v} = (\mathcal{S} C^n_b f^b) \Delta \psi + C^n_b \Delta f_b \\\\
\Delta \dot{\psi} = \Delta \omega
$$
With

- $\Delta f^b$ and $\Delta \omega$ -- accelerometer and angular rate measurement errors
- $S \coloneqq \begin{bmatrix} 0 & -1 \\\\ 1 & 0 \end{bmatrix}$

## Modified error parametrization

Here the velocity error is measured relative to the "platform frame" (where accelerometer readings are projected):
$$
\delta v = \hat{v} - R(\Delta \psi) v^n \approx \Delta v - (\mathcal{S} v^n) \Delta \psi
$$

The error equations:
$$
\Delta \dot{r} = \delta v + (\mathcal{S} v^n) \Delta \psi \\\\
\delta \dot{v} = -(S v^n) \Delta \omega + C^n_b \Delta f^b \\\\
\Delta \dot{\psi} = \Delta \omega
$$

# Observability analysis

Let's assume that velocity measurements in the world frame are available (like from GNSS).
From the standard error model we see that the heading error is coupled into the velocity error when there is acceleration.
Thus we may assume that the heading is observable only when the vehicle undergoes acceleration.

On the other hand, if we consider the modified error model, the heading error is coupled into the true velocity error as 
$$
\Delta v = \delta v + (\mathcal{S} v^n) \Delta \psi
$$
And in the similar fashion into the position error growth.
Thus it looks like that the heading may be observable from the velocity measurements when the velocity is nonzero.

These two conclusions contradict each other and therefore we want to carry out a formal observability analysis.
To do that for a time-invariant linear system with $n$ variables with the system matrix $F$ and the observation matrix $H$ one should check the rank of the matrix:
$$
M = \begin{bmatrix}
H \\\\
H F \\\\
\cdots \\\\
H F^{n-1}
\end{bmatrix}
$$
For a fully observable system this matrix must have $n$ linearly independent rows.

To use this result we must consider a time-invariant system for the INS errors with constant heading and velocity/acceleration.
For simplicity we assume $\psi = 0$ and error-free inertial measurements.

## Standard error model with constant acceleration

We can exclude the position error variables as irrelevant. 
Then write the error model componentwise:
$$
\Delta \dot{v}_1 = -f_2 \Delta \psi \\\\
\Delta \dot{v}_2 =  f_1 \Delta \psi \\\\
\Delta \dot{\psi} = 0
$$
The velocity error components are observed:
$$
z_1 = \Delta v_1 \\\\
z_2 = \Delta v_2
$$
Thus we have the following matrices:
$$
F = \begin{bmatrix}
0 & 0 & -f_2 \\\\
0 & 0 & f_1  \\\\
0 & 0 & 0
\end{bmatrix} \\\\
H = \begin{bmatrix}
1 & 0 & 0 \\\\
0 & 1 & 0
\end{bmatrix}
$$
The observability test matrix (the first 4 rows):
$$
M = \begin{bmatrix}
1 & 0 & 0 \\\\
0 & 1 & 0 \\\\
0 & 0 & -f_2 \\\\
0 & 0 & f_1
\end{bmatrix}
$$
We can clearly see that when $f \neq 0$ the system is fully observable. 
That is the velocity components and the heading angle will be determined exactly if there is no noise in the system.

## Modified error model with constant velocity

The error equations:
$$
\delta \dot{v}_1 = 0 \\\\
\delta \dot{v}_2 = 0 \\\\
\Delta \dot{\psi} = 0
$$
The full velocity error are observed:
$$
z_1 = \delta v_1 - v_2 \Delta \psi \\\\
z_2 = \delta v_2 + v_1 \Delta \psi 
$$
The system and observations matrices:
$$
F = \begin{bmatrix}
0 & 0 & 0 \\\\
0 & 0 & 0 \\\\
0 & 0 & 0
\end{bmatrix} \\\\
H = \begin{bmatrix}
1 & 0 & -v_2 \\\\
0 & 1 & v_1
\end{bmatrix}
$$
With $F$ being zero it is obvious that the system is not fully observable.
Only the two linear combinations as expressed by $z_1$ and $z_2$ are observable.

If we add measurements of the velocity in the body frame (by radar or odometer):
$$
z_3 = \delta v_1 \\\\
z_4 = \delta v_2
$$
We get the observation matrix
$$
H = \begin{bmatrix}
1 & 0 & -v_2 \\\\
0 & 1 & v_1 \\\\
1 & 0 & 0 \\\\
0 & 1 & 0
\end{bmatrix}
$$
It has 3 independent rows as long $v \neq 0$. 
That is the system is fully observable when there is GNSS and odometer measurements are available and the speed is not zero.

# Discussion

For GNSS aides INS, the heading is observable only when the system undergoes accelerations.
With the constant nonzero velocity the heading error is coupled into position error growth (or "true" velocity errors), however there is no way to separate the heading error and the "odometric" velocity errors.
If body velocity measurements are added, then the system becomes fully observable.

These results can be understood using the following principle: heading angle is observed by matching the same nonzero vector measured in the body and the world frames.

The acceleration vector is measured in the body frame by IMU and in the world frame it can be inferred by differentiation of the measured velocity.
Only in a navigation filter it is done implicitly, but the core principle is as stated.
If there is no acceleration, there is no vector to match and heading estimation is impossible.

If there are velocity measurements in the body frame, they can be matched with velocity measurements in the world frame and the heading is estimated from that.
Again, in real systems it is done implicitly by means of a Kalman filter of some sort.
In practice it means that velocity measurements in the body frame (along with GNSS velocity) are useful to enhance heading estimation.
