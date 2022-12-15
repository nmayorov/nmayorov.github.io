---
title: "INS error equations"
date: 2022-12-15
katex: true
---

In this note I want to derive Inertial Navigation System (INS) error equations from the first principles.
This is of course a well known subject, however it's useful to derive all the equations from scratch and revisit the subject from time to time.

# INS state variables and equations for them

An Inertial Navigation System computes position, velocity and attitude of a moving object using measurements from an Inertial Measurement Unit (IMU).

To formulate necessary equations frames of reference need to be defined (denoted by lowercase letters):

* $e$ -- ECEF frame 
* $n$ -- NED frame
* $b$ -- body frame with axes along inertial sensor axes

A vector $u$ expressed in a frame $a$ is denoted as $u^a$.

Variables computed by INS are:

* Position as latitude, longitude and altitude -- $\varphi, \lambda, h$
* Velocity relative to the ground (frame $e$) expressed in frame $n$ -- $V^n = \begin{bmatrix} V_N & V_E & V_D \end{bmatrix}^T$
* Attitude as rotation matrix from frame $n$ to frame $b$ -- $C^n_b$

The differential equations for the state variables are:
$$
\dot{\varphi} = \frac{V_N}{R_N + h} \\\\[1 pt]
\dot{\lambda} = \frac{V_E}{(R_E + h) \cos \varphi} \\\\
\dot{h} = -V_D \\\\
\dot{V}^n = -(2 \Omega^n + \rho^n) \times V^n + g^n + C^n_b f^b \\\\
\dot{C}^n_b = C^n_b (\omega^b \times) - ([\Omega^n + \rho^n] \times) C^n_b
$$
Where the following variables were introduce:

* $R_N$ and $R_E$ -- North and East curvature radii
* $\Omega^n$ -- Earth angular rate
* $\rho^n$ -- transport rate, i. e. angular rate of $n$ frame relative to $e$ frame due to object's velocity
* $g^n$ -- gravity vector
* $f^b$ -- accelerometer measurements of specific force
* $\omega^b$ -- gyro measurements of angular rate relative to an inertial frame
* $(\cdot \times)$ -- cross-product skew symmetric matrix defined as $(u \times) v \equiv u \times v$

The INS computer updates the state variables according to the equations using gyro and accel measurements containing errors:
$$
\tilde{\omega} = \omega^b + \Delta \omega \\\\
\tilde{f} = f^b + \Delta f \\\\
$$
The computed or estimated variables will be marked with hat.

# Error definitions

Errors must express deviation of estimated variables from true variables.
They don't need to be direct differences of the variables.

There are several conventional ways to introduce error variables for INS.
I will choose the most straightforward approach from which it's possible to derive equations for other parametrizations.

The selected parametrization are usually called "phi-angle error model".

## Position error 

Direct differences of latitude, longitude and altitude are not particularly convenient to interpret as well as to derive equations for.
Instead we can use correspondence of latitude, longitude, altitude and coordinates in ECEF, which we denote as $r^e$.

Define the position error as:
$$
\Delta r^n = C^n_e (\hat{r}^e - r^e)
$$
That is it ECEF coordinate errors projected into NED frame.

Relations between $\Delta r^n$ and geographic coordinates errors are
$$
\Delta r_N = (R_N + h) \Delta \varphi \\\\
\Delta r_E = (R_E + h) \cos \varphi \Delta \lambda \\\\
\Delta r_N = - \Delta h
$$

## Velocity error

The velocity error is defined as direct difference:
$$
\Delta V^n = \hat{V}^n - V^n
$$

## Attitude error

Attitude error is defined as a rotation vector $\phi$ associated with $n$ frame according to the following equation:
$$
C^n_b = \exp(\phi \times) \hat{C}^n_b
$$
Where the matrix exponential is used.
It has a closed form formula and to the first order in $\phi$ can be expressed as
$$
\exp(\phi \times) \approx I + (\phi \times)
$$

# Helper formulas

The curvature radii are computed as 
$$
R_N = \frac{a (1 - e^2)}{(1 - e^2 \sin^2 \varphi)^{3/2}} \\\\
R_E = \frac{a}{1 - e^2 \sin^2 \varphi}
$$
Where $a$ is Earth ellipsoid semi-major axis and $e$ is its eccentricity.
Considering that $e^2 \approx 7 \cdot 10^{-3}$ the dependency of the radii on the latitude is ignored for the error analysis.

The transport rate is computed as
$$
\rho^n =  \begin{bmatrix}
\dfrac{V_E}{R_E + h} \\\\
-\dfrac{V_N}{R_N + h} \\\\
-\dfrac{V_E \tan \varphi}{R_E + h}
\end{bmatrix} = F V^n
$$
The matrix $F$ relates a linear displacement in $n$ frame to its small angular rotation.
The error of the transport rate due to position and velocity errors is
$$
\Delta \rho^n = R \Delta r^n + T \Delta V^n
$$
Where matrices $T$ and $R$ are defined as
$$
T = \begin{bmatrix}
0 & \dfrac{1}{R_E + h} & 0 \\\\
-\dfrac{1}{R_N + h} & 0 & 0 \\\\
0 & -\dfrac{\tan \varphi}{R_E + h} & 0
\end{bmatrix} \\\\[2pt]
R = \begin{bmatrix}
0 & 0 & \dfrac{V_E}{(R_E + h)^2} \\\\
0 & 0 & -\dfrac{V_N}{(R_N + h)^2} \\\\
-\dfrac{V_E \sec^2 \varphi}{(R_E + h) (R_N + h)} & 0 & -\dfrac{V_E \tan \varphi}{(R_E + h)^2}
\end{bmatrix}
$$

The Earth rotation rate in $n$ frame has the form
$$
\Omega^n = \Omega \begin{bmatrix}
\cos \varphi \\\\
0 \\\\
-\sin \varphi
\end{bmatrix}
$$
Its perturbation due to the position error can be shown to be
$$
\Delta \Omega^n = (\Omega^n \times) T \Delta r^n
$$

The gravity model is
$$
g^n = \begin{bmatrix}
0 \\\\ 
0 \\\\
g(\varphi) \left(1 - \dfrac{2 h}{a} \right)
\end{bmatrix}
$$
The weak dependency in $g(\varphi)$ can be ignored in the error analysis.
The gravity vector perturbation due to the position error:
$$
\Delta g^n = G \Delta r^n
$$
with
$$
G = \begin{bmatrix}
0 & 0 & 0 \\\\
0 & 0 & 0 \\\\
0 & 0 & \dfrac{2 g}{a}  \\\\
\end{bmatrix}
$$


# Position error equation derivation

Differentiating the position error we get
$$
\begin{split}
\Delta \dot{r}^n = \dot{C}^n_e (\hat{r}^e - r^e) + C^n_e (\dot{\hat{r}}^e - \dot{r}^e) = -(\rho^n \times) C^n_e (\hat{r}^e - r^e) + C^n_e (\hat{V}^e - V^e) = \\\\
= -\rho^n \times \Delta r^n + C^n_e \hat{C}^e_n \hat{V}^n - V^n \approx -\rho^n \times \Delta r^n + (I + (\theta \times)) (V^n + \Delta V^n) - V^n \approx  \\\\
\approx -\rho^n \times \Delta r^n + \Delta V^n + \theta \times V^n
\end{split}
$$
Where we've used equations
$$
\dot{C}^e_n = C^e_n (\rho^n \times) \\\\
\dot{r}^e = V^e \\\\
\hat{C}^e_n \approx C^e_n (I + (\theta \times))
$$
Substituting $\theta = T \Delta r^n$ (by the same logic that the rotational perturbation $\theta$ of frame $n$ is caused by the linear perturbation $\Delta r^n$) we get:
$$
\Delta \dot{r}^n = -\left[ (\rho^n \times) + (V^n \times) T \right] \Delta r^n + \Delta V^n
$$

# Velocity error equation derivation

To derive the velocity error equation we subtract the equation for the true velocity from the equation for the INS estimated velocity keeping linear terms:
$$
\Delta \dot{V}^n = -(2 \Omega^n + \rho^n) \times \Delta V^n - (2 \Delta \Omega^n + \Delta \rho^n) \times V^n + \Delta g^n + (I - \phi \times) C^n_b (f^b + \Delta f) - C^n_b f^b
$$
Substituting expressions for all $\Delta$-s we get
$$
\Delta \dot{V}^n = -(2 \Omega^n + \rho^n) \times \Delta V^n + (V^n \times) (2 (\Omega^n \times) T \Delta r^n + R \Delta r^n + T \Delta V^n) + G \Delta r^n + (C^n_b f^b) \times \phi + C^n_b \Delta f
$$
And after the grouping the terms:
$$
\Delta \dot{V}^n = \left[2 (V^n \times) (\Omega^n \times) T + (V^n \times) R + G \right] \Delta r^n + \left[ -([2 \Omega^n + \rho^n]\times) + (V^n \times) T \right] \Delta V^n + (C^n_b f^b) \times \phi + C^n_b \Delta f
$$

# Attitude error equation derivation

Let's first derive a general equation.

## Differential equation for a small rotation vector as attitude error

Consider the rotation kinematics equation for a matrix $C$:
$$
\dot{C} = C (\omega \times) - (\xi \times) C
$$
Where $\omega$ and $\xi$ are angular rates of the frames between which $C$ transforms relative to some third frame.
And the same dynamics for its estimate:
$$
\dot{\hat{C}} = C (\hat{\omega} \times) - (\hat{\xi} \times) \hat{C} \text{ with } \hat{\omega} = \omega + \Delta \omega, \hat{\xi} = \xi + \Delta \xi
$$
Define the relation between the matrices to the first order in a rotation vector $\phi$ as 
$$
C \approx (I + \phi \times) \hat{C}
$$
From this we get the differential equation for $\phi$:
$$
\begin{split}
(\dot{\phi} \times) = \dot{C} \hat{C}^T + C \dot{\hat{C}}^T = [C (\omega \times) - (\xi \times) C ] \hat{C}^T + C[-(\hat{\omega} \times) \hat{C}^T + \hat{C}^T (\hat{\xi} \times)] = \\\\
= C ([\omega - \hat{\omega}] \times) \hat{C}^T - (\xi \times) C \hat{C}^T  + C \hat{C}^T (\hat{\xi} \times)
\end{split}
$$
The first term is already linear in the error, thus in it we can replace $\hat{C}$ with $C$. 
In the second and third term we substitute $C \hat{C}^T \approx I + (\phi \times)$:
$$
\begin{split}
(\dot{\phi} \times) \approx -C (\Delta \omega \times) C^T - (\xi \times) (I + \phi \times) + (I + \phi \times) ((\xi + \Delta \xi) \times) =  \\\\
\approx -C (\Delta \omega \times) C^T + (\Delta \xi \times) + (\phi \times) (\xi \times) - (\xi \times) (\phi \times)
\end{split}
$$
This equation can be transformed to a vector form using easy to proof identities:
$$
C (u \times) C^T = ((C u) \times) \text{ if $C$ is a rotation matrix} \\\\
(u \times ) (v \times) - (v \times) (u \times) = ((u \times v) \times)
$$
After applying them we get a relation for skew symmetric matrices, which can be written in an equivalent vector form:
$$
\dot \phi = -C \Delta \omega + \Delta \xi - \xi \times \phi
$$

## Equation for the INS attitude error vector

Applying the derived equation to the equation for $C^n_b$ we get for its error vector
$$
\dot{\phi} = -C^n_b \Delta \omega + \Delta \rho^n + \Delta \Omega^n - (\rho^n + \Omega^n) \times \phi
$$
Substituting the angular rate errors we get
$$
\dot{\phi} = (R + (\Omega^n \times) T) \Delta r^n + T \Delta V^n - (\Omega^n + \rho^n) \times \phi -C^n_b \Delta \omega
$$

# Summary of the equations

Let's write all the equations together:
$$
\Delta \dot{r}^n = -\left[ (\rho^n \times) + (V^n \times) T \right] \Delta r^n + \Delta V^n \\\\
\Delta \dot{V}^n = \left[2 (V^n \times) (\Omega^n \times) T + (V^n \times) R + G \right] \Delta r^n + \left[ -((2 \Omega^n + \rho^n)\times) + (V^n \times) T \right] \Delta V^n + (C^n_b f^b) \times \phi + C^n_b \Delta f \\\\
\dot{\phi} = (R + (\Omega^n \times) T) \Delta r^n + T \Delta V^n - (\Omega^n + \rho^n) \times \phi -C^n_b \Delta \omega
$$
They can be implemented exactly like this using corresponding and it's probably the best practical approach.
The explicit matrix expressions as the coefficients are convenient for the implementation.

Expanding the formula for the vertical position error we get
$$
\Delta \dot{r}_D = \Delta V_D
$$
Which agrees with the simple altitude kinematics equation $\dot{h} = -V_D$.

# Simplified formulas 
To reduce symbolic and conceptual clutter we can simplify the formulas if position errors are estimated and corrected (using GNSS etc.)
In this case it can be assumed that $V \Delta r / a \ll \Delta V$.
Then in all the formulas terms with $\Delta r^n$ can be neglected in comparison with terms with $\Delta V^n$.

With the stated simplification the formulas look like this:
$$
\Delta \dot{r}^n = \Delta V^n \\\\
\Delta \dot{V}^n = G \Delta r^n + \left[ -((2 \Omega^n + \rho^n)\times) + (V^n \times) T \right] \Delta V^n + (C^n_b f^b) \times \phi + C^n_b \Delta f \\\\
\dot{\phi} = T \Delta V^n - (\Omega^n + \rho^n) \times \phi - C^n_b \Delta \omega
$$

Let's briefly discuss these equations:

* Position error equation is straightforward: higher the velocity error -- quicker position error grows
* Rate of change of the vertical velocity error is influenced mainly by incorrectly computed gravity
* Rate of change of the horizontal velocity error is influenced by incorrect heading when horizontal acceleration is applied and always by incorrect roll and pitch (because $f \approx -g$)
* Projected accelerometer error also directly contributes to the rate of change of the velocity error
* Attitude error grows approximately with a rate equal to the projected gyro error

# Practical considerations

The equations was written as if using true values of the variables to compute its coefficients.
In practice estimated values (as only available) are used for that, it doesn't violate the validity of the equations to the first order in the errors.
However it might create consistency and false observability issues in estimation algorithms and some special care must be taken in practical systems.

# Conclusion 

I've presented a derivation of standard INS error equations with all coefficients given as explicit matrices, which is convenient for implementation.
From this standard "phi-angle" parametrization other parametrizations might be derived, which I can touch upon in a follow up post.
