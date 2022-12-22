---
title: "INS error equations"
date: 2022-12-20
katex: true
---

In this note I want to derive Inertial Navigation System (INS) error equations in concise and coherent fashion.
This is of course a well known subject, however it's useful to derive all the equations from scratch and revisit the subject from time to time.

# INS state variables and equations for them

An Inertial Navigation System computes position, velocity and attitude of a moving object using gyro and accelerometer measurements from an Inertial Measurement Unit (IMU).

To formulate necessary equations, frames of reference (denoted by lowercase letters) need to be defined:

* $e$ -- Earth-centered, Earth-fixed (ECEF) frame 
* $i$ -- inertial frame with axes coinciding with frame $e$ at the initial time
* $n$ -- locally horizontal North-East-Down (NED) frame
* $b$ -- body frame with axes along inertial sensor axes

A vector $u$ expressed in a frame $x$ is denoted as $u^x$.
A rotation matrix $C^x_y$ projects a vector as $u^x = C^x_y u^y$.

State variables computed by INS are:

* Geographical position as latitude, longitude and altitude -- $\varphi, \lambda, h$
* Velocity relative to the ground (frame $e$) expressed in frame $n$ -- $V^n = \begin{bmatrix} V_N & V_E & V_D \end{bmatrix}^T$
* Attitude as a rotation matrix from frame $n$ to frame $b$ -- $C^n_b$

The differential equations for the state variables are:
$$
\dot{\varphi} = \frac{V_N}{R_N + h} \\\\[1 pt]
\dot{\lambda} = \frac{V_E}{(R_E + h) \cos \varphi} \\\\
\dot{h} = -V_D \\\\
\dot{V}^n = -(\omega_{ie}^n + \omega_{in}^n) \times V^n + g^n + C^n_b f^b \\\\
\dot{C}^n_b = C^n_b [\omega_{ib}^b \times] - [\omega_{in}^n \times] C^n_b
$$
Where the following variables and notation is used:

* $R_N$ and $R_E$ -- North and East curvature radii
* $\omega_{xy}^z$ -- angular rate of frame $y$ relative to frame $x$ expressed in frame $z$:    
    * $\omega_{ie}^n$ -- Earth angular rate
    * $\omega_{in}^n = \omega_{ie}^n + \omega_{en}^n$ where $\omega_{en}^n$ is transport rate i. e. angular rate of frame $n$ relative to frame $e$ due to object's velocity
    * $\omega_{ib}^b$ -- angular rate of the body frame relative to the inertial frame measured by gyroscopes (error free)
* $g^n$ -- gravity vector
* $f^b$ -- specific force measured by accelerometers (error free)
* $[u \times]$ -- cross-product skew symmetric matrix, i. e. $[u \times] v \equiv u \times v$

The INS computer updates the state variables according to the equations using gyro and accelerometers measurements containing errors:
$$
\tilde{\omega}\_{ib}^b = \omega_{ib}^b + \epsilon^b \\\\
\tilde{f}^b = f^b + \nabla^b
$$
Where gyro and accelerometer errors are denoted by $\epsilon^b$ and $\nabla^b$ respectively.
The computed or estimated variables will be marked with hat.

# Three realizations of navigation frame $n$

In INS operation there are 3 realization of navigation frame $n$ (NED frame) which are slightly misaligned to each other:

1. True navigation frame $t$ -- NED frame at a true INS location
2. Computed navigation frame $c$ -- NED frame at a location computed by INS
3. Platform navigation frame $p$ -- frame where accelerometer measurements are projected using estimated matrix $\hat{C}^n_b$
   In gimbaled INS this frame is physically realized by a gyro-stabilized platform on which accelerometers are installed

From the mathematical point of view these frame can be defined as
$$
C^e_t \equiv C^e_n \\\\
C^e_c \equiv \hat{C}^e_n \\\\
C^p_b \equiv \hat{C}^n_b
$$

Small misalignments between the frames are described by rotation vectors according to the following scheme:
$$
\begin{array}{ccccc}
  & \theta          &   & \psi   &     \\\\
t & \longrightarrow & c & \longrightarrow & p \\\\
\end{array} \\\\
\overrightarrow{\phantom{xxx|}\phi = \theta + \psi\phantom{xxx|}}
$$

We can give the following interpretation of the misalignment vectors:

* $\phi$ is a full attitude error vector
* $\theta$ is a misalignment vector between true and computed $n$ frames due to the INS position error
* $\psi$ is a part of the attitude error ignoring computed and true $n$ frames misalignment

Rotation matrices between these system can be expressed to the first order of the introduced rotation vectors as:
$$
C^t_c \approx I + [\theta \times] \\\\
C^c_p \approx I + [\psi \times] \\\\
C^t_p \approx I + [\phi \times]
$$

When the distinction between frames $t$, $c$ and $p$ are not important (for example in coefficients of error equations), we will denote the frame as $n$.

# INS equations for an arbitrary frame

The velocity and attitude equations can be written for an arbitrary frame $x$ as 
$$
\dot{V}^x = -(\omega_{ie}^x + \omega_{ix}^x) \times V^x + g^x + C^x_b f_b \\\\
\dot{C}^x_b = C^x_b [\omega_{ib}^b \times] - [\omega_{ix}^x  \times] C^x_b
$$
Equations for the true velocity and attitude come out with $x = t$.
In the following these equations will be used for frames $c$ and $p$ as well.

# Error definitions

Errors must express deviation of estimated variables from true variables.
They don't need to be direct differences of the variables.

## Position error 

Direct differences of latitude and longitude are not particularly convenient to interpret as well as to derive equations for.
Instead we can use correspondence between latitude, longitude, altitude and coordinates in ECEF, which we denote as $r^e$.

Define the position error as
$$
\Delta r = C^n_e (\hat{r}^e - r^e)
$$
That is it ECEF coordinate errors projected into NED frame.

Relations between $\Delta r$ and geographic coordinates errors are
$$
\Delta r = \begin{bmatrix}
(R_N + h) \Delta \varphi \\\\
(R_E + h) \cos \varphi \Delta \lambda \\\\
-\Delta h
\end{bmatrix}
$$

## Velocity error

A velocity error can be defined relative to the true velocity expressed in true, computed or platform frames:
$$
\Delta V = \hat{V} - V^t \\\\
\partial V = \hat{V} - V^c \\\\
\delta V = \hat{V} - V^p
$$
Symbol $\partial$ was used simply for the lack of better distinguishable alternatives.

All these errors are related by a linear transform. 
For $\partial V$ we have:
$$
\partial V = \hat{V} - V^c = \hat{V} - C^c_t V^t \approx \hat{V} - [I - (\theta \times)] V^t = \Delta V + \theta \times V^n
$$
And analogously for $\delta V$:
$$
\delta V \approx \Delta V + \phi \times V^n
$$

## Attitude error

As an attitude error either vector $\phi$ or $\psi$ can be used.
When $\psi$ is used the full attitude error is computed as $\phi = \theta + \psi$ where $\theta$ is computed from the position error as explained next. 

# Error equations in general

Error equations are linear differential equations for position, velocity and attitude errors with time varying coefficients, which depend on position, velocity and attitude of INS.
Gyro and accelerometers error also appear in right hand side of the equations.
They are used in estimation algorithms for INS like linearized Kalman filter, Extended Kalman Filter, etc.

Linearity of the equations are achieved by a linearization process, meaning that the equations are approximate in nature and valid only as long as the errors are small.
The exact meaning of the errors being small will be given in the end.

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
\omega_{en}^n =  \begin{bmatrix}
\dfrac{V_E}{R_E + h} \\\\
-\dfrac{V_N}{R_N + h} \\\\
-\dfrac{V_E \tan \varphi}{R_E + h}
\end{bmatrix} = T V^n
$$
It is used in the kinematics equation for $C^e_n$:
$$
\dot{C}^e_n = C^e_n [\omega_{en}^n \times]
$$

The matrix $T$ relates a linear displacement in $n$ frame to its small angular rotation.
Perturbation of the transport rate due to position and velocity errors can be expressed as
$$
\Delta \omega_{en}^n \equiv \hat{\omega}\_{en}^n - \omega_{en}^n \equiv \omega_{ec}^c - \omega_{et}^t \approx R \Delta r + T \Delta V
$$


Where Jacobian matrices $T$ and $R$ are defined as
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

By the same logic we can relate the vector $\theta$ with the position error as
$$
\theta \approx T \Delta r
$$

The Earth rotation rate in $n$ frame has the form
$$
\omega_{ie}^n = C^n_e \omega_{ie}^e \text{ with } \omega_{ie}^e = \begin{bmatrix} 0 & 0 & \Omega \end{bmatrix}^T
$$
Its perturbation due to the position error can be expressed as 
$$
\Delta \omega_{ie}^n \equiv \hat{\omega}\_{ie}^n - \omega_{ie}^n \equiv \omega_{ie}^c - \omega_{ie}^t \approx \omega_{ie}^n \times \theta \approx (\omega_{ie}^n \times) T \Delta r^n
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
Perturbation of the gravity vector due to the position error is
$$
\Delta g^n \equiv \hat{g}^n - g^n \approx G \Delta r
$$
with the Jacobian matrix
$$
G = \begin{bmatrix}
0 & 0 & 0 \\\\
0 & 0 & 0 \\\\
0 & 0 & \dfrac{2 g}{a}  \\\\
\end{bmatrix}
$$

# Differential equation for a small rotation vector as attitude error

Consider the rotation kinematics equation for a matrix $C$:
$$
\dot{C} = C [\omega \times] - [\xi \times] C
$$
Where $\omega$ and $\xi$ are angular rates of the frames between which $C$ transforms relative to some third frame.
And the same dynamics for its estimate:
$$
\dot{\hat{C}} = \hat{C} [(\omega + \Delta \omega) \times] - [(\xi + \Delta \xi)  \times] \hat{C}
$$
Define the relation between the matrices to the first order in a rotation vector $\alpha$ as 
$$
C \approx (I + [\alpha \times]) \hat{C}
$$
From this we get the differential equation for $\alpha$:
$$
\begin{split}
[\dot{\alpha} \times] = \dot{C} \hat{C}^T + C \dot{\hat{C}}^T = (C [\omega \times] - [\xi \times] C ) \hat{C}^T + C(-[(\omega + \Delta \omega) \times] \hat{C}^T + \hat{C}^T [(\xi + \Delta \xi) \times]) = \\\\
= -C [\Delta \omega \times] \hat{C}^T - [\xi \times] C \hat{C}^T  + C \hat{C}^T [(\xi + \Delta \xi) \times]
\end{split}
$$
The first term is already linear in the error, thus in it we can replace $\hat{C}$ with $C$. 
In the second and third terms we substitute $C \hat{C}^T \approx I + [\alpha \times]$:
$$
\begin{split}
[\dot{\alpha} \times] \approx -C [\Delta \omega \times] C^T - [\xi \times] (I + [\alpha \times]) + (I + [\alpha \times]) [(\xi + \Delta \xi) \times] =  \\\\
\approx -C [\Delta \omega \times] C^T + [\Delta \xi \times] + [\alpha \times] [\xi \times] - [\xi \times] [\alpha \times]
\end{split}
$$
This equation can be transformed to a vector form using easy to proof identities:
$$
C [u \times] C^T = [(C u) \times] \text{ if $C$ is a rotation matrix} \\\\
$$
$$
[u \times][v \times] - [v \times][u \times] = [(u \times v) \times]
$$

After applying them we get a relation for skew symmetric matrices, which can be written in an equivalent vector form:
$$
\dot \alpha = -C \Delta \omega + \Delta \xi - \xi \times \alpha
$$

# $\phi$-angle error model

In this error model $\Delta r, \Delta V$ and $\phi$ are used as the error variables.
Conceptually this is the most straightforward approach where error equations are essentially derived as perturbation of nominal INS equations.

## Position error equation

Differentiating the position error we get
$$
\begin{split}
\Delta \dot{r} = \dot{C}^n_e (\hat{r}^e - r^e) + C^n_e (\dot{\hat{r}}^e - \dot{r}^e) = -[\omega_{en}^n \times] C^n_e (\hat{r}^e - r^e) + C^n_e (\hat{V}^e - V^e) = \\\\
= -\omega_{en}^n \times \Delta r + C^n_e \hat{C}^e_n \hat{V}^n - V^n \approx -\omega_{en}^n \times \Delta r + (I + [\theta \times]) (V^n + \Delta V^n) - V^n \approx  \\\\
\approx -\omega_{en}^n \times \Delta r + \Delta V + \theta \times V^n
\end{split}
$$

Substituting $\theta = T \Delta r$ we get:
$$
\Delta \dot{r} = -\left( [\omega_{en}^n \times] + [V^n \times] T \right) \Delta r + \Delta V
$$

## Velocity error equation

To derive the velocity error equation we subtract the equation for the true velocity from the equation for the INS estimated velocity keeping terms linear in the errors:
$$
\Delta \dot{V} = -(2 \omega_{ie}^n + \omega_{en}^n) \times \Delta V - (2 \Delta \omega_{ie}^n + \omega_{en}^n) \times V^n + \Delta g^n + (I - [\phi \times]) C^n_b (f^b + \nabla^b) - C^n_b f^b
$$
Substituting expressions for all $\Delta$-s we get
$$
\Delta \dot{V} = -(2 \omega_{ie}^n + \omega_{en}^n) \times \Delta V + [V^n \times] (2 [\omega_{ie}^n \times] T \Delta r^n + R \Delta r^n + T \Delta V^n) + G \Delta r^n + (C^n_b f^b) \times \phi + C^n_b \nabla^b
$$
And after the grouping the terms:
$$
\Delta \dot{V} = (2 [V^n \times] [\omega_{ie}^n \times] T + [V^n \times] R + G ) \Delta r + (-[(2 \omega_{ie}^n + \omega_{en}^n)\times] + [V^n \times] T) \Delta V + (C^n_b f^b) \times \phi + C^n_b \nabla^b
$$

## Attitude error equation

Applying the general equation for an attitude error vector to the matrix $C^n_b$ we get
$$
\dot{\phi} = -C^n_b \epsilon^b + \Delta \omega_{en}^n + \Delta \omega_{ie}^n - (\omega_{ie}^n + \omega_{en}^n) \times \phi
$$
And after Substituting the angular rate errors:
$$
\dot{\phi} = (R + [\omega_{ie}^n \times] T) \Delta r + T \Delta V - (\omega_{ie}^n + \omega_{en}^n) \times \phi -C^n_b \epsilon^b
$$

# $\psi$-angle error model

In this error model $\Delta r, \partial V$ and $\psi$ are used as the error variables.
This model is often considered as an alternative to $\phi$-angle model and thus I decided to describe it.

## Position error equation

To derive the position error equations we take the equation from the $\phi$-model:
$$
\Delta \dot{r} \approx -\omega_{en}^n \times \Delta r + \Delta V + \theta \times V^n \approx -\omega_{en}^n \times \Delta r + \partial V
$$
We can see that $\partial V$ is in some sense the most natural variable to describe the position error growth.

## Velocity error equation

Let's write differential equations for $\hat{V}$ and $V^c$ using the common notation:
$$
\dot{\hat{V}} = -(2 \omega_{ie}^c + \omega_{ec}^c) \times \hat{V} + \hat{g}^n + C^p_b (f^b + \nabla^b) \\\\
\dot{V}^c = -(2 \omega_{ie}^c + \omega_{ec}^c) \times V^c + g^c + C^c_b f^b
$$
Also consider relations
$$
C^p_b \approx (I - [\psi \times]) C^c_b \\\\
\hat{g}^n = g^n + \Delta g^n \\\\
g^c = C^c_t g^t \approx (I - [\theta \times]) g^n = g^n + g^n \times \theta \\\\
$$
Taking the difference of the equations and using the above relations we get to the first order of the errors:
$$
\partial \dot{V} = -(2 \omega_{ie}^n + \omega_{en}^n) \times \partial V + \Delta g^n - g^n \times \theta + (C^n_b f^b) \times \psi + C^n_b \nabla^b
$$
Substitute expressions for $\Delta g^n$ and $\theta$ to finally get
$$
\partial \dot{V} = (G - [g^n \times] T) \Delta r - (2 \omega_{ie}^n + \omega_{en}^n) \times \partial V + (C^n_b f^b) \times \psi + C^n_b \nabla^b
$$

## Attitude error equation

Let's write differential equations for $\hat{C}^n_b$ and $\hat{C}^c_b$ which are misaligned by vector $\psi$:
$$
\dot{\hat{C}}^n_b = \hat{C}^n_b [(\omega_{ib}^b + \epsilon^b) \times] - [\omega_{ic}^c \times] \hat{C}^n_b \\\\
\dot{C}^n_b = \hat{C}^n_b [\omega_{ib}^b \times] - [\omega_{ic}^c \times] C^c_b
$$
We see that the second angular velocity $\omega_{ic}^c$ is the same in both equations.
Thus using the general equation for an attitude error vector we get
$$
\dot{\psi} = -(\omega_{ie}^n + \omega_{en}^n) \times \psi - C^n_b \epsilon^b
$$

# Modified $\phi$-vector error model

In this error model $\Delta r, \delta V$ and $\phi$ are used as the error variables.
This model was not considered by early researchers, but was discovered later and found to be the most advantageous to use in practice.

## Formula for $\omega_{ip}^p$

In the following we will need an expression for $\omega_{ip}^p$.
To derive it we write the equation for $C^p_b$ in two different ways:

* As executed by the INS computer: $\dot{C}^p_b = C^p_b [(\omega_{ib}^b + \epsilon^b)] - [\omega_{ic}^c \times] C^p_b$
* Using the general equation for a frame $p$: $\dot{C}^p_b = C^p_b [\omega_{ib}^b \times] - [\omega_{ip}^p \times] C^p_b$

Equating right hand sides we get:
$$
C^p_b [(\omega_{ib}^b + \epsilon^b) \times] - [\omega_{ic}^c \times] C^p_b = C^p_b [\omega_{ib}^b \times] - [\omega_{ip}^p \times] C^p_b \\\\
C^p_b [\epsilon^b \times] (C^p_b)^T - [\omega_{ic}^c \times ] = -[\omega_{ip}^p \times] \\\\
\omega_{ip}^p = \omega_{ic}^c - C^n_b \epsilon^b
$$

## Position error equation derivation

Taking the equation from the $\phi$-error model and expressing $\Delta V$ via $\delta V$ we get
$$
\Delta \dot{r} = -\left( [\omega_{en}^n \times] + [V^n \times] T \right) \Delta r + \delta V + V^n \times \phi
$$
Here we see that the attitude and position errors are directly coupled.

# Velocity error derivation

Let's write differential equations for $\hat{V}$ and $V^p$ using the common notation:
$$
\dot{\hat{V}} = -(\omega_{ie}^c + \omega_{ic}^c) \times \hat{V} + \hat{g}^n + C^p_b (f^b + \nabla^b) \\\\
\dot{V}^p = -(\omega_{ie}^p + \omega_{ip}^p) \times V^p + g^p + C^p_b f^b
$$
Also consider relations
$$
\hat{g}^n = g^n + \Delta g^n \\\\
g^p = C^p_t g^t \approx (I - [\phi \times]) g^n = g^n + g^n \times \phi \\\\
\omega^p_{ie} = C^p_c \omega^c_{ie} \approx (I - [\psi \times]) \omega^c_{ie} = \omega^c_{ie} + \omega^c_{ie} \times \psi \\\\
\omega^p_{ip} = \omega^c_{ic} - C^n_b \varepsilon^b
$$
Taking the difference of the velocity equations and using the above relations we get to the first order in the errors:
$$
\delta V = V^n \times (-\omega_{ie}^n \times \psi + C^n_b \epsilon^b) - (2 \omega_{ie}^n + \omega_{en}^n) \times \delta V^n + \Delta g^n - g^n \times \phi + C^n_b \nabla^b
$$
Substituting $\psi = \phi - \theta \approx \phi - T \Delta r, \Delta g^n \approx G \Delta r$ and grouping the terms we finally get
$$
\delta \dot{V} = (G + [V^n \times][\omega_{ie}^n \times] T ) \Delta r - (2 \omega_{ie}^n + \omega_{en}^n) \times \delta V^n - ([g^n \times] + [V^n \times][\omega_{ie}^n \times]) \phi + C^n_b \nabla^b + [V^n \times] C^n_b \epsilon^b
$$
The most notable property of this equation is absence of specific force measurements $f^b$.

## Attitude error equations

Start from the equation already derived:
$$
\dot{\phi} = (R + [\omega_{ie}^n \times] T) \Delta r + T \Delta V - (\omega_{ie}^n + \omega_{en}^n) \times \phi -C^n_b \epsilon^b
$$
Expressing $\Delta V$ via $\delta V$ we get
$$
\dot{\phi} = (R + [\omega_{ie}^n \times] T) \Delta r + T \delta V + (T [V^n \times] - [(\omega_{ie}^n + \omega_{en}^n) \times]) \phi -C^n_b \epsilon^b
$$

# Validity of error models

The derivations were done using several first order approximations.
These approximations are valid if higher order truncated terms are insignificant compared to the first order (linear) parts.

For the attitude error vectors we've used linear approximations like $C^t_p \approx I + [\phi \times]$. 
Thus it is obvious that the requirement for them are
$$
\phi \ll 1 \\\\
\psi \ll 1
$$
A working limit on the attitude errors for algorithms like EKF is approximately 0.1 radian or 5 degrees.

We've also used an assumption that the vector $\theta$ is small. It means that 
$$
\Delta r \ll a
$$
Considering that Earth radius $a \approx 6400$ km, we can estimate the limit on acceptable position error as 500 km.
This value far exceeds useful range of INS accuracy.

The situation with the velocity error is more interesting.
In $\psi$-angle error model there is no truncation of higher order terms featuring $\partial V$.
This is an interesting property which could also be derived or understood by considering kinematic equations in ECEF frame.
However if we try to predict $\Delta V \approx \partial V + V^n \times \theta $ or $\delta V \approx \partial V - V^n  \times \psi$ errors (necessary for NED or body frame velocity processing in estimation algorithm), we will need to use the estimated velocity $\hat{V} = V^n + \Delta V$.
And then for $\hat{V} \times \theta$ and $\hat{V} \times \psi$ to be valid to the first order it must hold
$$
\Delta V \ll V^n
$$
For $\phi$-angle and its modified version other (but conceptually similar) arguments might be made to come up with the same conclusion.

# Summary and discussion

Let's summarize all the models, discuss their possible simplifications, properties and relations to each other.

## $\phi$-angle error model

$$
\Delta \dot{r} = -\left( [\omega_{en}^n \times] + [V^n \times] T \right) \Delta r + \Delta V \\\\
\Delta \dot{V} = (2 [V^n \times] [\omega_{ie}^n \times] T + [V^n \times] R + G ) \Delta r + (-[(2 \omega_{ie}^n + \omega_{en}^n)\times] + [V^n \times] T) \Delta V + (C^n_b f^b) \times \phi + C^n_b \Delta f \\\\
\dot{\phi} = (R + [\omega_{ie}^n \times] T) \Delta r + T \Delta V - (\omega_{ie}^n + \omega_{en}^n) \times \phi -C^n_b \epsilon^b
$$

Expanding the equation for the vertical position error we'll get
$$
\Delta r_D = \Delta V_D
$$
which agrees with the simple kinematics equations for the altitude $\dot{h} = -V_D$.

The formulas can be simplified if position errors are estimated and corrected (using GNSS, etc.)
In this case it can be assumed that $V \Delta r / a \ll \Delta V$.
Then in all the formulas terms with $\Delta r$ can be neglected in comparison with terms with $\Delta V$ (except $G \Delta r$!).

## $\psi$-angle error model

$$
\Delta \dot{r} = -\omega_{en}^n \times \Delta r + \partial V \\\\
\partial \dot{V} = (G - [g^n \times] T) \Delta r - (2 \omega_{ie}^n + \omega_{en}^n) \times \partial V + (C^n_b f^b) \times \psi + C^n_b \nabla^b \\\\
\dot{\psi} = -(\omega_{ie}^n + \omega_{en}^n) \times \psi - C^n_b \epsilon^b
$$

These equations already look simple enough which is probably the main advantage of this model.

## Modified $\phi$-angle error model

$$
\Delta \dot{r} = -\left( [\omega_{en}^n \times] + [V^n \times] T \right) \Delta r + \delta V + V^n \times \phi \\\\
\delta \dot{V} = (G + [V^n \times][\omega_{ie}^n \times] T ) \Delta r - (2 \omega_{ie}^n + \omega_{en}^n) \times \delta V^n - ([g^n \times] + [V^n \times][\omega_{ie}^n \times]) \phi + C^n_b \nabla^b + [V^n \times] C^n_b \epsilon^b \\\\
\dot{\phi} = (R + [\omega_{ie}^n \times] T) \Delta r + T \delta V + (T [V^n \times] - [(\omega_{ie}^n + \omega_{en}^n) \times]) \phi -C^n_b \epsilon^b
$$

If desired the equations can be simplified with assumption of $V \Delta r / a \ll \Delta V$ and $V \Omega \ll g$ as
$$
\Delta \dot{r} = \delta V + V^n \times \phi \\\\
\delta \dot{V} = (G + [V^n \times][\omega_{ie}^n \times] T ) \Delta r - (2 \omega_{ie}^n + \omega_{en}^n) \times \delta V^n - g^n \times \phi + C^n_b \nabla^b + [V^n \times] C^n_b \epsilon^b \\\\
\dot{\phi} = T \delta V + (T [V^n \times] - [(\omega_{ie}^n + \omega_{en}^n) \times]) \phi -C^n_b \epsilon^b
$$


## Discussion of $\phi$ and $\psi$-angle error models

$\phi$-angle $\psi$-angle error models are not substantially different and represent the following essential error dynamics:

* Position error equation is straightforward: higher the velocity error -- quicker position error grows
* Rate of change of the vertical velocity error is influenced mainly by incorrectly computed gravity because of the altitude error
* Roll and pitch errors create incorrect gravity projection and horizontal velocity error growth
* When horizontal acceleration is present, incorrect heading ($\phi_D, \psi_D$) causes horizontal error growth
* Projected accelerometer errors also directly contribute to the rate of change of the velocity error
* Attitude errors grow approximately with a rate equal to the projected gyro errors

In $\psi$ model roll and pitch errors are formed by horizontal $\psi$ components and horizontal position errors.
There is an oscillation mode of roll, pitch, horizontal position and velocity errors called Schuller oscillations with period of approximately 84 minutes.
Unlike the $\phi$ vector, the $\psi$ vector is free from these oscillations.

Equations of $\psi$-angle error model are simpler, however the indirect and tricky nature of $\psi$ vector makes them not appealing for use in practice in my opinion.
For example when the attitude of a stationary INS is estimated by averaging accelerometers (roll and pitch estimation) and gyroscopes (heading estimation by gyrocompassing), we can determine covariance of $\phi$ angle and its correlation with the inertial sensor biases.
But if $\psi$ angle is used we need to account for the interplay between $\phi$ and $\theta$ vectors to get the initial covariance of $\psi$ angle and correlations with other states.

## Discussion of modified $\phi$-angle error model

Modified $\phi$-angle error model differ substantially from the two other models.
In this model we can though of the estimated velocity as being expressed in the platform frame.
The misalignment of the platform frame from the true frame causes incorrect velocity projection when updating position and the position error growth as $V^n \times \phi$.
The velocity estimates in the platform frame also contains error $\delta V$, influenced only by accelerometer errors and incorrect projection of gravity (mainly).
Measurements of odometric type (odometer, doppler radar, etc.) allow estimation of only $\delta V$. 
That is in this model we get a nice separation of the full velocity error into odometric part $\delta V$ and misalignment part $V^n \times \phi$.

Another big advantage is that it doesn't require specific force measurements $f^b$ for its implementation.
It allows modelling of INS errors using only trajectory data and alleviates the question of required time step for error equations time propagation as $f^b$ can be a wide-band signal.

This model is also related to a modern viewpoint of the estimation on Lie groups. 
From this perspective the combination of velocity and attitude can be thought as elements of $SE(3)$ group and its errors as Lie algebra elements (defined on the left).
Using this perspective correction of velocity and attitude by the estimated errors can be done coherently using special Lie algebra elements exponentiation.

# Practical considerations

The equations was written as if using true values of the variables to compute its coefficients.
In practice estimated values (as only available) are used for that, it doesn't violate the validity of the equations to the first order in the errors.
However it might create consistency issues in estimation algorithms and some special care must be taken in practical systems.

For example when $V^n = 0$ then the small velocity error condition $\Delta V \ll V^n$ doesn't hold.
It means that we must assure that zero velocity is substituted in error propagation and measurements model coefficients in this case.
The problematic term in the modified $\phi$-angle error model is $V^n \times \phi$ in the position error equation.
We must assure that zero velocity is substituted in it to avoid creating false coupling between position and attitude error terms.

As for the model selection, the modified $\phi$-angle error model seems to be the most advantageous due to its useful conceptual and practical properties described above.

# Conclusion 

I've presented a derivation of 3 different INS error models in a uniform and concise fashion, discussed their essential properties and differences between each other.
Coefficients of the error equations were given as explicit matrices, which can be used directly for estimation algorithms implementations.
