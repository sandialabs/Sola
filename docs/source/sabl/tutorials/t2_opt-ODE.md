# 2: Optimization with Differential Constraints

This example considers a [constrained optimization problem](../../problems/optimization) of the form

$$
\begin{aligned}
    \min_{\y(t),\z} ~& \int_{0}^{T} g(\y(t),t) dt + R(\z)
    \\
    s.t. ~~& \frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t), \z, t), ~~ \y(0) = \h(\z),
\end{aligned}
$$

where the time-dependent state $\y(t)$ is the solution to an ordinary differential equation (ODE).

We will implement subclasses of {class}`Dynamic_Objective` and {class}`Dynamic_Constraint` by explicitly calculating the derivatives of $g$, $R$, $\f$, and $\h$, then show to how to solve the optimization problem with a {class}`Reduced_Space_Optimization`.
Afterward, we use {class}`Dynamic_Objective_AD` and {class}`Dynamic_Constraint_AD` to use automatic differentiation to calculate the derivatives.

:::{note}
This tutorial includes a few short exercises and their solutions.
The finished produced is included in the SABL souce code under `tutorials/Tutorial_2/`.
:::

## Problem Statement

The following equations, derived from Newton's laws, model in polar coordinates the two-dimensional motion of a satellite orbiting a celestial body under the influence of gravity.
<!-- SOURCE: https://personal.math.ubc.ca/~loew/m403/ode-intro.pdf -->

$$
\begin{aligned}
    \frac{\textrm{d}^{2}r}{\textrm{d}t^2} - r(t)\cdot\left(\frac{\textrm{d}\theta}{\textrm{d}t}\right)^2
    &= -\frac{k}{r(t)^{2}},
    \qquad
    r(t) \cdot \frac{\textrm{d}^{2}\theta}{\textrm{d}t^2} + 2 \frac{\textrm{d}r}{\textrm{d}t}\frac{\textrm{d}\theta}{\textrm{d}t}
    = 0.
\end{aligned}
$$

Here, $t$ is time, $r(t)$ is the radius of the orbit, and $\theta(t)$ is the angular coordinate.
Defining $\y(t) = (y_1(t),y_2(t),y_3(t),y_4(t))\trp = (r(t), \frac{dr}{dt}, \theta(t), \frac{d\theta}{dt})\trp$, we obtain a first-order system of ODEs:

$$
\begin{aligned}
    \frac{\textrm{d}}{\textrm{d}t}\y(t)
    % = \frac{\textrm{d}}{\textrm{d}t}
    % \left(\begin{array}{c}
    %     y_1(t) \\ y_2(t) \\ y_3(t) \\ y_4(t)
    % \end{array}\right)
    = \f(\y(t)) := \left(\begin{array}{c}
        y_2(t) \\
        y_1(t)y_4(t)^2 - k/y_1(t)^2 \\
        y_4(t) \\
        - 2 y_2(t) y_4(t) / y_1(t)
    \end{array}\right).
\end{aligned}
$$

Our goal is to choose initial conditions such that the satellite has an orbit with a given radius $\rho$ and angular momentum $\omega$.
This leads to the following optimization problem.

$$
\begin{aligned}
    \min_{\y(t),\z} ~& \int_{0}^{T} (y_1(t) - \rho)^2 + y_2(t)^2 + (y_3(t) - \omega t)^2 + (y_4(t) - \omega)^2 dt
    \\
    s.t. ~~& \frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t), t), ~~ \y(0) = (z_1, 0, 0, z_2),
\end{aligned}
$$

where $\z = (z_1,z_2)\trp$ describes the radius and angular velocity at time $0$.
This is the constrained optimization problem given earlier with

$$
\begin{aligned}
    \g(\y(t), t) &= (y_1(t) - \rho)^2 + y_2(t)^2 + (y_3(t) - \omega t)^2 + (y_4(t) - \omega)^2,
    \\
    R(\z) &= 0,
    \qquad
    \h(\z) = \z.
\end{aligned}
$$

Note that $\g$ can also be written as $\g(\y(t), t) = \|\y(t) - \boldsymbol{\alpha}(t)\|_2^2$ where $\boldsymbol{\alpha}(t) = (\rho, 0, \omega t, \omega)\trp$.

The dimension of the ODE state is $n_y = 4$ and the dimension of the control is $n_z = 2$.
The full state $\u$ of the optimization problem consists of the ODE state at a collection of time instances, but we never need form $\u$ explicitly in our implementation.

:::{warning}
The rest of this page is under construction, please check back later.
:::
