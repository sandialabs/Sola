# 3: Optimization with Time-dependent ODE Constraints

This example considers a [constrained optimization problem](../../problems/optimization) of the form

$$
\begin{aligned}
    \min_{u,z} ~& J(u,z)
    \\
    s.t. ~~& c(u,z) = 0
\end{aligned}
$$

where the state $u$ is the solution of a small system of ordinary differential equations.

We will implement subclasses of [`Objective`](sabl:optimization-objective) and [`Constraint`](sabl:optimization-constraints) and show how to solve the optimization problem with a [`Reduced_Space_Optimization`](sabl:optimizer-class).

## Problem Statement

Let $u(t)$ be a function of time.

$$
\begin{align*}
    \Delta u = f, ~~ \mathbf{x}\in\Omega,
    \qquad
    u(\mathbf{x}) = z(\mathbf{x}), ~~ \mathbf{x}\in\partial\Omega.
\end{align*}
$$

where $z(x)$ is the spatially-dependent control and $f(\mathbf{x}) = x_{1}^{2} + x_{2}^{3}$ is given.
We consider the objective function

$$
\begin{align*}
    J(u, z) = \int_{\Omega}\|\nabla u\|_{2} \:d\mathbf{x}.
\end{align*}
$$

This is the problem of choosing the Dirichlet boundary conditions that result in the state with the smallest gradient.
