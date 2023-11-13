# 2: Optimization with Steady-State PDE Constraints

This example considers a [constrained optimization problem](../../problems/optimization) of the form

$$
\begin{aligned}
    \min_{\u,\z} ~& J(\u,\z)
    \\
    s.t. ~~& \c(\u,\z) = 0
\end{aligned}
$$

where the state $\u$ is the solution of a steady-state partial differential equation.

We will implement subclasses of {class}`Objective` and {class}`Constraint` and show how to solve the optimization problem with a {class}`Reduced_Space_Optimization`.

## Problem Statement

Let $u(x)$ be a function of the two-dimensional spatial variable $\mathbf{x}\in\Omega = [-1, 1] \times [-1, 1]$.
Our goal is to choose Dirichlet boundary conditions for $u$ such that the solution to a Poisson equation has the smallest possible gradient over $\Omega$.
This is formalized as the constrained optimization problem

$$
\begin{align*}
    \min_{u, z} &~ J(u, z) = \int_{\Omega}\|\nabla u\|_{2} \:d\mathbf{x}
    \\
    s.t. ~~& \Delta u = f, ~~ \mathbf{x}\in\Omega,
    \qquad
    u(\mathbf{x}) = z(\mathbf{x}), ~~ \mathbf{x}\in\partial\Omega,
\end{align*}
$$

where $z(x)$ is the spatially-dependent control, $f(\mathbf{x}) = x_{1}^{2} + x_{2}^{3}$ is given, and $\partial\Omega$ is the boundary of $\Omega$.
