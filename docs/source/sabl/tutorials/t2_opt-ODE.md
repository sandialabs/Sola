# 2: Optimization with Differential Constraints

This tutorial considers two [constrained optimization problems](../../problems/optimization) of the form

$$
\begin{aligned}
    \min_{\y(t),\z} ~& \int_{0}^{T} g(\y(t),t) dt + R(\z)
    \\
    s.t. ~~& \frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t), \z, t), ~~ \y(0) = \h(\z),
\end{aligned}
$$

where the time-dependent state $\y(t)$ is the solution to an ordinary differential equation (ODE).
In the first problem, the control $\z$ is a low-dimensional vector; in the second, $\z$ gathers the values of a time-dependent function.

In each case, we will implement subclasses of {class}`Dynamic_Objective` and {class}`Dynamic_Constraint` by explicitly calculating the derivatives of $g$, $R$, $\f$, and $\h$.
<!-- Afterward, we use {class}`Dynamic_Objective_AD` and {class}`Dynamic_Constraint_AD` to use automatic differentiation to calculate the derivatives. -->

:::{note}
This tutorial includes a few short exercises and their solutions.
The finished produced is included in the SABL souce code under `tutorials/Tutorial_2/`.
:::

## Example 1: Selecting Initial Conditions

We start with a simple problem for demonstration purposes, then proceed to a more interesting setup later.

### Problem Statement

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
Defining

$$
\begin{aligned}
    \y(t)
    = \left(\begin{array}{c}
        y_1(t) \\ y_2(t) \\ y_3(t) \\ y_4(t)
    \end{array}\right)
    = \left(\begin{array}{c}
        r(t) \\ dr / dt \\ \theta(t) \\ d\theta / dt
    \end{array}\right),
\end{aligned}
$$

we obtain a first-order system of ODEs:

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

Our goal is to choose initial conditions such that the satellite attains an orbit with a given constant radius $\rho$ and constant angular momentum $\omega$.
This leads to the following optimization problem.

$$
\begin{aligned}
    \min_{\y(t),\z} ~& \int_{0}^{T} (y_1(t) - \rho)^2 + y_2(t)^2 + (y_3(t) - \omega t)^2 + (y_4(t) - \omega)^2 dt
    \\
    s.t. ~~& \frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t)), ~~ \y(0) = (z_1, 0, 0, z_2),
\end{aligned}
$$

where $\z = (z_1,z_2)\trp$ describes the radius and angular velocity at time $0$.
This is the constrained optimization problem given earlier with

$$
\begin{aligned}
    g(\y(t), t)
    &= \|\y(t) - \boldsymbol{\alpha}(t)\|_2^2,
    % = (y_1(t) - \rho)^2 + y_2(t)^2 + (y_3(t) - \omega t)^2 + (y_4(t) - \omega)^2,
    \qquad
    R(\z) = 0,
    \qquad
    \h(\z) = \z,
\end{aligned}
$$

where $\boldsymbol{\alpha}(t) = (\rho, 0, \omega t, \omega)\trp$.

The dimension of the ODE state is $n_y = 4$ and the dimension of the control is $n_z = 2$.
The full state $\u$ of the optimization problem consists of the ODE state at a collection of time instances, but we never need form $\u$ explicitly in our implementation.

### Implementing the Objective

To implement the abstract methods of {class}`Dynamic_Objective` for this problem, we need to calculate the $\y$ gradient and the action of the $(\y,\y)$ Hessian of $g$, as well as the $\z$ gradient and the $(\z,\z)$ Hessian action of $R$.
Since $R(\z) = 0$ in this problem, we immediately have $\grad{z}R(\z) = \0$ and $\grad{z,z}R(\z)\v_z = \0$ for all $\v_z\in\R^{n_z}$.
To form the gradient of $g$, we calculate the derivatives of $g$ with respect to the entries of $\y$.
Writing $\boldsymbol{\alpha}(t) = (\alpha_1(t),\ldots,\alpha_{n_y}(t))\trp$,

$$
\begin{aligned}
    \frac{\partial}{\partial y_j}g(\y, t)
    = \frac{\partial}{\partial y_j}\left[\sum_{k=1}^{n_y}(y_k - \alpha_k(t))^2\right]
    = 2(y_j - \alpha_j(t)).
\end{aligned}
$$

Hence, $\grad{y}g(\y(t), t) = 2(\y(t) - \boldsymbol{\alpha}(t))$.

::::{admonition} Exercise 1
:class: exercise

Calculate the Hessian action $\grad{y,y}g(\y(t), t)\v_y$ where $\v_y\in\R^{n_y}$.

:::{admonition} Solution
:class: solution dropdown

First, calculate the second derivatives of $g$ with respect to the entries of $\y$:

$$
\begin{aligned}
    \frac{\partial}{\partial y_i\partial y_j}g(\y(t), t)
    = \frac{\partial}{\partial y_i}\left[2(y_j - \alpha_j(t))\right]
    = 2\delta_{ij}
    = \begin{cases}
    2 & \textrm{if}~i = j, \\
    0 & \textrm{else}.
    \end{cases}
\end{aligned}
$$

This tells us that $\grad{y,y}g(\y(t), t) = 2\I$, where $\I\in\R^{n_y\times n_y}$ is the identity matrix.
Therefore,

$$
\begin{aligned}
    \grad{y,y}g(\y(t), t)\v_y
    = 2\v_y.
\end{aligned}
$$
:::
::::

We can now define a subclass of {class}`Dynamic_Objective` starting from its [inheritance template](optimization.Dynamic_Objective.template).
Our new class, `Tutorial_2_Objective`, has a constructor for setting the target radius $\rho$ and angular velocity $\omega$.
This new class must be defined in a file named `Tutorial_1_Objective.m`.
In the code, the input vectors $\v_y\in\R^{n_y}$ and $\v_z\in\R^{n_z}$ are called `y_in` and `z_in`, respectively, to remind the user of the vector size.
Similarly, the output vectors `y_out` and `z_out` have $n_y$ and $n_z$ rows, respectively.

```matlab
classdef Tutorial_2_Objective < Dynamic_Objective

    properties
        radius       % Target orbital radius.
        velocity     % Target angular velocity.
    end

    methods

        function this = Tutorial_2_Objective(T, n_t, radius, velocity)
            this = this@Dynamic_Objective(4, 2, T, n_t);
            this.radius = radius;
            this.velocity = velocity;
        end

        function [val, grad_y] = g(this, y, t)
            w = this.velocity;
            r = this.radius;
            a = [r; 0; w * t; w];
            y_minus_a = y - a;

            val = sum(y_minus_a.^2);
            grad_y = 2 * y_minus_a;
        end

        function [val, grad_z] = R(this, z)
            val = 0;
            grad_z = zeros(this.n_z, 1);
        end

        function [y_out] = g_yy_Apply(this, y_in, y, t)
            y_out = error('g_yy_Apply() not implemented');
        end

        function [z_out] = R_zz_Apply(this, z_in, z)
            z_out = zeros(this.n_z, size(z_in, 2));
        end

    end
end
```

:::{danger}
Note carefully that `R_zz_Apply()` is implemented in a vectorized fashion by treating the input `z_in` as a matrix with possibly more than one column.
:::

::::{admonition} Exercise 2
:class: exercise

Implement, in a vectorized fashion, `g_yy_Apply()` for computing $\grad{y,y}g(\y(t), t)\v_y$.

:::{admonition} Solution
:class: solution dropdown

```matlab
function [y_out] = g_yy_Apply(this, y_in, y, t)
    y_out = 2 * y_in;
end
```

:::
::::

### Implementing the Constraint

To implement a {class}`Dynamic_Constraint` subclass for this problem, we need the $\z$ Jacobian and the action of the $(\z,\z)$ Hessian of $\h$, as well as the $\y$ and $\z$ Jacobians of $\f$ and its Hessian actions.
Recall that $\h(\z) = (z_1, 0, 0, z_2)\trp$ prescribes the initial condition of the ODE system.
Its $\z$ Jacobian is given by

$$
\begin{aligned}
    \h_z(\z)
    = \left(\begin{array}{cc}
        1 & 0 \\ 0 & 0 \\ 0 & 0 \\ 0 & 1
    \end{array}\right).
\end{aligned}
$$

The Hessian $\h_{z,z}(\z)$ is zero, hence the Hessian action is $\bflambda\trp\h_{z,z}(\z)\v_z = 0$ for all $\bflambda\in\R^{n_y}$ and $\v_z\in\R^{n_z}$.
Furthermore, since $\f$ does not depend on $\z$ in this problem, we immediately have that $\f_z(\y(t), \z, t)$ is the $n_y\times n_z$ zero matrix and $\bflambda\trp\f_{y,z}(\y(t), \z, t)\v_z = \bflambda\trp\f_{z,y}(\y(t), \z, t)\v_y = \bflambda\trp\f_{z,z}(\y(t), \z, t)\v_z = 0$.

The $\y$ Jacobian of $\f$ is nontrivial, but easy to calculate:

$$
\begin{aligned}
    \f_y(\y, \z, t)
    = \left(\begin{array}{cccc}
        0 & 1 & 0 & 0 \\
        y_4^2 - \frac{2k}{y_1^3} & 0 & 0 & 2y_1 y_4 \\
        0 & 0 & 0 & 1 \\
        \frac{y_2 y_4}{y_1^2} & -\frac{y_4}{y_1} & 0 & -\frac{y_2}{y_1} \\
    \end{array}\right).
\end{aligned}
$$

::::{admonition} Exercise 3
:class: exercise

Calculate the Hessian action $\bflambda\trp\f_{y,y}(\y(t),\z,t)\v_y$.

:::{admonition} Solution
:class: solution dropdown

It's easiest to think of this in terms of the entries of the product, rather than forming the third-order tensor $\f_{y,y}$ directly.
The $j$-th entry of $\bflambda\trp\f_{y,y}\v_y$ is

$$
\begin{aligned}
    \sum_{i=1}^{n_y}\sum_{k=1}^{n_y}
    \lambda_{i}\frac{\partial^{2} f_{i}}{\partial y_{j}\partial y_{k}} v_k.
\end{aligned}
$$

Start by listing all $(i,j,k)$ tuples for which $\frac{\partial^{2} f_{i}(\y,\z,t)}{\partial y_{j}\partial y_{k}}$ is nonzero:

<!-- BUG: enclose with ```{div} ``` in future jupyterlab-myst version to shrink dropdown.
https://github.com/executablebooks/jupyter-book/issues/1928
-->

| $(i,j,k)$ | $\frac{\partial^{2}}{\partial y_{j}\partial y_{k}}f_{i}(\y,\z,t)$              |
| :-------: | :----------------------------------------------------------------------------- |
| $(2,1,1)$ | $\partial_{y_{1}}\partial_{y_{1}}\left[y_1y_4^2 - k/y_1^2\right] = -6k/y_1^4$  |
| $(2,1,4)$ | $\partial_{y_{1}}\partial_{y_{4}}\left[y_1y_4^2 - k/y_1^2\right] = 2y_4$       |
| $(2,4,1)$ | $\partial_{y_{4}}\partial_{y_{1}}\left[y_1y_4^2 - k/y_1^2\right] = 2y_4$       |
| $(2,4,4)$ | $\partial_{y_{4}}\partial_{y_{4}}\left[y_1y_4^2 - k/y_1^2\right] = 2y_1$       |
| $(4,1,1)$ | $\partial_{y_{1}}\partial_{y_{1}}\left[- 2 y_2 y_4 / y_1\right]  = -4y_2y_4/y_1^3$ |
| $(4,1,2)$ | $\partial_{y_{1}}\partial_{y_{2}}\left[- 2 y_2 y_4 / y_1\right]  = 2y_4/y_1^2$ |
| $(4,1,4)$ | $\partial_{y_{1}}\partial_{y_{4}}\left[- 2 y_2 y_4 / y_1\right]  = 2y_2/y_1^2$ |
| $(4,2,1)$ | $\partial_{y_{2}}\partial_{y_{1}}\left[- 2 y_2 y_4 / y_1\right]  = 2y_4/y_1^2$ |
| $(4,4,1)$ | $\partial_{y_{4}}\partial_{y_{1}}\left[- 2 y_2 y_4 / y_1\right]  = 2y_2/y_1^2$ |

Now gather the terms with like $j$ and multiply with the appropriate $\lambda_i$ and $v_k$'s:

$$
\begin{aligned}
    \bflambda\trp\f_{y,y}(\y(t),\z,t)\v_y
    = \left(\begin{array}{c}
        -6k\lambda_{2}v_{1}/y_1^4
        + 2\lambda_{2}y_{4}v_{4}
        - 4\lambda_{4}y_{2}y_{4}v_{1}/y_1^3
        + 2\lambda_{4}y_{4}v_{2}/y_1^2
        + 2\lambda_{4}y_{2}v_{4}/y_1^2
        \\
        2\lambda_{4}y_4v_{1}/y_1^2
        \\ 0 \\
        2\lambda_{2}y_4v_{1}
        + 2\lambda_{2}y_1v_{4}
        + 2\lambda_{4}y_2v_{1}/y_1^2
    \end{array}\right).
\end{aligned}
$$

:::
::::

### Implementation Verification

### Solving the Problem

### Automatic Differentiation

## Example 2: Trajectory Control Problem

:::{warning}
The rest of this page is under construction, please check back later.
:::
