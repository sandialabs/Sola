# 2: Optimization with Differential Constraints

This tutorial considers two [constrained optimization problems](../../problems/optimization) of the form

$$
\begin{aligned}
    \min_{\y(t),\z} ~& \int_{0}^{T} g(\y(t),t) dt + R(\z)
    \\
    s.t. ~~& \frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t), \z, t), \quad \y(0) = \h(\z),
\end{aligned}
$$

where the time-dependent state $\y(t)$ is the solution to an ordinary differential equation (ODE).
In the first problem, the control vector $\z$ specifies an initial condition for the state dynamics; in the second, $\z$ gathers the values of a time-dependent function that alters the state evolution at each time step.

In each case, we will implement subclasses of {class}`Dynamic_Objective` and {class}`Dynamic_Constraint` by explicitly calculating the derivatives of $g$, $R$, $\f$, and $\h$.
<!-- Afterward, we use {class}`Dynamic_Objective_AD` and {class}`Dynamic_Constraint_AD` to use automatic differentiation to calculate the derivatives. -->

:::{note}
This tutorial includes a few short exercises and their solutions.
The finished produced is included in the SABL souce code under `tutorials/Tutorial_2/`.
:::

## Example 1: Selecting Initial Conditions

This example sets up and solves a constrained optimization with a low-dimensional control vector.

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
    s.t. ~~& \frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t)), \quad \y(0) = (z_1, 0, 0, z_2),
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
    \frac{\partial}{\partial y_j}g(\y(t), t)
    = \frac{\partial}{\partial y_j}\left[\sum_{k=1}^{n_y}(y_k(t) - \alpha_k(t))^2\right]
    = 2(y_j - \alpha_j(t)).
\end{aligned}
$$

Hence, $\grad{y}g(\y(t), t) = 2(\y(t) - \boldsymbol{\alpha}(t))$.

:::{admonition} Exercise
:class: exercise

Calculate the Hessian action $\grad{y,y}g(\y(t), t)\v_y$ where $\v_y\in\R^{n_y}$.

:::

:::{admonition} Solution
:class: solution dropdown

First, calculate the second derivatives of $g$ with respect to the entries of $\y(t)$:

$$
\begin{aligned}
    \frac{\partial}{\partial y_i\partial y_j}g(\y(t), t)
    = \frac{\partial}{\partial y_i}\left[2(y_j(t) - \alpha_j(t))\right]
    = 2\delta_{ij}
    = \begin{cases}
    2 & \textrm{if}~i = j, \\
    0 & \textrm{else}.
    \end{cases}
\end{aligned}
$$

This tells us that $\grad{y,y}g(\y(t), t) = 2\I$, where $\I\in\R^{n_y\times n_y}$ is the identity matrix.
Therefore, for all $t$,

$$
\begin{aligned}
    \grad{y,y}g(\y(t), t)\v_y
    = 2\v_y.
\end{aligned}
$$
:::

We can now define a subclass of {class}`Dynamic_Objective` starting from its [inheritance template](optimization.Dynamic_Objective.template).
Our new class, `Tutorial_2_Objective`, has a constructor for setting the target radius $\rho$ and angular velocity $\omega$.
This new class must be defined in a file named `Tutorial_2_Objective.m`.
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
            this = this@Dynamic_Objective(4, 2, T, n_t);   % n_y = 4, n_z = 2.
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

:::{admonition} Exercise
:class: exercise

Implement, in a vectorized fashion, `g_yy_Apply()` for computing $\grad{y,y}g(\y(t), t)\v_y$.

:::

:::{admonition} Solution
:class: solution dropdown

```matlab
function [y_out] = g_yy_Apply(this, y_in, y, t)
    y_out = 2 * y_in;
end
```

:::

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
    \f_y(\y(t), \z, t)
    = \left(\begin{array}{cccc}
        0 & 1 & 0 & 0 \\
        y_4(t)^2 - \frac{2k}{y_1(t)^3} & 0 & 0 & 2y_1(t) y_4(t) \\
        0 & 0 & 0 & 1 \\
        \frac{y_2(t) y_4(t)}{y_1(t)^2} & -\frac{y_4(t)}{y_1(t)} & 0 & -\frac{y_2(t)}{y_1(t)} \\
    \end{array}\right).
\end{aligned}
$$

:::{admonition} Exercise
:class: exercise

Calculate the Hessian action $\bflambda\trp\f_{y,y}(\y(t),\z,t)\v_y$.

:::

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

Start by listing all $(i,j,k)$ tuples for which $\frac{\partial^{2}}{\partial y_{j}\partial y_{k}}f_{i}(\y(t),\z,t)$ is nonzero:

<!-- BUG: enclose with ```{div} ``` in future jupyterlab-myst version to shrink dropdown.
https://github.com/executablebooks/jupyter-book/issues/1928
-->

| $(i,j,k)$ | $\frac{\partial^{2}}{\partial y_{j}\partial y_{k}}f_{i}(\y(t),\z,t)$              |
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

We now define a new class, `Tutorial_2_Constraint`, that inherits from {class}`Dynamic_Constraint` and implements its abstract methods.
See the [inheritance template](optimization.Constraint.template) to get started.
The class must be defined in a file named `Tutorial_2_Constraint.m`.
Note that in the Jacobian methods, `y_in`, and `z_in` are always treated as matrices, not column vectors.

```matlab
classdef Tutorial_2_Constraint < Dynamic_Constraint

    properties
        k            % Constant of proportionality.
    end

    methods (Access = public)

        function this = Tutorial_2_Constraint(T, n_t, k)
            this = this@Dynamic_Constraint(4, 2, T, n_t);   % n_y = 4, n_z = 2.
            this.k = k;
        end

        function [f, f_y, f_z] = f(this, y, ~, ~)

            % Value of f.
            f = [y(2)
                 y(1) * y(4)^2 - (this.k / y(1)^2)
                 y(4)
                 -2 * y(2) * y(4) / y(1)];

            % y-Jacobian of f.
            f_y = [0, 1, 0, 0
                   y(4)^2 + 2 * (this.k / y(1)^3), 0, 0, 2 * y(1) * y(4)
                   0, 0, 0, 1
                   2 * y(2) * y(4) / y(1)^2, (-2 * y(4)) / y(1), 0, (-2 * y(2)) / y(1)];

            % z-Jacobian of f.
            f_z = zeros(this.n_y, this.n_z);
        end

        function [h, h_z] = h(~, z)
            % Value of h.
            h = [z(1); 0; 0; z(2)];

            % z-Jacobian of h.
            h_z = [1 0
                   0 0
                   0 0
                   0 1];
        end

        function [y_out] = f_yy_Apply(this, y_in, y, z, t, lambda)
            y_out = error('f_yy_Apply() not implemented');
        end

        % These Hessian actions are all zero.
        function [y_out] = f_yz_Apply(this, z_in, ~, ~, ~, ~)
            y_out = zeros(this.n_y, size(z_in, 2));
        end

        function [z_out] = f_zy_Apply(this, y_in, ~, ~, ~, ~)
            z_out = zeros(this.n_z, size(y_in, 2));
        end

        function [z_out] = f_zz_Apply(this, z_in, ~, ~, ~, ~)
            z_out = zeros(this.n_z, size(z_in, 2));
        end

        function [z_out] = h_zz_Apply(this, z_in, ~, ~)
            z_out = zeros(this.n_z, size(z_in, 2));
        end

    end
end
```

:::{admonition} Exercise
:class: exercise

Implement, in a vectorized fashion, `f_yy_Apply()` for computing $\bflambda\trp\f_{y,y}(\y(t), \z, t)\v_y$.

:::

:::{admonition} Solution
:class: solution dropdown

This potential implementation defines each term listed in the solution to Exercise 3.
Note that `y_in` is treated as a matrix.

```matlab
function [y_out] = f_yy_Apply(this, y_in, y, ~, ~, lambda)
    ijk211 = -6 * this.k * lambda(2) * y_in(1, :) / y(1)^4;
    ijk214 = 2 * lambda(2) * y(4) * y_in(4, :);
    ijk241 = 2 * lambda(2) * y(4) * y_in(1, :);
    ijk244 = 2 * lambda(2) * y(1) * y_in(4, :);
    ijk411 = -4 * lambda(4) * y(2) * y(4) * y_in(1, :) / y(1)^3;
    ijk412 = 2 * lambda(4) * y(4) * y_in(2, :) / y(1)^2;
    ijk414 = 2 * lambda(4) * y(2) * y_in(4, :) / y(1)^2;
    ijk421 = 2 * lambda(4) * y(4) * y_in(1, :) / y(1)^2;
    ijk441 = 2 * lambda(4) * y(2) * y_in(1, :) / y(1)^2;
    y_out = [ijk211 + ijk214 + ijk411 + ijk412 + ijk414
             ijk421
             zeros(1, size(y_in, 2))
             ijk241 + ijk244 + ijk441];
end
```

:::

<!-- TODO: ### Implementation Verification -->

### Solving the Problem

With all derivatives implemented, we can now form and solve the optimization problem.
The following script, `Tutorial_2.m`, instantiates `Tutorial_2_Objective` and `Tutorial_2_Constraint` objects and couples them with a {class}`Reduced_Space_Optimization` object.

```matlab
% Tutorial_2.m

%% Clear workspace and add the SABL source path.
clear;
close all;
clc;
run('~/Software/SABL/src/Set_Paths');   % MODIFY THIS TO MATCH YOUR PATH TO SABL.
rng(18);                                % Random seed for reproducing results (optional).

%% Instantiate the objective and constraints.
T = 4 * pi;
n_t = 200;

radius = 1;
velocity = 1;
proportionality = radius^3 * velocity^2;

objective = Tutorial_2_Objective(T, n_t, radius, velocity);
constraint = Tutorial_2_Constraint(T, n_t, proportionality);

%% Solve the optimization problem.
optimizer = Reduced_Space_Optimization(objective, constraint);
z0 = 1 + randn(objective.n_z, 1) / 3;   % Initial guess for the control.
[u, z] = optimizer.Optimize(z0);

disp("Control:");
disp(z);

objective.Plot(u);
```

```text
                                Norm of      First-order
 Iteration        f(x)          step          optimality   CG-iterations
     0            116.822                      3.49e+03
     1            116.822             10       3.49e+03           1
     2            116.822            2.5       3.49e+03           0
     3            116.822          0.625       3.49e+03           0
     4            116.822        0.15625       3.49e+03           0
     5            17.8439      0.0390625       1.21e+03           0
     6            5.77142      0.0390625            698           2
     7           0.190283       0.012319            131           2
     8         0.00325773     0.00934961           7.89           2
     9        0.000110624      0.0154862           2.73           2
    10        8.31968e-08       0.002003         0.0443           2
    11        1.12288e-13    8.76758e-05       8.75e-05           2

Optimization stopped because the relative objective function value is changing
by less than options.FunctionTolerance = 1.000000e-06.

Control:
    1.0000
    1.0000
```

The optimal solution to this problem, unsurprisingly, is to start with an initial radius of $\rho$ and angular velocity $\omega$, which achieves an objective function value of zero.
The following figure plots the state trajectory `u`.

:::{figure} ../../../img/t2-fig1.svg
:align: center
:width: 100 %

Figure 1: Results of `Tutorial_1.m`.
:::

:::{admonition} Plotting code
:class: tip dropdown

This is a method of the `Tutorial_2_Objective` class.

```matlab
function Plot(this, u)
    % Plot a trajectory two ways: x(t) and y(t) in time, and
    % (x(t),y(t)) as a 2D trajectory.
    %
    % Parameters
    % ----------
    % u
    %   Solution array, either an :math:`n_y \times n_t` matrix or
    %   an :math:`(n_y n_t, 1)` vector.

    t = this.t_mesh;
    x_goal = this.radius * cos(this.velocity .* t);
    y_goal = this.radius * sin(this.velocity .* t);

    y = reshape(u, this.n_y, this.n_t);
    r = y(1, :);
    theta = y(3, :);
    x1 = r .* cos(theta);
    x2 = r .* sin(theta);
    lim = 1.2 * max(r);

    % Plot the state trajectory coordinates in time.
    fig = figure;
    fig.Position(3:4) = [830, 300];
    subplot(1, 2, 1);
    plot(t, x1, '-', 'LineWidth', 2);
    hold on;
    plot(t, x2, '-', 'LineWidth', 2);
    plot(t, -ones(this.n_t, 1), 'k-', 'LineWidth', 0.1);
    plot(t, ones(this.n_t, 1), 'k-', 'LineWidth', 0.1);
    xlim([0, t(end)]);
    ylim([-lim, lim]);
    xlabel('$$t$$', 'Interpreter', 'latex');
    title('Position coordinates');
    legend({'$$x_1(t)$$', '$$x_2(t)$$', '', ''}, ...
           'Location', 'southeast', 'Interpreter', 'latex');

    % Plot the state trajectory in two-dimensional space.
    subplot(1, 2, 2);
    plot(x_goal, y_goal, 'k--', 'LineWidth', 2);
    hold on;
    plot(x1, x2, '-', 'LineWidth', 1);
    plot(x1(1), x2(1), '.', 'MarkerSize', 16);
    xlim([-lim, lim]);
    ylim([-lim, lim]);
    xlabel('$$x_1(t)$$', 'Interpreter', 'latex');
    ylabel('$$x_2(t)$$', 'Interpreter', 'latex');
    title('Position');
    axis('equal');
    legend({'target trajectory', 'realized trajectory', 'initial position'}, ...
           'Location', 'southeast', 'Interpreter', 'latex');
end
```

:::

## Example 2: Trajectory Control Problem

This example shows how to deal with controls that represent time-dependent functions.

### New Problem Statement

Consider now a version of the satellite problem in which the initial conditions are fixed and the satellite has thrusters that can be applied in the radial and tangential directions, represented by $q_1(t)$ and $q_2(t)$, respectively.
The new system dynamics are

$$
\begin{aligned}
    \frac{\textrm{d}}{\textrm{d}t}\y(t)
    = \f(\y(t),\z,t) := \left(\begin{array}{c}
        y_2(t) \\
        y_1(t)y_4(t)^2 - k/y_1(t)^2 + q_1(t)\\
        y_4(t) \\
        (q_2(t) - 2 y_2(t) y_4(t)) / y_1(t)
    \end{array}\right).
\end{aligned}
$$

The control $\z$ collects the values of $q_1(t)$ and $q_2(t)$ at all time mesh points where the state is recorded **except** for the initial time (this is because the time stepper uses a first-order backward Euler scheme to solve the ODE), i.e., $t_2,\ldots,t_{n_t}$.
Hence, the control dimension is $n_z = 2(n_t - 1)$ where $n_t$ is the number of points in the temporal discretization.
In this example, we order the entries of the control as

$$
\begin{aligned}
    \z = \left(\begin{array}{c}
        q_1(t_2) \\ q_2(t_2) \\ q_1(t_3) \\ q_2(t_3) \\ \vdots \\ q_1(t_{n_t}) \\ q_2(t_{n_t})
    \end{array}\right)
    \in \R^{n_z}.
\end{aligned}
$$

Our objective is to achieve an orbit of radius $\rho$ with constant angular velocity $\omega$ by guiding the satellite via the thrusters.
Since thruster fuel is expensive, we add a term to penalize the norm of the controls, resulting in the following constrained optimization:

$$
\begin{aligned}
    \min_{\y(t),\z} ~& \int_{0}^{T} \|\y(t) - \boldsymbol{\alpha}(t)\|_2^2 dt + \gamma \int_{t_2}^{T}\|\q(t)\|_2^2 dt
    \\
    s.t. ~~& \frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t), \z,t), \quad \y(0) = (\rho_0, 0, 0, \omega_0),
\end{aligned}
$$

where $\boldsymbol{\alpha}(t) = (\rho, 0, \omega t, \omega)\trp$ as before , $\rho_0 > 0$ is the initial position of the satellite (radius), $\omega_0 \in \R$ is the initial angular velocity of the satellite, $\gamma \ge 0$ is a regularization parameter, and $\q(t) = (q_1(t),q_2(t))\trp$ is the control at time $t$.
<!-- The new regularization term effectively penalizes the amount of fuel used for the satellite thrusters: a larger $\gamma$ should result in conservative control policies (small $\|\q(t)\|_2$, gentle adjustments) while a smaller $\gamma$ allows for more aggressive control policies (large $\|\q(t)\|_2$, strong adjustments). -->

### Implementation

The new objective function has the same $\g(\y(t),t)$ as in the previous example, but now we have

$$
\begin{aligned}
    R(\z)
    = \gamma \int_{t_2}^{T}\|\q(t)\|_2^2 \,dt
    = \gamma \int_{t_2}^{T} q_1(t)^2 + q_2(t)^2 \,dt.
\end{aligned}
$$

We need to represent this function with the discrete values recorded in $\z$.
Using the trapezoidal rule to estimate the integral, we have

$$
\begin{aligned}
    R(\z)
    \approx \gamma \sum_{j=2}^{n_t} w_{j} (q_1(t_j)^2 + q_2(t_j)^2),
\end{aligned}
$$

where $w_2 = w_{n_t} = \frac{1}{2}\delta t$ and $w_j = \delta t$, $j = 3,\ldots,n_t-1$, with temporal spacing $\delta t = (T - t_{2}) / (n_t - 2)$.
We can therefore write

$$
\begin{aligned}
    R(\z) \approx \gamma \w\trp(\z \ast \z),
    \qquad
    % \w = (w_2,w_2,w_3,w_3,\ldots,w_{n_t},w_{n_t})\trp.
    \w = \left(\begin{array}{c}
        w_2 \\ w_2 \\ w_3 \\ w_3 \\ \vdots \\ w_{n_t} \\ w_{n_t}
    \end{array}\right),
\end{aligned}
$$

where $\ast$ is the Hadamard (element-wise) product.

The following subclass of {class}`Dynamic_Objective` records the target radius $\rho$, target velocity $\omega$, and regularization parameter $\gamma$ in the constructor, then forms the weight vector $\w$ defined above.
The implementation of `g()` and `g_yy_Apply()` are the same as in `Tutorial_2_Objective`, but `R()` and `R_zz_Apply()` are different.

```matlab
classdef Tutorial_2B_Objective < Dynamic_Objective

    properties
        radius       % Target orbital radius.
        velocity     % Target angular velocity.
        regularizer  % Regularization hyperparameter.
        n_q          % Number of controls at each fixed time.
        w_z          % Quadrature weights for the control regularization.
    end

    methods

        function this = Tutorial_2B_Objective(T, n_t, radius, velocity, regularizer)
            n_q = 2;
            n_z = n_q * (n_t - 1);
            this = this@Dynamic_Objective(4, n_z, T, n_t);
            this.n_q = n_q;
            this.radius = radius;
            this.velocity = velocity;
            this.regularizer = regularizer;

            % Quadrature weights for the control regularization.
            w = ones(n_t - 1, 1);
            w(1) = 0.5;
            w(2) = 0.5;
            w = (T - this.t_mesh(2)) * w / (n_t - 2);
            this.w_z = repelem(w, this.n_q);
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
            val = this.regularizer * this.w_z' * (z.^2);
            grad_z = error('grad_z not implemented');
        end

        function [y_out] = g_yy_Apply(~, y_in, ~, ~)
            y_out = 2 * y_in;
        end

        function [z_out] = R_zz_Apply(this, z_in, z)
            z_out = error('R_zz_Apply() not implemented');
        end

    end
end
```

:::{admonition} Exercise
:class: exercise

Calculate the gradient $\grad{z}R(\z)$ and the Hessian action $\grad{z,z}R(\z)\v_z$.
Then, finish implementing `R()` and `R_zz_Apply()`.

:::

:::{admonition} Solution
:class: solution dropdown

Let us rewrite the entries of the weight vector as $\w = (\tilde{w}_1,\ldots,\tilde{w}_{n_z})\trp$.
Then we have

$$
\begin{aligned}
    R(\z) = \gamma\sum_{k=1}^{n_z}\tilde{w}_k z_k^2.
\end{aligned}
$$

The first and second derivatives of $R$ are then

$$
\begin{aligned}
    \frac{\partial R(\z)}{\partial z_j}
    = \frac{\partial}{\partial z_j}\left[\gamma\sum_{k=1}^{n_z}\tilde{w}_k z_k^2\right]
    = 2 \gamma \tilde{w}_{j} z_{j},
    \\
    \frac{\partial^{2} R(\z)}{\partial z_i \partial z_j}
    = \frac{\partial}{\partial z_i}\left[2 \gamma \tilde{w}_{j} z_{j}\right]
    = 2 \gamma \tilde{w}_{i} \delta_{ij}.
\end{aligned}
$$

Hence, $\grad{z}R(\z) = 2\gamma(\w\ast\z)$ and $\grad{z,z}R(\z) = 2\gamma\operatorname{diag}(\w)$, so $\grad{z,z}R(\z)\v_z = 2\gamma (\w\ast\v_z)$.

```matlab
function [val, grad_z] = R(this, z)
    val = this.regularizer * this.w_z' * (z.^2);
    grad_z = 2 * this.regularizer * this.w_z .* z;
end

function [z_out] = R_zz_Apply(this, z_in, ~)
    z_out = 2 * this.regularizer * this.w_z .* z_in;
end
```

:::

To simplify the implementation of the constraint, we use the Gauss--Newton Hessian approximation specified in {prf:ref}`alg:adjoint_gaussnewton`, so we do not need to implement the Hessian actions of $\f$ or $\h$.
Because $\h$ is a constant function, its $\z$ Jacobian is the $n_y \times n_z$ zero matrix.
The tricky part of implementing `f()` is extracting $(q_1(t),q_2(t))$ from $\z$ at the given time $t$.
The following subclass of {class}`Dynamic_Constraint` handles this with a private method, `Input_Indices()`.

```matlab
classdef Tutorial_2B_Constraint < Dynamic_Constraint

    properties
        k            % Constant of proportionality.
        n_q          % Number of controls at each fixed time.
        y0           % Initial condition.
    end

    methods (Access = public)

        function this = Tutorial_2B_Constraint(T, n_t, r0, w0, k)
            n_q = 2;
            n_z = n_q * (n_t - 1);
            this = this@Dynamic_Constraint(4, n_z, T, n_t);
            this.y0 = [r0; 0; 0; w0];
            this.n_q = n_q;
            this.k = k;
        end

        function [f, f_y, f_z] = f(this, y, z, t)
            I = this.Input_Indices(t);
            q = z(I);

            f = [y(2)
                 y(1) * y(4)^2 - (this.k / y(1)^2) + q(1)
                 y(4)
                 (q(2) - 2 * y(2) * y(4)) / y(1)];

            f_y = [0, 1, 0, 0
                   y(4)^2 + 2 * (this.k / y(1)^3), 0, 0, 2 * y(1) * y(4)
                   0, 0, 0, 1
                   (2 * y(2) * y(4) - q(2)) / y(1)^2, (-2 * y(4)) / y(1), 0, (-2 * y(2)) / y(1)];

            f_z = error('f_z not implemented');
        end

        function [h, h_z] = h(this, ~)
            h = this.y0;
            h_z = zeros(this.n_y, this.n_z);
        end

    end

    methods (Access = private)

        function [mask] = Input_Indices(this, t)
            % Get the indices of the control at time t, i.e.,
            % I = Input_Indices(t) --> q(t) = z(I).

            [~, t_index] = min(abs(t - this.t_mesh));
            if t_index == 1
                error('no control at initial time!');
            end
            idx = t_index - 1;
            mask = (this.n_q * (idx - 1) + 1):(this.n_q * idx);
        end

    end
end
```

:::{admonition} Exercise
:class: exercise

Calculate the Jacobian $\f_\z(\y(t), \z, t)$ and finish implementing `f()`.

:::

:::{admonition} Solution
:class: solution dropdown

The $\z$ Jacobian of $\f$ is the $n_y \times n_z$ matrix

$$
\begin{aligned}
    \f_\z(\y(t), \z, t)
    = \left(\begin{array}{cc|cc|c|cc}
        \frac{\partial f_1}{\partial q_1(t_2)} & \frac{\partial f_1}{\partial q_2(t_2)} &
        \frac{\partial f_1}{\partial q_1(t_3)} & \frac{\partial f_1}{\partial q_2(t_3)} &
        \cdots &
        \frac{\partial f_1}{\partial q_1(t_{n_t})} & \frac{\partial f_1}{\partial q_2(t_{n_t})}
        \\
        \vdots & \vdots & \vdots & \vdots & & \vdots & \vdots
        \\
        \frac{\partial f_{n_y}}{\partial q_1(t_2)} & \frac{\partial f_{n_y}}{\partial q_2(t_2)} &
        \frac{\partial f_{n_y}}{\partial q_1(t_3)} & \frac{\partial f_{n_y}}{\partial q_2(t_3)} &
        \cdots &
        \frac{\partial f_{n_y}}{\partial q_1(t_{n_t})} & \frac{\partial f_{n_y}}{\partial q_2(t_{n_t})}
    \end{array}\right).
\end{aligned}
$$

For $t = t_j$, every entry of this matrix is necessarily zero except for the block

$$
\begin{aligned}
    \left(\begin{array}{cc}
        \frac{\partial f_1}{\partial q_1(t_j)} & \frac{\partial f_1}{\partial q_2(t_j)}
        \\ \vdots & \vdots \\
        \frac{\partial f_{n_y}}{\partial q_1(t_j)} & \frac{\partial f_{n_y}}{\partial q_2(t_j)}
    \end{array}\right)
    = \left(\begin{array}{cc}
        0 & 0 \\
        1 & 0 \\
        0 & 0 \\
        0 & \frac{1}{y_1(t)}
    \end{array}\right).
\end{aligned}
$$

The `Input_Indices()` method provides the indices of the columns of this nonzero block.

```matlab
function [f, f_y, f_z] = f(this, y, z, t)
    I = this.Input_Indices(t);
    q = z(I);

    % ...

    f_z = zeros(this.n_y, this.n_z);
    f_z(:, I) = [0, 0
                 1, 0
                 0, 0
                 0, 1 / y(1)];
end
```

:::

### Optimization Solution

We can now solve this problem with a {class}`Reduced_Space_Optimization` object.
See `tutorials/Tutorial_2/Tutorial_2B_Objective.m` in the source code for the `Plot()` method.

```matlab
% Tutorial_2B.m

%% Clear workspace and add the SABL source path.
clear;
close all;
clc;
run('~/Software/SABL/src/Set_Paths');
rng(1342);                              % Random seed for reproducing results (optional).

%% Instantiate the objective and constraints.
T = 4 * pi;                             % Final simulation time.
n_t = 200;                              % Number of time steps.

radius = 1;
velocity = 1;
proportionality = radius^3 * velocity^2;

regularizer = 5;                        % Control regularization parameter.
r0 = 1.5;                               % Initial position (radius).
w0 = 2;                                 % Initial angular velocity.

objective = Tutorial_2B_Objective(T, n_t, radius, velocity, regularizer);
constraint = Tutorial_2B_Constraint(T, n_t, r0, w0, proportionality);

%% Solve the optimization problem with the Gauss-Newton Hessian approximation.
optimizer = Reduced_Space_Optimization(objective, constraint);
optimizer.Gauss_Newton_Hess = true;     % Use Algorithm 3, not Algorithm 2.

z0 = randn(objective.n_z, 1);           % Initial guess for the control.
[u, z] = optimizer.Optimize(z0);

%% Plot the results.
objective.Plot(u, z);
```

```text
                                Norm of      First-order
 Iteration        f(x)          step          optimality   CG-iterations
     0            6568.12                           279
     1            1229.31             10           49.5           8
     2            520.325        17.7446             15           7
     3            113.746        9.29213           9.62           9
     4            113.746        7.25114           9.62          13
     5            70.2078        1.81279           3.99           0
     6            43.6528        3.62557           6.36          14
     7            43.6528        3.03051           6.36          12
     8            36.4659       0.757627           2.27           0
     9            32.5519        1.51525           5.85          15
    10            29.8283        1.41729           1.77          13
    11            29.4634       0.434592            1.4          12
    12            29.4026      0.0778727        0.00827          11
    13            29.4025       0.025252        0.00329          18
    14            29.4025       0.012564        0.00165          19

Optimization stopped because the relative objective function value is changing
by less than options.FunctionTolerance = 1.000000e-06.
```

:::{figure} ../../../img/t2-fig2.svg
:align: center
:width: 100 %

Figure 2: Results of `Tutorial_2.m`.
:::

:::{admonition} Plotting code
:class: tip dropdown

This is a method of the `Tutorial_2B_Objective` class.

```matlab
function Plot(this, u, z)
    % Plot a trajectory two ways: x(t) and y(t) in time, and
    % (x(t),y(t)) as a 2D trajectory.
    %
    % Parameters
    % ----------
    % u
    %   Solution array, either an :math:`n_y \times n_t` matrix or
    %   an :math:`(n_y n_t, 1)` vector.

    t = this.t_mesh;
    x_goal = this.radius * cos(this.velocity .* t);
    y_goal = this.radius * sin(this.velocity .* t);

    y = reshape(u, this.n_y, this.n_t);
    r = y(1, :);
    theta = y(3, :);
    x1 = r .* cos(theta);
    x2 = r .* sin(theta);
    lim = 1.1 * max(r);
    t_q = t(2:end);
    q = reshape(z, this.n_q, []);

    % Plot the state trajectory coordinates in time.
    fig = figure;
    fig.Position(3:4) = [830, 535];
    subplot(2, 2, 1);
    plot(t, x1, '-', 'LineWidth', 2);
    hold on;
    plot(t, x2, '-', 'LineWidth', 2);
    plot(t, -ones(this.n_t, 1), 'k-', 'LineWidth', 0.1);
    plot(t, ones(this.n_t, 1), 'k-', 'LineWidth', 0.1);
    xlim([0, t(end)]);
    ylim([-lim, lim]);
    xlabel('$$t$$', 'Interpreter', 'latex');
    title('Position coordinates');
    legend({'$$x_1(t)$$', '$$x_2(t)$$', '', ''}, ...
           'Location', 'southeast', 'Interpreter', 'latex');

    % Plot the state trajectory in two-dimensional space.
    subplot(2, 2, 2);
    plot(x_goal, y_goal, 'k--', 'LineWidth', 2);
    hold on;
    plot(x1, x2, '-', 'LineWidth', 1);
    plot(x1(1), x2(1), '.', 'MarkerSize', 16);
    xlim([-lim, lim]);
    ylim([-lim, lim]);
    xlabel('$$x_1(t)$$', 'Interpreter', 'latex');
    ylabel('$$x_2(t)$$', 'Interpreter', 'latex');
    title('Position');
    axis('equal');
    legend({'target trajectory', 'realized trajectory', 'initial position'}, ...
           'Location', 'southeast', 'Interpreter', 'latex');

    % Plot the optimal control scheme as a function of time.
    subplot(2, 1, 2);
    plot(t_q, q(1, :), 'LineWidth', 2);
    hold on;
    plot(t_q, q(2, :), 'LineWidth', 2);
    xlim([0, t(end)]);
    xlabel('$$t$$', 'Interpreter', 'latex');
    title('Control profile');
    legend({'$$q_1(t)$$', '$$q_2(t)$$'}, ...
           'Location', 'southeast', 'Interpreter', 'latex');
end
```

:::

:::{admonition} Exercise
:class: exercise

Experiment with different values of $\gamma$ by changing `regularizer` in `Tutorial_2B.m`.
How do the controls and the resulting state trajectory differ for smaller $\gamma$ versus larger $\gamma$?
:::

:::{tip}
The full control vector $\z$ measures $q_1(t)$ and $q_2(t)$ at the same discrete time points (except for the initial time) as the state, the `t_mesh` property which is initialized in the constructor of the {class}`Dynamic_Objective` class.
However, it is also possible to have a distinct time domain for the controls and the state (for example, perhaps the satellite can only receive new commands from mission control every hour).
Future examples incorporate this approach.
:::
