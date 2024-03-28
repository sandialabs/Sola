# 1: Optimization with Algebraic Constraints

This example considers a simple [constrained optimization problem](../../problems/optimization) of the form

$$
\begin{aligned}
    \min_{\u,\z} ~& J(\u,\z)
    \\
    s.t. ~~& c(\u,\z) = 0
\end{aligned}
$$

where the state $\u$ and control $\z$ are a finite-dimensional vectors and the constraint $\c(\u,\z)$ represents a nonlinear system of algebraic equations.

We will implement subclasses of {class}`Objective` and {class}`Constraint` by explicitly calculating the derivatives of $J$ and $\c$, then show how to solve the optimization problem with a {class}`Reduced_Space_Optimization`.
Afterward, we use {class}`Objective_AD` and {class}`Constraint_AD` to use automatic differentiation to calculate the derivatives.

:::{note}
This tutorial includes a few short exercises and their solutions.
The finished product is included in the SABL source code under `tutorials/Tutorial_1/`.
:::

## Problem Statement

Let $n_{u} = 3$ be the state dimension, $n_{z} = 2$ be the control dimension, and denote the state and control by $\u = (u_{1},u_{2},u_{3})\trp$ and $\z = (z_{1},z_{2})\trp$, respectively.
We consider the objective function

$$
\begin{aligned}
    J(\u, \z)
    &= (u_{1} - \alpha_{1})^{2} + (u_{2} - \alpha_{2})^{2} + (u_{3} - \alpha_{3})^{2} \\
    &\quad+ (z_{1} - \alpha_{4})^{2} + (z_{2} - \alpha_{5})^{2} + (u_{1}z_{1} - \alpha_{1}\alpha_{4})^{2},
\end{aligned}
$$

where $\alpha_{1},\alpha_{2},\alpha_{3},\alpha_{4},\alpha_{5}\in\mathbb{R}$ are known constants.
The $n_{u} = 3$ constraints are encoded by the function

$$
\begin{aligned}
    \c(\u, \z)
    &= \left(\begin{array}{c}
        u_1 + u_2 - z_1 \\
        z_1u_2 - z_2 \\
        u_3^3 - z_2^2
    \end{array}\right).
\end{aligned}
$$

When the constraint equation $\c(\u, \z) = \0$ is satisfied, we have the equations

$$
\begin{aligned}
    u_1 + u_2 &= z_1,
    &
    z_1u_2 &= z_2,
    &
    u_3^3 &= z_2^2.
\end{aligned}
$$

:::{note}
Without constraints, the choice of $\u$ and $\z$ that minimizes $J(\u,\z)$ is $\u^{*} = (\alpha_{1},\alpha_{2},\alpha_{3})\trp$, $\z^{*} = (\alpha_{4},\alpha_{5})\trp$.
Depending on $\alpha_1,\ldots,\alpha_5$, this solution may not satisfy the constraints $\c(\u,\z) = \0$.
However, for any choice of $\alpha_1$ and $\alpha_2$, setting

$$
\begin{aligned}
    \alpha_4 &= \alpha_1 + \alpha_2,
    &
    \alpha_5 &= \alpha_2 \alpha_4,
    &
    \alpha_3 &= (\alpha_2 + \alpha_4)^{2/3}
\end{aligned}
$$

guarantees that $\c(\u^{*},\z^{*}) = \0$.
We will use this observation to design a sanity check for our implementation.
:::

To minimize the objective $J$ while satisfing the constraints $\c$, we need to represent $J$ and $\c$ _and their derivatives_ with SABL classes.
We will start by calculating and implementing the derivatives by hand, but later we will show how to leverage automatic differentiation.

## Implementing the Objective

### Objective Value and Gradients

To implement the abstract methods of {class}`Objective` for this particular choice of $J$, we need to calculate the $\u$ and $\z$ gradients of $J$.
The $i$th entry of the gradient of $J$ with respect to a vector $\x = (x_1,\ldots,x_n)\trp$ is $\frac{\partial J(\u,\z)}{\partial x_{i}}$.
To begin, we have

$$
\begin{aligned}
    \frac{\partial J(\u,\z)}{\partial u_{1}}
    &= \frac{\partial}{\partial u_{1}}\left[(u_{1} - \alpha_{1})^{2} + (u_{2} - \alpha_{2})^{2} + (u_{3} - \alpha_{3})^{2} + (z_{1} - \alpha_{4})^{2} + (z_{2} - \alpha_{5})^{2} + (u_{1}z_{1} - \alpha_{1}\alpha_{4})^{2}\right]
    \\
    &= \frac{\partial}{\partial u_{1}}\left[(u_{1} - \alpha_{1})^{2} + (u_{1}z_{1} - \alpha_{1}\alpha_{4})^{2}\right]
    \\
    &= 2(u_{1} - \alpha_{1}) + 2(u_{1}z_{1} - \alpha_{1}\alpha_{4})z_{1}.
\end{aligned}
$$

Similarly calculating $\frac{\partial J(\u,\z)}{\partial u_{2}}$ and $\frac{\partial J(\u,\z)}{\partial u_{3}}$, we obtain the $\u$ gradient of $J$,

$$
\begin{aligned}
    \grad{u}J(\u,\z)
    &= \left(\begin{array}{c}
        \frac{\partial J(\u,\z)}{\partial u_{1}} \\
        \frac{\partial J(\u,\z)}{\partial u_{2}} \\
        \frac{\partial J(\u,\z)}{\partial u_{3}}
    \end{array}\right)
    = \left(\begin{array}{c}
        2(u_{1} - \alpha_{1}) + 2(u_{1}z_{1} - \alpha_{1}\alpha_{4})z_{1}
        \\
        2(u_{2} - \alpha_{2})
        \\
        2(u_{3} - \alpha_{3})
    \end{array}\right).
\end{aligned}
$$

::::{admonition} Exercise 1
:class: exercise

Calculate the gradient $\grad{z}J(\u,\z)$.
:::{admonition} Solution
:class: solution dropdown

$$
\begin{aligned}
    \grad{z}J(\u,\z)
    &= \left(\begin{array}{c}
        \frac{\partial J(\u,\z)}{\partial z_{1}} \\
        \frac{\partial J(\u,\z)}{\partial z_{2}}
    \end{array}\right)
    = \left(\begin{array}{c}
        2(z_{1} - \alpha_{4}) + 2(u_{1}z_{1} - \alpha_{1}\alpha_{4})u_{1}
        \\
        2(z_{2} - \alpha_{5})
    \end{array}\right).
\end{aligned}
$$
:::
::::

### Objective Hessians

Next, we need to calculate the action of the Hessians of $J$.
The $(i,j)$th entry of the $(\x,\y)$ Hessian of $J$ is given by $\frac{\partial^{2} J(\u,\z)}{\partial x_{i}\partial y_{j}}$.
Hence, the $(\u,\u)$ Hessian $\grad{u,u}J(\u,\z)$ is a $n_u\times n_u$ matrix.

$$
\begin{aligned}
    \grad{u,u}J(\u,\z)
    &= \left(\begin{array}{ccc}
        \frac{\partial^{2} J(\u,\z)}{\partial u_1^{2}} & \frac{\partial^{2} J(\u,\z)}{\partial u_1\partial u_2} & \frac{\partial^{2} J(\u,\z)}{\partial u_1\partial u_3}
        \\
        \frac{\partial^{2} J(\u,\z)}{\partial u_2\partial u_1} & \frac{\partial^{2} J(\u,\z)}{\partial u_2^{2}} & \frac{\partial^{2} J(\u,\z)}{\partial u_2\partial u_3}
        \\
        \frac{\partial^{2} J(\u,\z)}{\partial u_3\partial u_1} & \frac{\partial^{2} J(\u,\z)}{\partial u_3\partial u_2} & \frac{\partial^{2} J(\u,\z)}{\partial u_3\partial u_3}
    \end{array}\right)
    \\
    &= \left(\begin{array}{ccc}
        \frac{\partial}{\partial u_1}\left[2(u_{1} - \alpha_{1}) + 2(u_{1}z_{1} - \alpha_{1}\alpha_{4})z_{1}\right] & \frac{\partial}{\partial u_1}\left[2(u_{2} - \alpha_{2})\right] & \frac{\partial}{\partial u_1}\left[2(u_{3} - \alpha_{3})\right]
        \\
        \frac{\partial}{\partial u_2}\left[2(u_{1} - \alpha_{1}) + 2(u_{1}z_{1} - \alpha_{1}\alpha_{4})z_{1}\right] & \frac{\partial}{\partial u_2}\left[2(u_{2} - \alpha_{2})\right] & \frac{\partial}{\partial u_2}\left[2(u_{3} - \alpha_{3})\right]
        \\
        \frac{\partial}{\partial u_3}\left[2(u_{1} - \alpha_{1}) + 2(u_{1}z_{1} - \alpha_{1}\alpha_{4})z_{1}\right] & \frac{\partial}{\partial u_3}\left[2(u_{2} - \alpha_{2})\right] & \frac{\partial}{\partial u_3}\left[2(u_{3} - \alpha_{3})\right]
    \end{array}\right)
    \\
    &= \left(\begin{array}{ccc}
        2 + 2z_{1}^{2} & 0 & 0
        \\
        0 & 2 & 0
        \\
        0 & 0 & 2
    \end{array}\right).
\end{aligned}
$$

The _action_ of the Hessian is the matrix-vector product $\grad{u,u}J(\u,\z)\v_u$ where $\v_u \in \R^{n_u}$.
In this case, writing $\v_u = (v_1,v_2,v_3)\trp$, we have

$$
\begin{aligned}
    \grad{u,u}J(\u,\z)\v_u
    &= \left(\begin{array}{ccc}
        2 + 2z_{1}^{2} & 0 & 0
        \\
        0 & 2 & 0
        \\
        0 & 0 & 2
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2} \\ v_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        (2 + 2z_{1}^{2})v_{1}
        \\
        2v_{2}
        \\
        2v_{3}
    \end{array}\right).
\end{aligned}
$$

<!-- :::{tip}
Component notation may be helpful for calculating these products.
Letting $[\![\x]\!]_{i}$ denote the $i$th entry of a vector $\x$ and $[\![\A]\!]_{ij}$ denote the $(i,j)$th entry of a matrix $\A$, the entries of the Hessian-vector product $\grad{u,u}J(\u,\z)\v_u$ are given by

$$
\begin{aligned}
    {[}\![\grad{u,u}J(\u,\z)\v_u]\!]_{i}
    &= \sum_{j=1}^{n_{u}}{[}\![\grad{u,u}J(\u,\z)]\!]_{ij} [\![\v_u]\!]_{j}
    = \sum_{j=1}^{n_{u}}\frac{\partial^{2}}{\partial u_{i}\partial u_{j}}J(\u,\z) v_j.
\end{aligned}
$$

For this problem, there are few nonzero terms in the sum.

$$
\begin{aligned}
    {[}\![\grad{u,u}J(\u,\z)\v_u]\!]_{1}
    &= \frac{\partial^{2}}{\partial u_{1}\partial u_{1}}J(\u,\z) v_{1}
    + \frac{\partial^{2}}{\partial u_{1}\partial u_{2}}J(\u,\z) v_{2}
    + \frac{\partial^{2}}{\partial u_{1}\partial u_{3}}J(\u,\z) v_{3}
\end{aligned}
$$

We will use this strategy for computing the constraint Hessian actions.
::: -->

::::{admonition} Exercise 2
:class: exercise

Calculate the Hessian actions $\grad{u,z}J(\u,\z)\v_u$, $\grad{z,u}J(\u,\z)\v_z$, and $\grad{z,z}J(\z,\z)\v_z$, where $\v_z\in\R^{n_z}$.

:::{admonition} Solution
:class: solution dropdown

First, calculate each Jacobian.

$$
\begin{aligned}
    \grad{u,z}J(\u,\z)
    &= \left(\begin{array}{cc}
        \frac{\partial^{2} J(\u,\z)}{\partial u_1 \partial z_1} & \frac{\partial^{2} J(\u,\z)}{\partial u_1\partial z_2}
        \\
        \frac{\partial^{2} J(\u,\z)}{\partial u_2 \partial z_1} & \frac{\partial^{2} J(\u,\z)}{\partial u_2\partial z_2}
        \\
        \frac{\partial^{2} J(\u,\z)}{\partial u_3 \partial z_1} & \frac{\partial^{2} J(\u,\z)}{\partial u_3\partial z_2}
    \end{array}\right)
    = \left(\begin{array}{cc}
        4u_{1}z_{1} - 2\alpha_{1}\alpha_{4} & 0
        \\ 0 & 0
        \\ 0 & 0
    \end{array}\right)
    \\ \\
    \grad{z,u}J(\u,\z)
    &= \left(\begin{array}{ccc}
        \frac{\partial^{2} J(\u,\z)}{\partial z_1 \partial u_1} & \frac{\partial^{2} J(\u,\z)}{\partial z_1\partial u_2} & \frac{\partial^{2} J(\u,\z)}{\partial z_1\partial u_3}
        \\
        \frac{\partial^{2} J(\u,\z)}{\partial z_2\partial u_1} & \frac{\partial^{2} J(\u,\z)}{\partial z_2 \partial u_2} & \frac{\partial^{2} J(\u,\z)}{\partial z_2\partial u_3}
    \end{array}\right)
    = \left(\begin{array}{ccc}
        4u_{1}z_{1} - 2\alpha_{1}\alpha_{4} & 0 & 0
        \\
        0 & 0 & 0
    \end{array}\right),
    \\ \\
    \grad{z,z}J(\u,\z)
    &= \left(\begin{array}{cc}
        \frac{\partial^{2} J(\u,\z)}{\partial z_1^{2}} & \frac{\partial^{2} J(\u,\z)}{\partial z_1\partial z_2}
        \\
        \frac{\partial^{2} J(\u,\z)}{\partial z_2\partial z_1} & \frac{\partial^{2} J(\u,\z)}{\partial z_2^{2}}
    \end{array}\right)
    = \left(\begin{array}{cc}
        2 + 2u_{1}^{2} & 0
        \\
        0 & 2
    \end{array}\right).
\end{aligned}
$$

Writing $\v_z = (v_1,v_2,v_3)\trp$, we have

$$
\begin{aligned}
    \grad{u,z}J(\u,\z)\v_z
    &= \left(\begin{array}{c}
        (4u_{1}z_{1} - 2\alpha_{1}\alpha_{4})v_{1} \\ 0 \\ 0
    \end{array}\right),
    \\ \\
    \grad{z,u}J(\u,\z)\v_u
    &= \left(\begin{array}{c}
        (4u_{1}z_{1} - 2\alpha_{1}\alpha_{4})v_{1}
        \\
        0
    \end{array}\right),
    \\ \\
    \grad{z,z}J(\u,\z)\v_z
    &= \left(\begin{array}{c}
        (2 + 2u_{1}^{2})v_{1}
        \\ 2v_{2}
    \end{array}\right).
\end{aligned}
$$

:::
::::

### Objective Implementation

We can now define a subclass of {class}`Optimize` starting from its [inheritance template](optimization.Objective.template).
Our new class, `Tutorial_1_Objective`, has a constructor for setting the constants $\alpha_{1},\ldots,\alpha_{5}$, and implements the abstract methods of {class}`Optimize` using the derivatives calculated above.
This new class must be defined in a file named `Tutorial_1_Objective.m`.
In the code, the input vectors $\v_u\in\R^{n_u}$ and $\v_z\in\R^{n_z}$ are called `u_in` and `z_in`, respectively, to remind the user of the vector size.
Similarly, the output vectors `u_out` and `z_out` have $n_u$ and $n_z$ rows, respectively.

```matlab
% Tutorial_1_Objective.m

classdef Tutorial_1_Objective < Objective

    properties
        a   % alpha constants.
    end

    methods (Access = public)

        % Constructor: set the alpha constants.
        function this = Tutorial_1_Objective(a)
            this.a = a;
        end

        function [val, grad_u, grad_z] = J(this, u, z)
            u1z1_minus_a1a4 = (u(1) * z(1)) - (this.a(1) * this.a(4));

            % Compute the value of J(u, z).
            val = sum((u - this.a(1:3)).^2);
            val = val + sum((z - this.a(4:5)).^2);
            val = val + u1z1_minus_a1a4^2;

            % Compute the u gradient of J.
            grad_u = 2 * (u - this.a(1:3));
            grad_u(1) = grad_u(1) + 2 * u1z1_minus_a1a4 * z(1);

            % Compute the z gradient of J.
            grad_z = 2 * (z - this.a(4:5));
            grad_z(1) = grad_z(1) + 2 * u1z1_minus_a1a4 * u(1);
        end

        function [u_out] = J_uu_Apply(this, u_in, u, z)
            u_out = 2 * u_in;
            u_out(1, :) = u_out(1, :) .* (1 + z(1)^2);
        end

        function [u_out] = J_uz_Apply(this, z_in, u, z)
            u_out = zeros(length(u), size(z_in, 2));
            u_out(1, :) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) .* z_in(1, :);
        end

        function [z_out] = J_zu_Apply(this, u_in, u, z)
            error('J_zu_Apply() not implemented');
        end

        function [z_out] = J_zz_Apply(this, u_in, u, z)
            error('J_zz_Apply() not implemented');
        end

    end
end
```

:::{danger}
Note carefully that `J_uu_Apply()` and `J_uz_Apply()` are implemented in a _vectorized_ fashion by treating the input `u_in` or `z_in` as a matrix with possibly more than one column.
If `u_in` and `z_in` were always column vectors, the following implementation would be valid:

```matlab
function [u_out] = J_uu_Apply(this, u_in, u, z)
    u_out = 2 * u_in;
    u_out(1) = u_out(1) * (1 + z(1)^2);
end

function [u_out] = J_uz_Apply(this, z_in, u, z)
    u_out = zeros(length(u), 1);
    u_out(1) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) * z_in(1);
end
```

However, MATLAB's `fminunc()` function may require multiple simultaneous Hessian-vector evaluations, so we cannot assume `u_in` has only one column.
Compare the above code carefully with the correct implementation.

```matlab
function [u_out] = J_uu_Apply(this, u_in, u, z)
    u_out = 2 * u_in;
    u_out(1, :) = u_out(1, :) .* (1 + z(1)^2);
end

function [u_out] = J_uz_Apply(this, z_in, u, z)
    u_out = zeros(length(u), size(z_in, 2));
    u_out(1, :) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) .* z_in(1, :);
end
```

:::

::::{admonition} Exercise 3
:class: exercise

Implement, in a vectorized fashion, `J_zu_Apply()` for computing $\grad{zu}J(\u,\z)\v_u$ and `J_zz_Apply()` for computing $\grad{zz}J(\u,\z)\v_z$.

:::{admonition} Solution
:class: solution dropdown

This is one possible implementation.

```matlab
function [z_out] = J_zu_Apply(this, u_in, u, z)
    z_out = zeros(length(z), size(u_in, 2));
    z_out(1, :) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) .* u_in(1, :);
end

function [z_out] = J_zz_Apply(this, u_in, u, z)
    z_out = 2 * u_in;
    z_out(1, :) = z_out(1, :) .* (1 + u(1)^2);
end
```

:::
::::

## Implementing the Constraint

### Solution Operator

To implement a {class}`Constraint` subclass for this problem, we first need a solution operator $\S:\z\mapsto \u$ such that $\c(\S(\z),\z) = \0$ for all $\z$.
In this case, we can directly solve the constraint equations

$$
\begin{aligned}
    u_1 + u_2 &= z_1,
    &
    z_1u_2 &= z_2,
    &
    u_3^3 &= z_2^2
\end{aligned}
$$

for $\u$.
That is,

$$
\begin{aligned}
    \left(\begin{array}{c}
        u_1 \\ u_2 \\ u_3
    \end{array}\right)
    = \S(\z)
    := \left(\begin{array}{c}
        z_{1} -  (z_{2} / z_{1})\\
        z_{2} / z_{1} \\
        z_{2}^{2/3}
    \end{array}\right).
\end{aligned}
$$

:::{note}
For a more complicated system of equations, we may not be able to write the solution operator $\S$ explicitly.
In fact, evaluating $\S$ often includes complex procedures such as Newton's method or numerical integration schemes.
In outer-loop problems, the cost of evaluating $\S$ is typically the computational bottleneck, especially when $\S$ involves solving a discretized partial differential equation.
:::

### Constraint Jacobians

Next, we calculate the derivatives of $\c$.
We start with $\c_u$ and $\c_z$, which are Jacobian matrices.
Writing

$$
\begin{aligned}
    \c(\u,\z)
    = \left(\begin{array}{c}
        c_{1}(\u,\z) \\ c_{2}(\u,\z) \\ c_{3}(\u,\z)
    \end{array}\right)
    = \left(\begin{array}{c}
        u_1 + u_2 - z_1 \\
        z_1u_2 - z_2 \\
        u_3^3 - z_2^2
    \end{array}\right),
\end{aligned}
$$

the $(i,j)$th entry of the $\u$-Jacobian of $\c$ is $\frac{\partial c_{i}}{\partial u_{j}}$.
Hence, $\c_u(\u,\z)\in\R^{n_u\times n_u}$ and $\c_z(\u,\z)\in\R^{n_u\times n_z}$.

$$
\begin{aligned}
    \c_{u}(\u, \z)
    &= \left(\begin{array}{ccc}
        \frac{\partial c_1(\u,\z)}{\partial u_1} & \frac{\partial c_1(\u,\z)}{\partial u_2} & \frac{\partial c_1(\u,\z)}{\partial u_3} \\
        \frac{\partial c_2(\u,\z)}{\partial u_1} & \frac{\partial c_2(\u,\z)}{\partial u_2} & \frac{\partial c_2(\u,\z)}{\partial u_3} \\
        \frac{\partial c_3(\u,\z)}{\partial u_1} & \frac{\partial c_3(\u,\z)}{\partial u_2} & \frac{\partial c_3(\u,\z)}{\partial u_3}
    \end{array}\right)
    \\
    &= \left(\begin{array}{lll}
        \frac{\partial}{\partial u_1}\left[u_1 + u_2 - z_1\right] & \frac{\partial}{\partial u_2}\left[u_1 + u_2 - z_1\right] & \frac{\partial}{\partial u_3}\left[u_1 + u_2 - z_1\right] \\
        \frac{\partial}{\partial u_1}\left[z_1u_2 - z_2\right] & \frac{\partial}{\partial u_2}\left[z_1u_2 - z_2\right] & \frac{\partial}{\partial u_3}\left[z_1u_2 - z_2\right] \\
        \frac{\partial}{\partial u_1}\left[u_3^3 - z_2^2\right] & \frac{\partial}{\partial u_2}\left[u_3^3 - z_2^2\right] & \frac{\partial}{\partial u_3}\left[u_3^3 - z_2^2\right]
    \end{array}\right)
    = \left(\begin{array}{ccc}
        1 & 1 & 0 \\
        0 & z_{1} & 0 \\
        0 & 0 & 3u_{3}^{2}
    \end{array}\right).
\end{aligned}
$$

::::{admonition} Exercise 4
:class: exercise

Calculate the Jacobian $\c_z(\u,\z)$.

:::{admonition} Solution
:class: solution dropdown

$$
\begin{aligned}
    \c_{z}(\u, \z)
    = \left(\begin{array}{ccc}
        \frac{\partial c_1(\u,\z)}{\partial z_1} & \frac{\partial c_1(\u,\z)}{\partial z_2} \\
        \frac{\partial c_2(\u,\z)}{\partial z_1} & \frac{\partial c_2(\u,\z)}{\partial z_2} \\
        \frac{\partial c_3(\u,\z)}{\partial z_1} & \frac{\partial c_3(\u,\z)}{\partial z_2}
    \end{array}\right)
    &= \left(\begin{array}{cc}
        -1 & 0 \\
        u_{2} & -1 \\
        0 & -2z_{2}
    \end{array}\right).
\end{aligned}
$$

:::
::::

To facilitate {prf:ref}`alg:adjoint_gradient`, we need implementations for solving the linear system $\c_u(\u,\z)\trp\bflambda = \v_u$ (equivalently, computing the inverse-transpose Jacobian-vector product $\c_u(\u,\z)^{-\mathsf{T}}\v_u$) and computing $\c_z(\u,\z)\trp\v_z$, where $\v_u\in\R^{n_u}$ and $\v_z\in\R^{n_z}$ as before.
Furthermore, {prf:ref}`alg:adjoint_hessvec` and {prf:ref}`alg:adjoint_gaussnewton` require computing $\c_z(\u,\z)\v_z$ and solving the linear system $\c_u(\u,\z)\bfmu = \v_u$ (computing $\c_u(\u,\z)^{-1}\v_u$).
Explicitly calculating the inverse of $\c_u$, we have

$$
\begin{aligned}
    \c_{u}(\u, \z)^{-\mathsf{T}} \v_u
    &= \left(\begin{array}{ccc}
        1 & 0 & 0 \\
        -1/z_{1} & 1/z_{1} & 0 \\
        0 & 0 & 1/3u_{3}^{2}
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2} \\ v_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        v_{1} \\ (v_{2} - v_{1})/z_{1} \\ v_{3}/3u_{3}^{2}
    \end{array}\right).
\end{aligned}
$$

::::{admonition} Exercise 5
:class: exercise

Calculate $\c_z(\u,\z)\trp\v_u$, $\c_u(\u,\z)^{-1}\v_u$, and $\c_z(\u,\z)\v_z$.

:::{admonition} Solution
:class: solution dropdown

$$
\begin{aligned}
    \c_{z}(\u, \z)\trp \v_u
    &= \left(\begin{array}{ccc}
        -1 & u_{2} & 0 \\
        0 & -1 & -2z_{2}
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2} \\ v_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        -v_{1} + u_{2}v_{2}
        \\
        -v_{2} - 2z_{2}v_{3}
    \end{array}\right),
    \\ \\
    \c_u(\u,\z)^{-1}\v_u
    &= \left(\begin{array}{ccc}
        1 & -1/z_{1} & 0 \\
        0 & 1/z_{1} & 0 \\
        0 & 0 & 1/3u_{3}^{2}
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2} \\ v_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        v_{1} - v_{2}/z_{1} \\ v_{2}/z_{1} \\ v_{3}/3u_{3}^{2}
    \end{array}\right),
    \\ \\
    \c_z(\u,\z)\v_z
    &= \left(\begin{array}{cc}
        -1 & 0 \\
        u_{2} & -1 \\
        0 & -2z_{2}
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2}
    \end{array}\right)
    =
    \left(\begin{array}{c}
        -v_{1} \\ u_{2}v_{1} - v_{2} \\ -2z_{2}v_{3}
    \end{array}\right).
\end{aligned}
$$

:::
::::

### Constraint Hessians

With these derivatives, we can now solve the constrained optimization problem using {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_gaussnewton`.
However, using {prf:ref}`alg:adjoint_hessvec` instead of {prf:ref}`alg:adjoint_gaussnewton` also requires Hessian actions of the constraint, $\bflambda\trp\c_{u,u}(\u,\z)\v_{u}$, $\bflambda\trp\c_{u,z}(\u,\z)\v_{z}$, $\bflambda\trp\c_{z,u}(\u,\z)\v_{u}$, and $\bflambda\trp\c_{z,z}(\u,\z)\v_{z}$.
The vector $\bflambda = (\lambda_1,\ldots,\lambda_{n_u})\trp\in\R^{n_u}$ is called the _adjoint_ of the state.

The Hessian $\c_{u,u}$ and its siblings are third-order tensors ("three-dimensional matrices") and are difficult to write explicitly.
Instead, we write the equation for the components of the Hessian action: the $j$th entry of the vector-Hessian-vector product $\bflambda\trp\c_{u,u}(\u,\z)\v_u$ is given by

$$
\begin{aligned}
    {[}\![\bflambda\trp\c_{u,u}(\u,\z)\v_u]\!]_{j}
    % &= \sum_{i=1}^{n_{u}}\sum_{k=1}^{n_u}{[}\![\bflambda]\!]_{i}
    % {[}\![\c_{u,u}(\u,\z)]\!]_{ijk} [\![\v_u]\!]_{k}
    = \sum_{i=1}^{n_{u}}\sum_{k=1}^{n_u}
        \lambda_{i}\frac{\partial^{2} c_{i}(\u,\z)}{\partial u_{j}\partial u_{k}} v_k.
\end{aligned}
$$

It is easy to see that $\frac{\partial^{2} c_{i}(\u,\z)}{\partial u_{j}\partial u_{k}} = 0$ except when $i = j = k = 3$, so we have

$$
\begin{aligned}
    \bflambda\trp\c_{u,u}(\u,\z)\v_u
    = \left(\begin{array}{c}
        0 \\ 0 \\ 6\lambda_{3}u_{3}v_{3}
    \end{array}\right).
\end{aligned}
$$

::::{admonition} Exercise 6
:class: exercise

Calculate the Hessian actions $\bflambda\trp\c_{u,z}(\u,\z)\v_z$, $\bflambda\trp\c_{z,u}(\u,\z)\v_u$, and $\bflambda\trp\c_{z,z}(\u,\z)\v_z$.

:::{admonition} Solution
:class: solution dropdown

$$
\begin{aligned}
    \bflambda\trp \c_{u,z}(\u,\z)\v_{z}
    &= \left(\begin{array}{c}
        0 \\ \lambda_{2} v_{1} \\ 0
    \end{array}\right),
    \\ \\
    \bflambda\trp \c_{z,u}(\u,\z)\v_{u}
    &= \left(\begin{array}{c}
        0 \\ \lambda_{1}v_{2}
    \end{array}\right),
    \\ \\
    \bflambda\trp \c_{z,z}(\u,\z)\v_{z}
    &= \left(\begin{array}{c}
        0 \\ -2\lambda_{2}v_{2}
    \end{array}\right).
\end{aligned}
$$

:::
::::

### Constraint Implementation

We now define a new class, `Tutorial_1_Constraint`, inheriting from {class}`Constraint` that implements the abstract methods, see the [inheritance template](optimization.Constraint.template) to get started.
The class must be defined in a file named `Tutorial_1_Constraint.m`.
Note that in the Jacobian methods, `u_in`, and `z_in` are always treated as matrices, not column vectors.

```matlab
% Tutorial_1_Constraint.m

classdef Tutorial_1_Constraint < Constraint

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u2 = z(2) / z(1);
            u = [z(1) - u2; u2; z(2)^(2 / 3)];
        end

        % Jacobian actions.
        function [u_out] = c_u_Transpose_Inverse_Apply(this, u_in, u, z)
            u_out = [u_in(1, :)
                     (u_in(2, :) - u_in(1, :)) / z(1)
                     u_in(3, :) / (3 * u(3)^2)];
        end

        function [z_out] = c_z_Transpose_Apply(this, u_in, u, z)
            error('c_z_Transpose_Apply() not implemented');
        end

        function [u_out] = c_u_Inverse_Apply(this, u_in, u, z)
            error('c_u_Inverse_Apply() not implemented');
        end

        function [u_out] = c_z_Apply(this, z_in, u, z)
            error('c_z_Apply() not implemented');
        end

        % Hessian actions. These methods are not required when
        % Reduced_Space_Optimization.Gauss_Newton_Hess=true.
        function [u_out] = c_uu_Apply(this, u_in, u, z, lambda)
            u_out = [0; 0; 6 * lambda(3) * u(3) * u_in(3)];
        end

        function [u_out] = c_uz_Apply(this, z_in, u, z, lambda)
            error('c_uz_Apply() not implemented');
        end

        function [z_out] = c_zu_Apply(this, u_in, u, z, lambda)
            error('c_zu_Apply() not implemented');
        end

        function [z_out] = c_zz_Apply(this, z_in, u, z, lambda)
            error('c_zz_Apply() not implemented');
        end

        % This method is required for finite difference checks.
        function [con] = c(this, u, z)
            con = [u(1) + u(2) - z(1)
                   z(1) * u(2) - z(2)
                   u(3)^3 - z(2)^2];
        end

    end
end
```

::::{admonition} Exercise 7
:class: exercise

Implement the following Jacobian methods in a vectorized fashion.

- `c_z_Transpose_Apply()` for computing $\c_z(\u,\z)\trp\v_u$
- `c_u_Inverse_Apply()` for computing $\c_u(\u,\z)^{-1}\v_u$
- `c_z_Apply()` for computing $\c_z(\u,\z)\v_z$

In addition, implement the following Hessian methods.
These need not be vectorized.

- `c_uz_Apply()` for computing $\bflambda\trp\c_{u,z}(\u,\z)\v_z$
- `c_zu_Apply()` for computing $\bflambda\trp\c_{z,u}(\u,\z)\v_u$
- `c_zz_Apply()` for computing $\bflambda\trp\c_{z,z}(\u,\z)\v_z$

:::{admonition} Solution
:class: solution dropdown

Jacobian methods:

```matlab
function [z_out] = c_z_Transpose_Apply(this, u_in, u, z)
    z_out = [-u_in(1, :) + u_in(2, :) * u(2)
             -u_in(2, :) - 2 * u_in(3, :) * z(2)];
end

function [u_out] = c_u_Inverse_Apply(this, u_in, u, z)
    u2z = u_in(2, :) / z(1);
    u_out = [u_in(1, :) - u2z
             u2z
             u_in(3, :) / (3 * u(3)^2)];
end

function [u_out] = c_z_Apply(this, z_in, u, z)
    u_out = [-z_in(1, :)
             z_in(1, :) * u(2) - z_in(2, :)
             -2 * z_in(2, :) * z(2)];
end
```

Hessian methods:

```matlab
function [u_out] = c_uz_Apply(this, z_in, u, z, lambda)
    u_out = [0; lambda(2) * z_in(1); 0];
end

function [z_out] = c_zu_Apply(this, u_in, u, z, lambda)
    z_out = [0; lambda(1) * u_in(2)];
end

function [z_out] = c_zz_Apply(this, z_in, u, z, lambda)
    z_out = [0; -2 * lambda(2) * z_in(2)];
end
```

:::
::::

## Implementation Verification

When extending {class}`Objective` and {class}`Constraint` with derivatives calculated by hand, it is a good idea to perform finite difference checks to verify that the derivatives of $J$ and $\c$ were implemented consistently.
We begin by instantiating our `Tutorial_1_Objective` and `Tutorial_1_Constraint` classes in a new file named `Tutorial_1.m`, then call the finite difference check methods for each object.
Here we are use the constants $\alpha_{1} = 7$, $\alpha_{2} = 1$, $\alpha_{3} = 4$, and $\alpha_{4} = \alpha_{5} = 8$.
In this case, the solution to the optimization problem is $\u = (\alpha_1,\alpha_2,\alpha_3)\trp$ and $\z = (\alpha_4,\alpha_5)\trp$.

```matlab
% Tutorial.m

%% Clear workspace and add the SABL optimization source path.
clear;
close all;
clc;
addpath('~/SABL/src/');     % MODIFY THIS TO MATCH YOUR PATH TO SABL.
rng(1342);                  % Random seed for reproducing exact results (optional).

%% Instantiate the objective and constraints.
alphas = [7; 1; 4; 8; 8];
obj = Tutorial_1_Objective(alphas);
con = Tutorial_1_Constraint();

%% Run finite difference checks.
n_u = 3;
n_z = 2;
u0 = rand(n_u, 1) + 2;
z0 = rand(n_z, 1) + 1;

obj.Finite_Difference_Gradient_Check(u0, z0);
obj.Finite_Difference_Hessian_Check(u0, z0);
con.Finite_Difference_Constraint_Check(u0, z0);
```

At this point, running the file should produce a report like the following.

```text
>> run Tutorial_1.m

u gradient finite difference check
h = 0.01 and error = 0.00018433
h = 0.001 and error = 1.8433e-05
h = 0.0001 and error = 1.8433e-06
h = 1e-05 and error = 1.8471e-07
h = 1e-06 and error = 2.0459e-08

z gradient finite difference check
h = 0.01 and error = 0.00026405
h = 0.001 and error = 2.6405e-05
h = 0.0001 and error = 2.6405e-06
h = 1e-05 and error = 2.6409e-07
h = 1e-06 and error = 2.8522e-08

uu Hessian finite difference check
h = 0.01 and error = 1.8186e-13
h = 0.001 and error = 3.4563e-13
h = 0.0001 and error = 5.5373e-12
h = 1e-05 and error = 2.1222e-10
h = 1e-06 and error = 8.2187e-10

uz Hessian finite difference check
h = 0.01 and error = 0.00031617
h = 0.001 and error = 3.1617e-05
h = 0.0001 and error = 3.1617e-06
h = 1e-05 and error = 3.1619e-07
h = 1e-06 and error = 3.1613e-08

zu Hessian finite difference check
h = 0.01 and error = 1.0473e-05
h = 0.001 and error = 1.0473e-06
h = 0.0001 and error = 1.0472e-07
h = 1e-05 and error = 9.416e-09
h = 1e-06 and error = 3.8897e-09

zz Hessian finite difference check
h = 0.01 and error = 2.6526e-13
h = 0.001 and error = 3.8852e-13
h = 0.0001 and error = 3.7054e-11
h = 1e-05 and error = 1.0642e-09
h = 1e-06 and error = 4.7359e-09

Constraint z Jacobian finite difference check
h = 0.01 and error = 5.4324e-05
h = 0.001 and error = 5.4324e-06
h = 0.0001 and error = 5.4323e-07
h = 1e-05 and error = 5.4299e-08
h = 1e-06 and error = 5.4361e-09

Constraint u Jacobian Inverse finite difference check
h = 0.01 and error = 4.2273e-05
h = 0.001 and error = 4.2279e-06
h = 0.0001 and error = 4.228e-07
h = 1e-05 and error = 4.2286e-08
h = 1e-06 and error = 4.2168e-09
```

The finite difference verification methods perform finite differences with varying step sizes.
The error should be small for some or all step sizes, which is what we see in this case.

## Solving the Problem

Now that we are confident that our derivatives are implemented correctly, we couple the objective and constraint by instantiating a {class}`Reduced_Space_Optimization` object.

```matlab
% Tutorial_1.m (continued)

%% Instantiate the optimization problem, solve it, and report the results.
opt = Reduced_Space_Optimization(obj, con);
[u, z] = opt.Optimize(z0);

disp("Experiment 1: alphas = [" + num2str(alphas') + "]");
disp(" ");
disp("State:");
disp(u);
disp("Control:");
disp(z);
disp("Objective:");
disp(obj.J(u, z));
```

```text
>> run Tutorial_1.m

 Iteration        f(x)          step          optimality   CG-iterations
     0            3215.61                           355
Objective function returned complex; trying a new point...
     1            3215.61             10            355           1
     2            3010.95            0.5            439           0
     3            2487.43              1            585           1
     4               1178              2            658           1
     5             547.49        4.00167            776          50
     6             46.374         4.1229            186           2
     7           0.445696        4.01851           11.3           2
     8           0.240943        1.00008           1.07           7
     9          0.0505641           0.25          0.477           5
    10        4.45777e-07       0.212295        0.00579           2
    11        8.76926e-14    0.000606593       6.35e-07           2

Optimization stopped because the relative objective function value is changing
by less than options.FunctionTolerance = 1.000000e-06.

alphas = [7  1  4  8  8]

State:
    7.0000
    1.0000
    4.0000

Control:
    8.0000
    8.0000

Objective:
   8.7693e-14
```

We have successfully found the minimizer $\u^{*} = (7, 1, 4)$, $\z^{*} = (8, 8)$.

:::{admonition} Exercise 8
:class: exercise

By default, {meth}`Reduced_Space_Optimization.Optimize` uses {prf:ref}`alg:adjoint_hessvec` for estimating the Hessian action of the reduced-space objective $\hat{J}(\z) = J(\S(\z),\z)$.
Repeat the above experiment with {prf:ref}`alg:adjoint_gaussnewton` by setting
`opt.Gauss_Newton_Hess = true` and calling `opt.Optimize()` again.
This strategy only uses the Jacobian actions of the constraint, not its Hessian actions, which is often computationally less expensive (but sometimes less accurate).
:::

Let us now solve this problem for a choice of constants $(\alpha_1,\ldots,\alpha_2)$ where the solution is not as obvious.
Consider

$$
\begin{aligned}
    \alpha_1 &= 1,
    &
    \alpha_2 &= 2,
    &
    \alpha_3 &= 3,
    &
    \alpha_4 &= 4,
    &
    \alpha_5 &= 5.
\end{aligned}
$$

We need to re-instantiate the objective and the optimizer, but the constraint is the same as before.

```matlab
% Tutorial_1.m (continued)

%% Instantiate and solve the problem for different alphas.
alphas = [1; 2; 3; 4; 5];
obj = Tutorial_1_Objective(alphas);
opt = Reduced_Space_Optimization(obj, con);

[u, z] = opt.Optimize(z0);

%% Report the results.
disp(" ");
disp("alphas = [" + num2str(alphas') + "]");
disp(" ");
disp("OPTIMIZATION");
disp("------------");
disp("State:");
disp(u);
disp("Control:");
disp(z);
disp("Objective:");
disp(obj.J(u, z));
```

```text
>> run Tutorial_1.m

alphas = [1  2  3  4  5]

OPTIMIZATION
------------
State:
    1.3279
    1.7193
    3.0165

Control:
    3.0472
    5.2391

Objective:
    1.1537
```

:::{admonition} Exercise 7
:class: exercise

Since the control $\z$ has only entries, we can roughly verify the solution to this problem with a dense grid search of the two-dimensional control space.
Loop through $500$ equally spaced points $z_{1}^{(i)}\in[2, 6]$ and again through $500$ equally spaced points $z_{2}^{(j)}\in[2, 6]$.
In the inner loop, set $\z^{(i,j)} = (z_{1}^{(i)},z_{2}^{(j)})$ and solve the constraint to get $\u^{(i,j)} = \S(\z^{(i,j)})$, then evaluate the objective function $J(\u^{(i,j)}, \z^{(i,j)})$.
Compare the control $\z^{(i,j)}$ that approximately minimizes the objective function to the control produced by the constrained optimization, as well as the corresponding objective function values.
:::

## Automatic Differentiation

Calculating and implementing the derivatives of $J$ and $\c$ by hand results in code that executes quickly, but the process can be time consuming and errors can be difficult to locate.
In some cases, automatic differentation tools can be used to implement the abstract methods of {class}`Objective` and {class}`Constraint` by only implementing $J$ and $\c$.

First, we define new classes `Tutorial_1_Objective_AD` and `Tutorial_1_Constraint_AD` that inherit from {class}`Objective_AD` and {class}`Constraint_AD`, respectively.
As before, each class must be defined in a file whose name matches the class name.

```matlab
% Tutorial_1_Objective_AD.m

classdef Tutorial_1_Objective_AD < Objective_AD

    properties
        a   % alpha constants.
    end

    methods (Access = public)

        % Constructor: set the alpha constants.
        function this = Tutorial_1_Objective_AD(a, n_u, n_z)
            this = this@Objective_AD(n_u, n_z);
            this.a = a;
        end

        function [val] = J_AD(this, u, z)
            val = sum((u - this.a(1:3)).^2);
            val = val + sum((z - this.a(4:5)).^2);
            val = val + (u(1) * z(1) - this.a(1) * this.a(4))^2;
        end

    end
end
```

```matlab
% Tutorial_1_Constraint_AD.m

classdef Tutorial_1_Constraint_AD < Constraint_AD

    methods (Access = public)

        function [c] = c_AD(this, u, z)
            c = [u(1) + u(2) - z(1)
                 z(1) * u(2) - z(2)
                 u(3)^3 - z(2)^2];
        end

    end
end
```

These classes are instantiated with the state and control dimensions $n_u$ and $n_z$.
Before the optimization can occur, we need to call each object's `AD_Initialization()` method.
Below, we use `evalc()` to suppress additional output from the automatic differentiation tool.

```matlab
% Tutorial_1.m (continued)

%% Use automatic differentiation classes to solve the optimization problem.
alphas = [7; 1; 4; 8; 8];
obj_AD = Tutorial_1_Objective_AD(alphas, n_u, n_z);
con_AD = Tutorial_1_Constraint_AD(n_u, n_z);
opt_AD = Reduced_Space_Optimization(obj_AD, con_AD);

evalc("obj_AD.AD_Initialization()");
evalc("con_AD.AD_Initialization()");

[u, z] = opt_AD.Optimize(z0);

% Report the results.
disp("alphas = [" + num2str(alphas') + "]");
disp(" ");
disp("State:");
disp(u);
disp("Control:");
disp(z);
disp("Objective:");
disp(obj_AD.J(u, z));
```

```text
>> run Tutorial_1.m

                                Norm of      First-order
 Iteration        f(x)          step          optimality   CG-iterations
     0            3215.61                           355
     1            3215.61             10            355           1
     2            1840.14            2.5            674           0
     3            548.493              5            798           1
     4            7.33805        5.97466           85.4           2
     5         0.00357153       0.690386           1.62           2
     6        3.73168e-09      0.0307208        0.00182           2
     7        4.21046e-21    2.23334e-05       1.86e-09           2

Optimization completed: The first-order optimality measure, 1.864731e-09,
is less than options.OptimalityTolerance = 1.000000e-08, and no negative/zero
curvature is detected in the trust-region model.

alphas = [7  1  4  8  8]

State:
    7.0000
    1.0000
    4.0000

Control:
    8.0000
    8.0000

Objective:
   4.2105e-21
```

SABL uses [ADiGator](https://github.com/matt-weinstein/adigator) for automatic differentiation, which generates several auxiliary files under a new folder `AdiGator_Files/` in the current directory.
The file generation incurs a small up-front cost, but it is usually much faster (when it works) than implementing derivatives by hand.
However, if the objective or constraint are changed, these files need to be regenerated by calling `AD_Initialization()` again.

The next tutorial focuses on another low-dimensional problem, this time with time-dependent differential equation constraints.
