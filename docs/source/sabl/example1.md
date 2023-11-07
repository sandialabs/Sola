# Example 1: System of ODEs

This example considers a simple [constrained optimization problem](../optimization) of the form

$$
\begin{aligned}
    \min_{u,z} ~& J(u,z)
    \\
    s.t. ~~& c(u,z) = 0.
\end{aligned}
$$

We will implement subclasses of [`Objective`](sabl:optimization-objective) and [`Constraint`](sabl:optimization-constraints) and show how to solve the optimization problem with a [`Reduced_Space_Optimization`](sabl:optimizer-class).

## Problem Statement

Let $n_{u} = 3$ be the state dimension, $n_{z} = 2$ be the control dimension, and denote the state and control by $u = (~u_{1}~~u_{2}~~u_{3}~)^{\mathsf{T}}$ and $z = (~z_{1}~~z_{2}~)^{\mathsf{T}}$, respectively.
We consider the objective function

::::{margin}
:::{note}
By inspection, the solution to this problem without any constraints on $u$ and $z$ is $u^{*} = (~\alpha_{1}~~\alpha_{2}~~\alpha_{3}~)^{\mathsf{T}}$, $z^{*} = (~\alpha_{4}~~\alpha_{5}~)^{\mathsf{T}}$.
:::
::::

$$
\begin{aligned}
    J(u, z)
    &= (u_{1} - \alpha_{1})^{2} + (u_{2} - \alpha_{2})^{2} + (u_{3} - \alpha_{3})^{2} \\
    &\quad+ (z_{1} - \alpha_{4})^{2} + (z_{2} - \alpha_{5})^{2} + (u_{1}z_{1} - \alpha_{1}\alpha_{4})^{2},
\end{aligned}
$$

where $\alpha_{1},\alpha_{2},\alpha_{3},\alpha_{4},\alpha_{5}\in\mathbb{R}$ are known constants.
The $n_{c} = 3$ constraints are encoded by the function

$$
\begin{aligned}
    c(u, z)
    &= \left(\begin{array}{c}
        u_1 + u_2 - z_1 \\
        z_1u_2 - z_2 \\
        u_3^3 - z_2^2
    \end{array}\right).
\end{aligned}
$$

When the constraint equation $c(u, z) = 0$ is satisfied, we have the equations

$$
\begin{aligned}
    u_1 + u_2 &= z_1,
    &
    z_1u_2 &= z_2,
    &
    u_3^3 &= z_2^2.
\end{aligned}
$$

## Implementing the Objective

To implement the [abstract methods](tab:objective_abstract) of `Objective` for this particular choice of $J$, we need to calculate the derivatives of $J$.
Recall that the $i$th entry of the $x$-gradient of $J$ is $\frac{\partial J}{\partial x_{1}}$.

$$
\begin{aligned}
    \nabla_{u}J(u, z)
    &= \left(\begin{array}{c}
        \partial J / \partial u_{1} \\
        \partial J / \partial u_{2} \\
        \partial J / \partial u_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        2(u_{1} - \alpha_{1}) + 2(u_{1}z_{1} - \alpha_{1}\alpha_{4})z_{1}
        \\
        2(u_{2} - \alpha_{2})
        \\
        2(u_{3} - \alpha_{3})
    \end{array}\right),
    \\ \\
    \nabla_{z}J(u, z)
    &= \left(\begin{array}{c}
        \partial J / \partial z_{1} \\
        \partial J / \partial z_{2}
    \end{array}\right)
    = \left(\begin{array}{c}
        2(z_{1} - \alpha_{4}) + 2(u_{1}z_{1} - \alpha_{1}\alpha_{4})u_{1}
        \\
        2(z_{2} - \alpha_{5})
    \end{array}\right).
\end{aligned}
$$

Next, we need the action of the Hessians of $J$.
Recall that the $(i,j)$th entry of the $x,y$ Hessian of $J$ is given by $\frac{\partial^{2} J}{\partial x_{i}\partial y_{j}}$.

$$
\begin{aligned}
    \nabla_{u,u}J(u, z)v
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
    \end{array}\right),
    \\ \\
    \nabla_{u,z}J(u, z)v
    &= \left(\begin{array}{cc}
        4u_{1}z_{1} - 2\alpha_{1}\alpha_{4} & 0
        \\ 0 & 0
        \\ 0 & 0
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2}
    \end{array}\right)
    = \left(\begin{array}{c}
        (4u_{1}z_{1} - 2\alpha_{1}\alpha_{4})v_{1} \\ 0 \\ 0
    \end{array}\right),
    \\ \\
    \nabla_{z,u}J(u, z)v
    &= \left(\begin{array}{ccc}
        4u_{1}z_{1} - 2\alpha_{1}\alpha_{4} & 0 & 0
        \\
        0 & 0 & 0
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2} \\ v_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        (4u_{1}z_{1} - 2\alpha_{1}\alpha_{4})v_{1}
        \\
        0
    \end{array}\right),
    \\ \\
    \nabla_{z,z}J(u, z)v
    &= \left(\begin{array}{cc}
        2 + 2u_{1}^{2} & 0
        \\
        0 & 2
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2}
    \end{array}\right)
    = \left(\begin{array}{c}
        (2 + 2u_{1}^{2})v_{1}
        \\ 2v_{2}
    \end{array}\right).
\end{aligned}
$$

We can now begin with the `Optimize` [template](sabl:optimization-template), add a constructor for setting the constants $\alpha_{1},\ldots,\alpha_{5}$, and implement the abstract methods.
This new class must be defined in a file named `Example_1_Objective.m`.

```matlab
classdef Example_1_Objective < Objective

    methods (Access = public)

        % Constructor: set the alpha constants.
        function this = Example_1_Objective(a1, a2, a3, a4, a5)
            this.a = [a1; a2; a3; a4; a5];
        end

        function [val, grad_u, grad_z] = J(this, u, z)
            u1z1_minus_a1a4 = (u(1) * z(1)) - (this.a(1) * this.a(4));

            % Calculate the value of J(u, z).
            val = sum((u - this.a(1:3)).^2);
            val = val + sum((z - this.a(4:5)).^2);
            val = val + u1z1_minus_a1a4^2;

            % Calculate the u gradient of J.
            grad_u = 2 * (u - this.a(1:3));
            grad_u(1) = grad_u(1) + 2 * u1z1_minus_a1a4 * z(1);

            % Calculate the z gradient of J.
            grad_z = 2 * (z - this.a(4:5));
            grad_z(1) = grad_z(1) + 2 * u1z1_minus_a1a4 * u(1);
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = 2 * v;
            Mv(1, :) = Mv(1, :) .* (1 + z(1)^2);
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = zeros(length(u), size(v, 2));
            Mv(1, :) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) .* v(1, :);
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = zeros(length(z), size(v, 2));
            Mv(1, :) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) .* v(1, :);
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = 2 * v;
            Mv(1, :) = Mv(1, :) * (1 + u(1)^2);
        end

    end
end
```

:::{danger}
Note carefully that the `J_xx_Apply()` methods are implemented in a _vectorized_ fashion by treating the input `v` as a matrix with possibly more than one column.
If `v` were always a column vector, the following implementation would be valid:

```matlab
function [Mv] = J_uz_Apply(this, v, u, z)
    Mv = zeros(length(u), 1);
    Mv(1) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) * v(1);
end
```

However, MATLAB's `fminunc()` function may require multiple simultaneous Hessian-vector evaluations, so we cannot assume `v` has only one column.
Compare the above code carefully with the correct implementation.

```matlab
function [Mv] = J_uz_Apply(this, v, u, z)
    Mv = zeros(length(u), size(v, 2));
    Mv(1, :) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) .* v(1, :);
end
```

:::

:::{tip}
Even though we can easily write out the Hessians of $J$ in this case, our implementation avoids explicitly forming $\nabla_{u,u}J$ and its siblings.
Avoid forming full Hessian matrices computationally whenever possible.

To calculate the Hessian action, it may be helpful to write the Hessian-vector multiplication in component notation.
Letting $[\![x]\!]_{i}$ denote the $i$th entry of a vector $x$ and $[\![A]\!]_{ij}$ denote the $(i,j)$th entry of a matrix $A$, we have

$$
\begin{aligned}
    {[}\![\nabla_{u,u}J(u, z)]\!]_{ij}
    &= \frac{\partial^{2}}{\partial u_{i}\partial u_{j}}J(u, z),
    &
    [\![\nabla_{u,z}J(u, z)]\!]_{ij}
    &= \frac{\partial^{2}}{\partial u_{i}\partial z_{j}}J(u, z),
    \\
    [\![\nabla_{z,u}J(u, z)]\!]_{ij}
    &= \frac{\partial^{2}}{\partial z_{i}\partial u_{j}}J(u, z),
    &
    [\![\nabla_{z,z}J(u, z)]\!]_{ij}
    &= \frac{\partial^{2}}{\partial z_{i}\partial z_{j}}J(u, z),
\end{aligned}
$$

so that the Hessian-vector products are given by

$$
\begin{aligned}
    {[}\![\nabla_{u,u}J(u, z)v]\!]_{i}
    &= \sum_{j=1}^{n_{u}}\frac{\partial^{2}}{\partial u_{i}\partial u_{j}}J(u, z) [\![v]\!]_{j},
    &
    [\![\nabla_{u,z}J(u, z)v]\!]_{i}
    &= \sum_{j=1}^{n_{z}}\frac{\partial^{2}}{\partial u_{i}\partial z_{j}}J(u, z) [\![v]\!]_{j},
    \\
    [\![\nabla_{z,u}J(u, z)v]\!]_{i}
    &= \sum_{j=1}^{n_{u}}\frac{\partial^{2}}{\partial z_{i}\partial u_{j}}J(u, z) [\![v]\!]_{j},
    &
    [\![\nabla_{z,z}J(u, z)v]\!]_{i}
    &= \sum_{j=1}^{n_{z}}\frac{\partial^{2}}{\partial z_{i}\partial z_{j}}J(u, z) [\![v]\!]_{j}.
\end{aligned}
$$

We will use this strategy for computing the constraint Hessian actions.
:::

## Implementing the Constraint

To implement a `Constraint` subclass for this problem, we first need a solution operator $S:z\mapsto u$ such that $c(S(z),z) = 0$ for all $z$.
In this case, we can simply solve the equations

$$
\begin{aligned}
    u_1 + u_2 &= z_1,
    &
    z_1u_2 &= z_2,
    &
    u_3^3 &= z_2^2
\end{aligned}
$$

for $u = (~u_1~~u_2~~u_3~)^\mathsf{T}$:

$$
\begin{aligned}
    u_1 &= z_1 - u_2 = z_1 - (z_2 / z_1),
    &
    u_2 &= z_2 / z_1,
    &
    u_3 &= z_2^{2/3}.
\end{aligned}
$$

Hence, the solution operator $S$ is given by

$$
\begin{aligned}
    S(z)
    = \left(\begin{array}{c}
        z_{1} -  (z_{2} / z_{1})\\
        z_{2} / z_{1} \\
        z_{2}^{2/3}
    \end{array}\right).
\end{aligned}
$$

Next, we calculate the derivatives of $c$.
We start with $c_u$ and $c_z$, which are Jacobian matrices.
Writing

$$
\begin{aligned}
    c(u,z)
    = \left(\begin{array}{c}
        c_{1}(u,z) \\ \vdots \\ c_{n_{c}}(u,z)
    \end{array}\right)
    = \left(\begin{array}{c}
        u_1 + u_2 - z_1 \\
        z_1u_2 - z_2 \\
        u_3^3 - z_2^2
    \end{array}\right),
\end{aligned}
$$

the $(i,j)$th entry of the $x$ Jacobian of $c$ is $\frac{\partial c_{i}}{\partial x_{j}}$.
For this problem, we have

$$
\begin{aligned}
    c_{u}(u, z)
    &= \left(\begin{array}{ccc}
        1 & 1 & 0 \\
        0 & z_{1} & 0 \\
        0 & 0 & 3u_{3}^{2}
    \end{array}\right),
    &
    c_{z}(u, z)
    &= \left(\begin{array}{cc}
        -1 & 0 \\
        u_{2} & -1 \\
        0 & -2z_{2}
    \end{array}\right).
\end{aligned}
$$

Hence, for a vector $\lambda = (~\lambda_{1}~~\lambda_{2}~~\lambda_{3}~)^{\mathsf{T}}$, we have

$$
\begin{aligned}
    c_{u}(u, z)^{-\mathsf{T}} \lambda
    &= \left(\begin{array}{ccc}
        1 & 0 & 0 \\
        -1/z_{1} & 1/z_{1} & 0 \\
        0 & 0 & 1/3u_{3}^{2}
    \end{array}\right)
    \left(\begin{array}{c}
        \lambda_{1} \\ \lambda_{2} \\ \lambda_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        \lambda_{1} \\ (\lambda_{2} - \lambda_{1})/z_{1} \\ \lambda_{3}/3u_{3}^{2}
    \end{array}\right),
    \\ \\
    c_{z}(u, z)^{\mathsf{T}} \lambda
    &= \left(\begin{array}{ccc}
        -1 & u_{2} & 0 \\
        0 & -1 & -2z_{2}
    \end{array}\right)
    \left(\begin{array}{c}
        \lambda_{1} \\ \lambda_{2} \\ \lambda_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        -\lambda_{1} + \lambda_{2}u_{2}
        \\
        -\lambda_{2} - 2\lambda_{3}z_{2}
    \end{array}\right).
\end{aligned}
$$

For the Hessian action of $\hat{J}$, we also need the following for a vector $v = (~v_{1}~~v_{2}~)^{\mathsf{T}}$.

$$
\begin{aligned}
    c_{u}(u, z)^{-1} \lambda
    &= \left(\begin{array}{ccc}
        1 & -1/z_{1} & 0 \\
        0 & 1/z_{1} & 0 \\
        0 & 0 & 1/3u_{3}^{2}
    \end{array}\right)
    \left(\begin{array}{c}
        \lambda_{1} \\ \lambda_{2} \\ \lambda_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        \lambda_{1} - \lambda_{2}/z_{1} \\
        \lambda_{2}/z_{1} \\
        \lambda_{3}/3u_{3}^{2}
    \end{array}\right),
    \\
    c_{z}(u, z) v
    &= \left(\begin{array}{cc}
        -1 & 0 \\
        u_{2} & -1 \\
        0 & -2z_{2}
    \end{array}\right)
    \left(\begin{array}{c}
        v_{1} \\ v_{2}
    \end{array}\right)
    = \left(\begin{array}{c}
        -v_{1} \\
        v_{1}u_{2} - v_{2} \\
        -2v_{2}z_{2}
    \end{array}\right).
\end{aligned}
$$

In addition, we must calculate the action of $c_{u,u}$, $c_{u,z}$, $c_{z,u}$, and $c_{z,z}$.
These Hessians are third-order tensors ("three-dimensional matrices") and are difficult to write out explicitly.
Instead, we write the equation for the components of the Hessian action.
Note that, in this case, most entries of these Hessian actions will be zero: for instance, nonzero entries of $c_{u,z}(u,z)$ must correspond to $\frac{\partial^{2}c_{2}}{\partial u_{2} \partial z_{1}}$.

$$
\begin{aligned}
    {[}\![\lambda^{\mathsf{T}} c_{u,u}(u,z)\mu]\!]_{j}
    = \sum_{i,j}\frac{\partial^{2}c(u,z)_{i}}{\partial u_{j}\partial u_{k}}\lambda_{i}\mu_{k}
    \quad\Longrightarrow\quad
    \lambda^{\mathsf{T}} c_{u,u}(u,z)\mu
    &= \left(\begin{array}{c}
        0 \\ 0 \\ 6\lambda_{3}\mu_{3}u_{3}
    \end{array}\right)
    \\ \\
    [\![\lambda^{\mathsf{T}} c_{u,z}(u,z)v]\!]_{j}
    = \sum_{i,j}\frac{\partial^{2}c(u,z)_{i}}{\partial u_{j}\partial z_{k}}\lambda_{i}v_{k}
    \quad\Longrightarrow\quad
    \lambda^{\mathsf{T}} c_{u,z}(u,z)v
    &= \left(\begin{array}{c}
        0 \\ \lambda_{2} v_{1} \\ 0
    \end{array}\right)
    \\ \\
    [\![\lambda^{\mathsf{T}} c_{z,u}(u,z)\mu]\!]_{j}
    = \sum_{i,j}\frac{\partial^{2}c(u,z)_{i}}{\partial z_{j}\partial u_{k}}\lambda_{i}\mu_{k}
    \quad\Longrightarrow\quad
    \lambda^{\mathsf{T}} c_{z,u}(u,z)\mu
    &= \left(\begin{array}{c}
        0 \\ \lambda_{1}\mu_{2}
    \end{array}\right)
    \\ \\
    [\![\lambda^{\mathsf{T}} c_{z,z}(u,z)v]\!]_{j}
    = \sum_{i,j}\frac{\partial^{2}c(u,z)_{i}}{\partial z_{j}\partial z_{k}}\lambda_{i}v_{k}
    \quad\Longrightarrow\quad
    \lambda^{\mathsf{T}} c_{z,z}(u,z)v
    &= \left(\begin{array}{c}
        0 \\ -2\lambda_{2}v_{2}
    \end{array}\right).
\end{aligned}
$$

With these derivatives calculated, we may now define a `Constraint` subclass that implements the abstract methods by filling out the `Constraint` [template](sabl:constraint-template).
This new class must be defined in a file named `Example_1_Constraint.m`.

```matlab
classdef Example_1_Constraint < Constraint

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u2 = z(2) / z(1);
            u = [z(1) - u2; u2; z(2)^(2 / 3)];
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = [v(1, :)
                  (v(2, :) - v(1, :)) / z(1)
                  v(3, :) / (3 * u(3)^2)];
        end

        function [Mv] = c_z_Transpose_Apply(this, lambda, u, z)
            Mv = [-lambda(1, :) + lambda(2, :) * u(2)
                  -lambda(2, :) - 2 * lambda(3, :) * z(2)];
        end

        function [Mv] = c_u_Inverse_Apply(this, lambda, u, z)
            Mv2 = lambda(2, :) / z(1);
            Mv = [lambda(1, :) - Mv2
                  Mv2
                  lambda(3, :) / (3 * u(3)^2)];
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            Mv = [-v(1, :)
                  v(1, :) * u(2) - v(2, :)
                  -2 * v(2, :) * z(2)];
        end

        function [Mv] = c_uu_Apply(this, lambda, u, z, v)
            Mv = [0; 0; 6 * lambda(3) * u(3) * v(3)];
        end

        function [Mv] = c_uz_Apply(this, lambda, u, z, v)
            Mv = [0; lambda(2) * v(1); 0];
        end

        function [Mv] = c_zu_Apply(this, lambda, u, z, v)
            Mv = [0; lambda(1) * v(2)];
        end

        function [Mv] = c_zz_Apply(this, lambda, u, z, v)
            Mv = [0; -2 * lambda(2) * v(2)];
        end

    end
end
```

## Implementation Verification

We can now instantiate our `Example_1_Objective` and `Example_1_Constraint` classes and couple them with a `Reduced_Space_Optimization` object.
For this demonstration, we use the constants $\alpha_{1} = 7$, $\alpha_{2} = 1$, $\alpha_{3} = 4$, and $\alpha_{4} = \alpha_{5} = 8$.
Save the following in a new file called `Example_1_Driver.m`.

```matlab
%% Clear workspace and add the SABL optimization source path.
clear;
close all;
clc;
addpath('~/SABL/src/optimization/');             % MODIFY THIS TO MATCH YOUR PATH TO SABL.

%% Instantiate the optimization problem.
obj = Example_1_Objective(7, 1, 4, 8, 8);
con = Example_1_Constraint();
opt = Reduced_Space_Optimization(obj, con);
```

Before solving the optimization problem, it is a good idea to run finite difference checks to verify that the derivatives of $J$ and $c$ were implemented consistently.
The `Objective` and `Reduced_Space_Optimization` classes both have these verification methods.

```matlab
%% Run finite difference checks for the objective function only.
u0 = rand(3, 1) + 2;
z0 = rand(2, 1) + 1;
obj.Finite_Difference_Gradient_Check(u0, z0);
obj.Finite_Difference_Hessian_Check(u0, z0);

%% Run finite difference checks for the optimization problem as a whole.
opt.verbose = true;
z0 = rand(2, 1) + 1;
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);
```

Running this file should produce a report like the following.

```text
>> run Driver_Example_1

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

Gradient finite difference check
h = 0.01 and error = 0.0012923
h = 0.001 and error = 0.00012917
h = 0.0001 and error = 1.2916e-05
h = 1e-05 and error = 1.2915e-06
h = 1e-06 and error = 1.2831e-07

Hessian finite difference check
h = 0.01 and error = 0.97609
h = 0.001 and error = 0.97611
h = 0.0001 and error = 0.97611
h = 1e-05 and error = 0.97611
h = 1e-06 and error = 0.97611
```

:::{warning}
The rest of this page is under construction, please check back later.
:::

:::{admonition} **TODO**
Discuss finite difference analysis, what successful / unsuccessful checks look like.
:::

## Solving the Problem

```matlab
%% Do the optimization.
[u, z] = opt.Optimize(z0);

disp("STATE:");
disp(u);
disp("CONTROL:");
disp(z);
```

```text
>> run Driver_Example_1

                                Norm of      First-order
 Iteration        f(x)          step          optimality   CG-iterations
     0            3005.53                           443
Objective function returned complex; trying a new point...
     1            3005.53             10            443           1
     2            2759.54            0.5            522           0
     3            2168.38              1            638           1
     4            846.255              2            610           1
     5            846.255        4.07401            610          50
     6             319.29        1.00753            415           0
     7            88.3561        2.22474            233          12
     8            15.4327        2.02985           21.7           2
     9            2.92094        2.00047           3.55           2
    10          0.0348141        1.42759          0.422           5
    11        5.33041e-05       0.182098         0.0156           3
    12        3.99466e-12     0.00687426       4.77e-06           2
    13        4.89872e-25    1.87781e-06       1.55e-12           2

Optimization completed: The first-order optimality measure, 1.546467e-12,
is less than options.OptimalityTolerance = 1.000000e-08, and no negative/zero
curvature is detected in the trust-region model.

STATE:
    7.0000
    1.0000
    4.0000

CONTROL:
    8.0000
    8.0000
```

We have successfully found the minimizer $u^{*} = (7, 1, 4)$, $z^{*} = (8, 8)$!

The next example explores a time-dependent problem.
