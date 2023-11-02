# main.md

## Problem Formulation

The SABL library is designed to solve constrained optimization problems of the form

$$
\begin{aligned}
    \begin{aligned}
    \min_{u,z} ~& J(u,z)
    \\
    s.t. ~~& c(u,z) = 0
    \end{aligned}

\end{aligned}
$$ (eqn:opt_prob)

via a reduced space optimization strategy.
Here, $u$ indicates the state variable, $z$ is the control variable, $J$ is the objective function, and $c$ specifies the constraints.
If $c(u,z)=0$ admits a unique solution for any given $z$, then there exists a solution operator $S(z)$ such that $c(S(z),z)=0$ $\forall z$. We may reformulate~{eq}`eqn:opt_prob` in reduced space and solve the equivalent optimization problem

$$
\begin{aligned}
    \min_{z } \hat{J}(z)=J(S(z),z),

\end{aligned}
$$ (eqn:rs_opt_prob)

which will be the focus of this document and the corresponding MATLAB codes. Along with solving~{eq}`eqn:rs_opt_prob`, we would also like the code to facilitate the use of hyper-differential sensitivity analysis (HDSA) for the solution of~{eq}`eqn:rs_opt_prob`, as well as for other outer-loop problems.

Note that~{eq}`eqn:opt_prob` and {eq}`eqn:rs_opt_prob` may be posed as infinite-dimensional optimizations over function spaces.
After discretization for computation, we have finite-dimensional $u\in\mathbb{R}^{n_u}$ and $z\in\mathbb{R}^{n_z}$.
With these dimensions, we have $J:\mathbb{R}^{n_{u}}\times\mathbb{R}^{n_{z}}\to\mathbb{R}$, $c:\mathbb{R}^{n_{u}}\times\mathbb{R}^{n_{z}}\to\mathbb{R}^{n_{c}}$, $S:\mathbb{R}^{n_{z}}\to\mathbb{R}^{n_{u}}$, and $\hat{J}:\mathbb{R}^{n_{z}}\to\mathbb{R}$, where $n_{c}$ is the number of equations needed to express the constraints.
In problems where the constraints are described by a partial differential equation, the state dimension $n_{u}$ and/or the control dimension $n_{z}$ can be very large.

## Adjoint-based Optimization

The unconstrained problem {eq}`eqn:rs_opt_prob` can be solved with off-the-shelf minimizers as long as the gradient $\nabla_z \hat{J}(z)\in\mathbb{R}^{n_{z}}$ and Hessian-vector products $v\mapsto \nabla_{z,z} \hat{J}(z) v$ can be computed efficiently.
To that end, we utilize adjoint-based derivative formulas.
{prf:ref}`alg:adjoint_gradient` shows how to efficiently calculate the gradient, and {prf:ref}`alg:adjoint_hessvec` shows how to compute Hessian-vector products using incremental state and incremental adjoint equations, which avoids explicitly forming the (very large) Hessian matrix $\nabla_{z,z}\hat{J}(z)\in\mathbb{R}^{n_{z}\times n_{z}}$.

In the algorithms and throughout this document, we use $c_u(u,z)\in\mathbb{R}^{n_{c}\times n_{u}}$ and $c_z(u,z)\in\mathbb{R}^{n_{c}\times n_{z}}$ to denote the Jacobians of $c$ with respect to $u$ and $z$, respectively, $\nabla_u J(u, z)\in\mathbb{R}^{n_{u}}$ and $\nabla_z J(u, z)\in\mathbb{R}^{n_{z}}$ to denote the gradients (where $u$ and $z$ are independent of one another) of $J$ with respect to $u$ and $z$, respectively, and a superscript $\mathsf{T}$ to denote a matrix transpose.
We also use $\lambda^{\mathsf{T}} c_{u,u}(u,z) \in \mathbb{R}^{n_{u}\times n_{u}}$ and $\nabla_{u,u} J(u,z)\in\mathbb{R}^{n_{u}\times n_{u}}$ to denote the $(u,u)$ Hessian of the scalar functions $\lambda^{\mathsf{T}} c$ and $J$, respectively, and similar expressions for the $(u,z)$, $(z,u)$, and $(z,z)$ Hessians.

\begin{algorithm}[H]
    \caption{Adjoint-based gradient calculation}
    \label{alg:adjoint_gradient}
    **Input: ** $\overline{z}\in\mathbb{R}^{n_{z}}$\\
    **1: ** Solve the state equation to determine $\overline{u}\in\mathbb{R}^{n_{u}}$ such that $c(\overline{u},\overline{z})=0$\\
    **2: ** Solve the adjoint equation $c_u(\overline{u},\overline{z})^{\mathsf{T}} \overline{\lambda} = - \nabla_u J(\overline{u},\overline{z})$ for $\overline{\lambda}\in\mathbb{R}^{n_{c}}$ (a linear solve) \\
    **3: ** Compute $\nabla_z \hat{J}(\overline{z}) = c_z(\overline{u},\overline{z})^{\mathsf{T}} \overline{\lambda} + \nabla_z J(\overline{u},\overline{z})$ \\
    **Return: ** $\nabla_z \hat{J}(\overline{z}) \in \mathbb{R}^{n_{z}}$
\end{algorithm}

\begin{algorithm}[H]
\caption{Adjoint-based Hessian-vector product calculation}
\label{alg:adjoint_hessvec}
    **Input: ** $v\in\mathbb{R}^{n_{z}}$, $\overline{u}\in\mathbb{R}^{n_{u}}$, $\overline{z}\in\mathbb{R}^{n_{z}}$, $\overline{\lambda}\in\mathbb{R}^{n_{c}}$
    \\
    **1: ** Compute $w = c_z(\overline{u},\overline{z}) v \in\mathbb{R}^{n_{c}}$
    \\
    **2: ** Solve the incremental state equation $c_u(\overline{u},\overline{z}) \overline{\mu} = - w$ for $\overline{\mu} \in \mathbb{R}^{n_{u}}$ (a linear solve)
    \\
    **3: ** Compute $y_J = \nabla_{u,u} J(\overline{u},\overline{z}) \overline{\mu} + \nabla_{u,z} J(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{u}}$
    \\
    **4: ** Compute $y_c = \overline{\lambda}^{\mathsf{T}} c_{u,u}(\overline{u},\overline{z}) \overline{\mu} +  \overline{\lambda}^{\mathsf{T}} c_{u,z}(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{u}}$
    \\
    **5: ** Solve the incremental adjoint equation $c_u(\overline{u},\overline{z})^{\mathsf{T}} \overline{\gamma} = -(y_J + y_c)$ for $\overline{\gamma} \in \mathbb{R}^{n_{c}}$ (a linear solve)
    \\
    **6: ** Compute $x_J= \nabla_{z,u}J(\overline{u},\overline{z}) \overline{\mu} +  \nabla_{z,z}J(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{z}}$
    \\
    **7: ** Compute $x_c = c_z(\overline{u},\overline{z})^{\mathsf{T}} \overline{\gamma} + \overline{\lambda}^{\mathsf{T}} c_{z,u}(\overline{u},\overline{z}) \overline{\mu} + \overline{\lambda}^{\mathsf{T}} c_{z,z}(\overline{u},\overline{z})v \in\mathbb{R}^{n_{z}}$
    \\
    **8: ** Compute $\nabla_{z,z} \hat{J}(\overline{z}) v = x_J + x_c$
    \\
    **Return: ** $\nabla_{z,z} \hat{J}(\overline{z}) v \in \mathbb{R}^{n_{z}}$
\end{algorithm}

## Object-oriented Implementation

To implement Algorithms \ref{alg:adjoint_gradient}--\ref{alg:adjoint_hessvec} and solve {eq}`eqn:opt_prob` via {eq}`eqn:rs_opt_prob`, we adopt an object-oriented design in MATLAB by defining base classes for the optimization objective $J$, the constraints $c$, and the optimization problem as a whole.
Each of these classes defines pure-virtual functions that must be implemented in derived classes to solve specific instances of {eq}`eqn:opt_prob`.

### Objective Function

The `Objective` class is designed to encode the objective functional $J$ and its derivatives.
To represent a particular objective function with the `Objective` class, define a new class that inherits from `Objective` and implements the pure functions listed in [Table TODO](tab:objective_virtuals).

\begin{table}[!h]
\centering
\begin{tabular}{|l|l|}
    \hline
    Function Signature & Mathematical Description
    \\ \hline
    `[val, grad_u, grad_z] = J(u, z)`
    & Evaluate $J$, $\nabla_u J$, and $\nabla_z J$
    \\ \hline
    `[Mv] = J_uu_Apply(v, u, z)`
    & Compute the product $\nabla_{u,u} J(u, z) v$
    \\
    `[Mv] = J_uz_Apply(v, u, z)`
    & Compute the product $\nabla_{u,z} J(u, z) v$
    \\
    `[Mv] = J_zu_Apply(v, u, z)`
    & Compute the product $\nabla_{z,u} J(u, z) v$
    \\
    `[Mv] = J_zz_Apply(v, u, z)`
    & Compute the product $\nabla_{z,z} J(u, z) v$
    \\ \hline
\end{tabular}
\caption{Pure virtual functions of the `Objective` class. The `J()` function is used in {prf:ref}`alg:adjoint_gradient` for calculating the gradient of $\hat{J}$; the others are used in {prf:ref}`alg:adjoint_hessvec` for Hessian-vector products of $\hat{J}$.}
\label{tab:objective_virtuals}
\end{table}

{\color{red}
WARNING! The `J_xx_Apply()` functions must be implemented in a vectorized fashion! i.e., we cannot assume that `v` is $n \times 1$, it may be $n \times X$ for some $X > 1$ because of how MATLAB's `fminunc()` does things.
}

Here is a template to start with a new `Objective` class.

```matlab
classdef MyObjective < Objective

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            error('J() not implemented');
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            error('J_uu_Apply() not implemented');
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            error('J_uz_Apply() not implemented');
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            error('J_zu_Apply() not implemented');
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            error('J_zz_Apply() not implemented');
        end

    end
end
```

#### Example: Objective

Consider the following minimization objective for a problem with state dimension $n_{u} = 3$ and control dimension $n_{z} = 2$,

$$
\begin{aligned}
    J(u, z)
    &= (u_{1} - \alpha_{1})^{2} + (u_{2} - \alpha_{2})^{2} + (u_{3} - \alpha_{3})^{2}
    + (z_{1} - \alpha_{4})^{2} + (z_{2} - \alpha_{5})^{2} + (u_{1}z_{1} - \alpha_{1}\alpha_{4})^{2},
\end{aligned}
$$

where $u = (~u_{1}~~u_{2}~~u_{3}~)^{\mathsf{T}}$ is the state, $z = (~z_{1}~~z_{2}~)^{\mathsf{T}}$ is the control, and $\alpha_{1},\alpha_{2},\alpha_{3},\alpha_{4},\alpha_{5}\in\mathbb{R}$ are known constants.
By inspection, the solution to this problem (if there are no constraints on $u$ and $z$) is $u^{*} = (~\alpha_{1}~~\alpha_{2}~~\alpha_{3}~)^{\mathsf{T}}$, $z^{*} = (~\alpha_{4}~~\alpha_{5}~)^{\mathsf{T}}$.
    % \min_{z \in R^2} J(u,z) = || S(z) - (7,1,4)^T ||^2 + (z_1-8)^2 + (z_2-8)^2 + (u(1)z(1)-56)^2
    %
    % s.t
    %
    % u = S(z)
    %
    % solves
    %
    % u_1 + u_2     + 0     = z_1
    %
    % 0   + z_1u_2 + 0     = z_2
    %
    % 0   +  0      + u_3^3 = z_2^2
The gradients of $J$ with respect to $u$ and $z$ are, respectively,

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
    \\
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

The products of the Hessians of $J$ are given by

$$
\begin{aligned}
    \nabla_{u,u}J(u, z)\mu
    &= \left(\begin{array}{ccc}
        2 + 2z_{1}^{2} & 0 & 0
        \\
        0 & 2 & 0
        \\
        0 & 0 & 2
    \end{array}\right)
    \left(\begin{array}{c}
        \mu_{1} \\ \mu_{2} \\ \mu_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        (2 + 2z_{1}^{2})\mu_{1}
        \\
        2\mu_{2}
        \\
        2\mu_{3}
    \end{array}\right),
    \\
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
    \\
    \nabla_{z,u}J(u, z)\mu
    &= \left(\begin{array}{ccc}
        4u_{1}z_{1} - 2\alpha_{1}\alpha_{4} & 0 & 0
        \\
        0 & 0 & 0
    \end{array}\right)
    \left(\begin{array}{c}
        \mu_{1} \\ \mu_{2} \\ \mu_{3}
    \end{array}\right)
    = \left(\begin{array}{c}
        (4u_{1}z_{1} - 2\alpha_{1}\alpha_{4})\mu_{1}
        \\
        0
    \end{array}\right),
    \\
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

\noindent**Remark:**
In this example we are able to write the Hessians of $J$ explicitly, in part because $n_{u}$ and $n_{z}$ are small.
In more general cases, it may be helpful to write the equation for the components of the Hessian actions.
Let $[\![x]\!]_{i}$ denote the $i$th entry of a vector $x$ and $[\![A]\!]_{ij}$ denote the $(i,j)$th entry of a matrix $A$.
We have

$$
\begin{aligned}
    [\![\nabla_{u,u}J(u, z)]\!]_{ij}
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
    [\![\nabla_{u,u}J(u, z)\mu]\!]_{i}
    &= \sum_{j=1}^{n_{u}}\frac{\partial^{2}}{\partial u_{i}\partial u_{j}}J(u, z)[\![\mu_{j}]\!],
    &
    [\![\nabla_{u,z}J(u, z)v]\!]_{i}
    &= \sum_{j=1}^{n_{z}}\frac{\partial^{2}}{\partial u_{i}\partial z_{j}}J(u, z)[\![v_{j}]\!],
    \\
    [\![\nabla_{z,u}J(u, z)\mu]\!]_{i}
    &= \sum_{j=1}^{n_{u}}\frac{\partial^{2}}{\partial z_{i}\partial u_{j}}J(u, z)[\![\mu_{j}]\!],
    &
    [\![\nabla_{z,z}J(u, z)v]\!]_{i}
    &= \sum_{j=1}^{n_{z}}\frac{\partial^{2}}{\partial z_{i}\partial z_{j}}J(u, z) [\![v_{j}]\!].
\end{aligned}
$$

In the next section we use _Einstein notation_, which omits the summation signs.

\vspace{.25in}

The following class implements this example objective function as a subclass of `Objective`.

\lstinputlisting[style=Matlab-editor,frame=single,numbers=left]{../../tests/optimization/Example_1/Example_1_Objective.m}

In addition to the functions listed in [Table TODO](tab:objective_virtuals), the `Objective` class is equipped with finite difference checkers to verify the consistency between the gradient and Hessian functions.
% \begin{table}[!ht]
% \centering
% \begin{tabular}{|l|l|}
%     \hline
%     Function Signature & Description
%     \\ \hline
%     `[diffs_u, diffs_z] = Finite_Difference_Gradient_Check(u, z)`
%     & Check gradients of $J$.
%     \\
%     `[diffs_uu, diffs_uz, diffs_zu, diffs_zz] = Finite_Difference_Hessian_Check(u, z)`
%     & Check Hessian-vector products of $J$.
%     \\ \hline
% \end{tabular}
% \caption{Verification functions in the `Objective` class.}
% \label{tab:objective-checkers}
% \end{table}
Below, we instantiate our example `Objective` listed above and run the finite difference checks.

```matlab
TODO;
```

### Constraints

The `Constraint` class encodes the constraint function $c(u, z) = 0$, its derivatives, and the solution operator $S:z\mapsto u$.
To represent a particular set of constraints with the `Constraint` class, define a new class that inherits from `Constraint` and implements the pure virtual functions listed in [Table TODO](tab:constraint_virtuals).

\begin{table}[!ht]
\centering
\begin{tabular}{|l|l|}
    \hline
    Function Signature & Mathematical Description
    \\ \hline
    `[u] = State_Solve(z)`
    & Solve $c(u, z)=0$ for $u$
    \\
    `[Mv] = c_u_Transpose_Inverse_Apply(v, u, z)`
    & Compute the product $c_u(u, z)^{-\mathsf{T}} v$
    \\
    `[Mv] = c_z_Transpose_Apply(v, u, z)`
    & Compute the product $c_z(u, z)^{\mathsf{T}} v$
    \\ \hline
    `[Mv] = c_u_Inverse_Apply(v, u, z)`
    & Compute the product $c_u(u, z)^{-1} v$
    \\
    `[Mv] = c_z_Apply(v, u, z)`
    & Compute the product $c_z(u, z) v$
    \\
    `[Mv] = c_uu_Apply(v, u, z, lambda)`
    & Compute the product $\lambda^{\mathsf{T}} c_{u, u}(u, z) v$
    \\
    `[Mv] = c_uz_Apply(v, u, z, lambda)`
    & Compute the product $\lambda^{\mathsf{T}} c_{u, z}(u, z) v$
    \\
    `[Mv] = c_zu_Apply(v, u, z, lambda)`
    & Compute the product $\lambda^{\mathsf{T}} c_{z, u}(u, z)v$
    \\
    `[Mv] = c_zz_Apply(v, u, z, lambda)`
    & Compute the product $\lambda^{\mathsf{T}} c_{z, z}(u, z) v$
    \\ \hline
\end{tabular}
\caption{Pure virtual functions of the `Constraint` class. The top block of functions are used in {prf:ref}`alg:adjoint_gradient` for calculating the gradient of $\hat{J}$; the lower block of functions are used in {prf:ref}`alg:adjoint_hessvec` for Hessian-vector products of $\hat{J}$.}
\label{tab:constraint_virtuals}
\end{table}

Here is a template to start with a new `Constraint` class.

```matlab
classdef MyConstraint < Constraint

    methods (Access = public)

        function [u] = State_Solve(this, z)
            error('StateSolve() not implemented');
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            error('c_u_Transpose_Inverse_Apply() not implemented');
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            error('c_z_Transpose_Apply() not implemented');
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            error('c_u_Inverse_Apply() not implemented');
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            error('c_z_Apply() not implemented');
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            error('c_uu_Apply() not implemented');
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            error('c_uz_Apply() not implemented');
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            error('c_zu_Apply() not implemented');
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            error('c_zz_Apply() not implemented');
        end

    end
end
```

#### Example: Constraint

Consider again a problem with state dimension $n_{u} = 3$ and control dimension $n_{z} = 2$ with $n_{c}=3$ constraints encoded by the function $c : \mathbb{R}^{n_{u}}\times\mathbb{R}^{n_{z}}\to\mathbb{R}^{n_{c}}$ given by

$$
\begin{aligned}
    c(u, z)
    &= \left(\begin{array}{c}
        u_1 + u_2 - z_1 \\
        z_1u_2 - z_2 \\
        u_3^3 - z_2^2
    \end{array}\right),
\end{aligned}
$$

where $u = (~u_{1}~~u_{2}~~u_{3}~)^{\mathsf{T}}$ is the state and $z = (~z_{1}~~z_{2}~)^{\mathsf{T}}$ is the control.
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

For a given $z$, the solution operator $S:\mathbb{R}^{n_{z}}\to\mathbb{R}^{n_{u}}$ produces the state such that $c(S(z), z) = 0$.
In this case, we can construct $S$ by solving the preceding equations directly for $u$:

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

Next, we determine the derivatives of $c$.
We start with $c_u$ and $c_z$, which are Jacobians.

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

Hence, a vector $\lambda = (~\lambda_{1}~~\lambda_{2}~~\lambda_{3}~)^{\mathsf{T}}$, we have

$$
\begin{aligned}
    c_{u}(u, z)^{-\mathsf{T}} \lambda
    % &=
    % \left(\begin{array}{ccc}
    %     1 & 0 & 0 \\
    %     1 & z_{1} & 0 \\
    %     0 & 0 & 3u_{3}^{2}
    % \end{array}\right)^{-1}
    % \left(\begin{array}{c}
    %     \lambda_{1} \\ \lambda_{2} \\ \lambda_{3}
    % \end{array}\right)
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
    \\
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

We can now implement `c_u_Transpose_Inverse_Apply()` and
`c_z_Transpose_Apply()` for the computation of the gradient of $\hat{J}$.
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
These Hessians are third-order tensors (``three-dimensional matrices'') and are difficult to write out explicitly.
Instead, here we use Einstein summation notation.
Note that, in this case, most entries of these Hessian actions will be zero: for instance, nonzero entries of $c_{u,z}(u,z)$ must correspond to $\frac{\partial^{2}}{\partial u_{2} \partial z_{1}}[\![c]\!]_{2}$.

$$
\begin{aligned}
    [\![\lambda^{\mathsf{T}} c_{u,u}(u,z)\mu]\!]_{j}
    = \frac{\partial}{\partial u_{j}\partial u_{k}}[\![c(u,z)]\!]_{i}\lambda_{i}\mu_{k}
    \quad\Longrightarrow\quad
    \lambda^{\mathsf{T}} c_{u,u}(u,z)\mu
    &= \left(\begin{array}{c}
        0 \\ 0 \\ 6\lambda_{3}\mu_{3}u_{3}
    \end{array}\right)
    \\
    [\![\lambda^{\mathsf{T}} c_{u,z}(u,z)v]\!]_{j}
    = \frac{\partial}{\partial u_{j}\partial z_{k}}[\![c(u,z)]\!]_{i}\lambda_{i}v_{k}
    \quad\Longrightarrow\quad
    \lambda^{\mathsf{T}} c_{u,z}(u,z)v
    &= \left(\begin{array}{c}
        0 \\ \lambda_{2} v_{1} \\ 0
    \end{array}\right)
    \\
    [\![\lambda^{\mathsf{T}} c_{z,u}(u,z)\mu]\!]_{j}
    = \frac{\partial}{\partial z_{j}\partial u_{k}}[\![c(u,z)]\!]_{i}\lambda_{i}\mu_{k}
    \quad\Longrightarrow\quad
    \lambda^{\mathsf{T}} c_{z,u}(u,z)\mu
    &= \left(\begin{array}{c}
        0 \\ \lambda_{1}\mu_{2}
    \end{array}\right)
    \\
    [\![\lambda^{\mathsf{T}} c_{z,z}(u,z)v]\!]_{j}
    = \frac{\partial}{\partial z_{j}\partial z_{k}}[\![c(u,z)]\!]_{i}\lambda_{i}v_{k}
    \quad\Longrightarrow\quad
    \lambda^{\mathsf{T}} c_{z,z}(u,z)v
    &= \left(\begin{array}{c}
        0 \\ -2\lambda_{2}v_{2}
    \end{array}\right).
\end{aligned}
$$

With these computed, we can finally implement this constraint.

\lstinputlisting[style=Matlab-editor,frame=single,numbers=left]{../../tests/optimization/Example_1/Example_1_Constraint.m}

{\color{red}
WARNING! The `c_x_Apply()` functions must be implemented in a vectorized fashion! i.e., we cannot assume that `v` is $n \times 1$, it may be $n \times X$ for some $X > 1$ because of how MATLAB's `fminunc()` does things.
}

### Optimization Problem

The `Reduced_Space_Optimization` class combines an `Objective` and a `Constraint` to represent and solve an optimization problem of the form~{eq}`eqn:rs_opt_prob`.
Unlike the previous classes, the user does not need to subclass `Reduced_Space_Optimization`---it's ready to be used.
[Table TODO](tab:optimization_functions) lists the functions defined in this class.

\begin{table}[!ht]
\centering
\begin{tabular}{|l|l|}
    \hline
    Function Signature & Mathematical Description
    \\ \hline
    `Reduced_Space_Optimization(obj, con)`
    & Constructor taking an objective and constraints.
    \\
    `[u, z] = Optimize(z0)`
    & Solve the optimization problem with an initial control guess.
    \\ \hline
    `[val, grad, hessian_data] = Jhat(z)`
    & Use {prf:ref}`alg:adjoint_gradient` to compute $\hat{J}(z)$ and its gradient.
    \\
    `[Hv] = Jhat_hessVec(hessian_data, v)`
    & Use {prf:ref}`alg:adjoint_hessvec` to compute the action of $\nabla_{z,z}\hat{J}(z)$.
    \\ \hline
    `[diffs] = Finite_Difference_Gradient_Check(z)`
    & Compare `Jhat()` to finite differences.
    \\
    `[diffs] = Finite_Difference_Hessian_Check(z)`
    & Compare `Jhat_hessVec()` to finite differences.
    \\ \hline
\end{tabular}
\caption{Functions of the `Reduced_Space_Optimization` class.}
\label{tab:optimization_functions}
\end{table}

We implement {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_hessvec` in `Constrained_Optimization`'s private functions `Jhat()` and `Jhat_hessVec()`, respectively,  and interface them with MATLAB's solver `fminunc()` in the public function `Optimize()`. In addition, `Constrained_Optimization` has the public functions \\ `Finite_Difference_Gradient_Check()` and `Finite_Difference_Hessian_Check()` to test the accuracy of our calculation of $\nabla_z \hat{J}(z)$ and $\nabla_{z,z} \hat{J}(z)v$.

## Differential Equations

We specialize the general optimization problem~{eq}`eqn:opt_prob` and consider the case where $c(u, z)=0$ corresponds to the discretization of a system of ordinary differential equations (ODEs)

$$
\begin{aligned}
&\frac{dy}{dt} = f(y,z,t) \\
&y(0) = h(z)
\end{aligned}
$$

where $y:[0,T] \to \mathbb{R}^m$, $T>0$, $z \in \mathbb{R}^n$, $f:\mathbb{R}^m \times \mathbb{R}^n \times [0,T]  \to \mathbb{R}^m$, $h:\mathbb{R}^n \to \mathbb{R}^m$, and $u$ collects $y(t)$ at all monitored times $t$ (we make this precise shortly).
The objective function is

$$
\begin{aligned}
J(u,z) = \int_0^{T} g(y(t),t) dt + R(z)
\end{aligned}
$$

where $g:\mathbb{R}^m \times \mathbb{R} \to \mathbb{R}$ and $R:\mathbb{R}^n \to \mathbb{R}$ are user-specified functions.
This scenario is encapsulated in the subclass `Constrained_ODE_Optimization`, derived from the base class `Constrained_Optimization`.
The subclass inputs the state dimension $m$, the final time $T$, and the number of time steps $N$, and requires the user to specify the objective function and ODE system through the pure virtual functions summarized in [Table TODO](tab:con_ode_opt_pure_virtual_funs_1).

\begin{table}[!ht]
\centering
\begin{tabular}{|l|l|}
    \hline
    Function Signature & Mathematical Description
    \\ \hline
    `[val,grad_y] = Time_Instance_Objective(y,t)`
    & Evaluate $g(y,t)$ and $\nabla_y g(y,t)$
    \\
    `[val,grad_z] = Regularization_Objective(z)`
    & Evaluate $R(z)$ and $\nabla_z R(z)$
    \\
    `[f,f_y,f_z] = Time_Instance_RHS(y,z,t)`
    & Evaluate $f(y,z,t)$, $f_y(y,z,t)$, and $f_z(y,z,t)$
    \\
    `[h, h_z] = Initial_Condition(z)`
    & Evaluate $h(z)$ and $h_z(z)$
    \\
    `[Mv] = Time_Instance_Objective_yy_Apply(v,y,t)`
    & Compute the product $\nabla_{y,y} g(y,t) v$
    \\
    `[Mv] = Regularization_Objective_zz_Apply(v,z)`
    & Compute the product $\nabla_{z,z} R(z) v$
    \\
    `[Mv] = Time_Instance_RHS_yy_Apply(v,y,z,t,lambda)`
    & Compute the product $\lambda^{\mathsf{T}} f_{y,y}(y,z,t) v$
    \\
    `[Mv]  = Time_Instance_RHS_yz_Apply(v,y,z,t,lambda)`
    & Compute the product $\lambda^{\mathsf{T}} f_{y,z}(y,z,t) v$
    \\
    `[Mv]  = Time_Instance_RHS_zy_Apply(v,y,z,t,lambda)`
    & Compute the product $\lambda^{\mathsf{T}} f_{z,y}(y,z,t) v$
    \\
    `[Mv]  = Time_Instance_RHS_zz_Apply(v,y,z,t,lambda)`
    & Compute the product $\lambda^{\mathsf{T}} f_{z,z}(y,z,t) v$
    \\
    `[Mv]  = Initial_Condition_zz_Apply(v,z,lambda)`
    & Compute the product $\lambda^{\mathsf{T}} h_{z,z}(z) v$
    \\ \hline
\end{tabular}
\caption{Pure virtual function in `Constrained_ODE_Optimization`.}
\label{tab:con_ode_opt_pure_virtual_funs_1}
\end{table}

In order to implement the pure virtual functions of `Constrained_Optimization` (the functions listed in [Table TODO](tab:con_opt_pure_virtual_funs_1) and [Table TODO](tab:con_opt_pure_virtual_funs_2)), we must fix a time discretization scheme and write a finite dimensional optimization problem in the form of~{eq}`eqn:opt_prob`. To this end, let $t_k=T\frac{k-1}{N-1}$, $k=1,2,\dots,N$ be equally spaced nodes in the time interval $[0,T]$ and $y_{i,k}$ denote the value of the $i^{th}$ component of $y(t_k)\in\mathbb{R}^{m}$, i.e. the $i^{th}$ state variable (of $m$ variables) at the $k^{th}$ time (of $N$ nodes). Concatenating these discrete evaluations yields the space-time state vector

$$
\begin{aligned}
    u
    &= \left(\begin{array}{c}
        y_1 \\ y_2 \\ \vdots \\ y_N
    \end{array}\right)
    \in \mathbb{R}^{mN},
    &&\text{where}&
    y_k
    &= \left(\begin{array}{c}
        y_{1,k} \\ y_{2,k} \\ \vdots \\ y_{m,k}
    \end{array}\right) \in \mathbb{R}^{m}.
\end{aligned}
$$

The full state $u\in\mathbb{R}^{mN}$ therefore consists of the values of the discrete ODE state $y(t)\in\mathbb{R}^{m}$ at all $N$ times $t_{1},\ldots,t_{N}$.
Thus, in this case we have $n_{u} = mN$ and $n_{z} = n$.
Next, applying the backward Euler time integration scheme, we have

$$
\begin{aligned}
    c(u,z) =
    \left(\begin{array}{c}
        y_1-h(z) \\
        y_2 - y_1 - \Delta t f(y_2,z,t_2) \\
        y_3 - y_2 - \Delta t f(y_3,z,t_3) \\
        y_4 - y_3 - \Delta t f(y_4,z,t_4) \\
        \vdots \\
        y_N - y_{N-1} - \Delta t f(y_N,z,t_N) \\
    \end{array}\right) \in \mathbb{R}^{mN},

\end{aligned}
$$ (eqn:back_euler)

where $\Delta t = \frac{T}{N-1}$. Hence, $n_{c} = n_{u} = mN$. Finally, applying the trapezoid rule for integration, the discretized objective function is

$$
\begin{aligned}
    J(u,z) = \sum\limits_{k=1}^N w_k g(y_k,t_k) + R(z),

\end{aligned}
$$ (eqn:dis_trap_objective)

where $w_k=\Delta T$, $k=2,3,\dots,N-1$, and $w_1=w_N=0.5 \Delta T$, are integration weights.

Given an implementation of $g$, $R$, and their derivatives, the pure virtual functions corresponding to the objective $J$ and its derivatives are easily implemented based on~{eq}`eqn:dis_trap_objective`. The pure virtual functions corresponding to the constraint $c$ require more careful consideration.
To solve the equation $c(u,z)=0$, we exploit the block structure in~{eq}`eqn:back_euler` and solve $m \times m$ systems of equations at each of the $N$ time points sequentially (i.e., perform time integration rather than solving $c(u,z)=0$ as a system of $mN$ coupled equations). However, the structure of $c(u,z)$ is critical for ensuring that the adjoints are computed correctly. The Jacobian of $c$ with respect to $u$, $c_u(u,z) $, is the $mN \times mN$ matrix

$$
\begin{aligned}
    \left( \begin{array}{cccccc}
    I_m & 0 & 0 & 0 & \dots & 0 \\
    -I_m & I_m - \Delta t f_y(y_2,z,t_2) & 0 & 0 & \dots & 0 \\
    0 & -I_m & I_m - \Delta t f_y(y_3,z,t_3) & 0 & \dots & 0 \\
    0 & 0 & -I_m & I_m - \Delta t f_y(y_4,z,t_4) & \dots & 0 \\
    \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
    0 & 0  & \dots & 0 & -I_m & I_m - \Delta t f_y(y_N,z,t_N) \\
    \end{array} \right).
\end{aligned}
$$

Hence the implementation of $c_u(u,z)^{-\mathsf{T}} v$ for the adjoint solve in {prf:ref}`alg:adjoint_gradient` must consider the system of linear equations $c_u(u,z)^{\mathsf{T}} x = v$ written as

$$
\begin{aligned}
    \left( \begin{array}{cccccc}
    I_m & -I_m & 0 & 0 & \dots & 0 \\
    0 & I_m-\Delta t f_y(y_2,z,t_2) & -I_m & 0 & \dots & 0 \\
    0 & 0 & I_m - \Delta t f_y(y_3,z,t_3) & -I_m & \dots & 0 \\
    0 & 0 & 0 & I_m - \Delta t f_y(y_4,z,t_4) & -I_m & \dots  \\
    \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
    0 & 0  & \dots & 0 & 0 & I_m - \Delta t f_y(y_N,z,t_N) \\
    \end{array} \right)
    \left( \begin{array}{c}
    x_1 \\
    x_2 \\
    x_3 \\
    x_4 \\
    \vdots \\
    x_N \\
    \end{array} \right) \\
    =
    \left( \begin{array}{c}
    v_1 \\
    v_2 \\
    v_3 \\
    v_4 \\
    \vdots \\
    v_N \\
    \end{array} \right)
\end{aligned}
$$

and compute $x=c_u(u,z)^{-\mathsf{T}} v$ via the backward-in-time stepping algorithm

$$
\begin{aligned}
    & x_N = (I_m - \Delta t f_y(y_N,z,t_N))^{-1} v_N \\
    & x_k = (I_m - \Delta t f_y(y_k,z,t_k))^{-1} (v_k + x_{k+1}) \qquad k=N-1,n-2,\dots,2 \\
    & x_1 = v_1 + x_2 .
\end{aligned}
$$

This ensures that the time discretization of the adjoint equation is consistent with the time integration of the state equation at the discrete level, thus providing an accurate gradient computation. Similar care must be taken when implementing all other matrix-vector products and linear solves involving derivatives of the constraint $c$.
