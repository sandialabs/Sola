# Code Structure

This page describes the software design paradigm of SABL.

## Constrained Optimization

:::{admonition} Summary
:class: dropdown

SABL defines the following MATLAB classes.

- [`Objective`](sabl:optimization-objective), representing $J(u,z)$.
- [`Constraint`](sabl:optimization-constraints), representing $c(u,z)$.
- [`Reduced_Space_Optimization`](sabl:optimizer-class), representing $\min_{u,z}J(u,z)$ subject to $c(u,z) = 0$.

The user must define subclasses of `Objective` and `Constraint`, implement their abstract methods, and use `Reduced_Space_Optimization.Optimize()` to solve an optimization problem.
<!-- For time-dependent problems, the specialized  -->
:::

This section focuses on solving constrained optimization problems of the form

$$
\begin{aligned}
    \min_{u,z} ~& J(u,z)
    \qquad\Longleftrightarrow\qquad\min_{z}\hat{J}(z)=J(S(z),z),
    \\
    s.t. ~~& c(u,z) = 0,
\end{aligned}
$$ (sabl:opt_prob)

where $u$ is the state variable, $z$ is the control variable, $J$ is the objective function, $c$ specifies the constraints, $S(z)$ is the solution operator corresponding to the constraints (i.e., $c(S(z),z) = 0$ for all $z$), and $\hat{J}$ is the objective function in reduced space.

To implement {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_hessvec` and solve {eq}`sabl:opt_prob`, SABL adopts an object-oriented design in MATLAB by defining [abstract classes](https://www.mathworks.com/help/matlab/matlab_oop/abstract-classes-and-interfaces.html) [`Objective`](sabl:optimization-objective) for the optimization objective $J$ and [`Constraint`](sabl:optimization-constraints) for the constraints $c$.
Each of these classes defines abstract methods that must be implemented in user-defined subclasses to solve specific instances of {eq}`sabl:opt_prob`.
The [`Reduced_Space_Optimization`](sabl:optimizer-class) class couples `Objective` and `Constraint` objects and solves the corresponding constrained optimization problem.

(sabl:optimization-objective)=
### Optimization Objective

The `Objective` class is designed to encode the objective functional $J$ and its derivatives.
To represent a particular objective function with the `Objective` class, define a new class that inherits from `Objective` and implements the pure functions listed below in [Table 1](tab:objective_abstract).

:::{table} Abstract methods of the `Objective` class.
:align: center
:name: tab:objective_abstract

| Function Signature | Mathematical Description |
| :----------------- | :----------------------- |
| `[val, grad_u, grad_z] = J(u, z)` | Evaluate $J(u,z)$, $\nabla_u J(u,z)$, and $\nabla_z J(u,z)$ |
| `[Mv] = J_uu_Apply(v, u, z)` | Compute the product $\nabla_{u,u} J(u, z) v$ |
| `[Mv] = J_uz_Apply(v, u, z)` | Compute the product $\nabla_{u,z} J(u, z) v$ |
| `[Mv] = J_zu_Apply(v, u, z)` | Compute the product $\nabla_{z,u} J(u, z) v$ |
| `[Mv] = J_zz_Apply(v, u, z)` | Compute the product $\nabla_{z,z} J(u, z) v$ |
:::

:::{danger}
Because of how MATLAB's `fminunc()` function is designed, the `J_xx_Apply()` functions must be implemented in a _vectorized_ fashion, i.e., assuming that `v` is a matrix where each column is a test direction.
:::

#### Objective Verification

In addition to the functions listed in [Table 1](tab:objective_abstract), the `Objective` class is equipped with the following methods to verify the consistency between the objective function itself and its gradient and Hessian functions.

:::{table} Verification functions in the `Objective` class.
:align: center
:name: tab:objective_checkers

| Function Signature | Mathematical Description |
| :----------------- | :----------------------- |
| `Finite_Difference_Gradient_Check(u, z)` | Check gradients of $J$ |
| `Finite_Difference_Hessian_Check(u, z)` | Check Hessian-vector products of $J$ |
:::

(sabl:optimization-template)=
#### Objective Template

The following template can be used to start a new `Objective` subclass.

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

(sabl:optimization-constraints)=
### Optimization Constraints

:::{attention}
The rest of this page is under construction, please check back later.
:::

The `Constraint` class encodes the constraint function $c(u, z) = 0$, its derivatives, and the solution operator $S:z\mapsto u$.
To represent a particular set of constraints with the `Constraint` class, define a new class that inherits from `Constraint` and implements the abstract methods listed in [Table 3](tab:constraint_abstract).

:::{table} Abstract methods of the `Constraint` class.
:align: center
:name: tab:constraint_abstract

| Function Signature                    | Mathematical Description |
| :------------------------------------ | :----------------------- |
| `[u] = State_Solve(z)`                | Given $z$, solve $c(u, z)=0$ for $u$ |
| `[Mv] = c_u_Transpose_Inverse_Apply(v, u, z)` | Compute the product $c_u(u, z)^{-\mathsf{T}} v$ |
| `[Mv] = c_z_Transpose_Apply(v, u, z)` | Compute the product $c_z(u, z)^{\mathsf{T}} v$ |
| `[Mv] = c_u_Inverse_Apply(v, u, z)`   | Compute the product $c_u(u, z)^{-1} v$ |
| `[Mv] = c_z_Apply(v, u, z)`           | Compute the product $c_z(u, z) v$ |
| `[Mv] = c_uu_Apply(v, u, z, lambda)`  | Compute the product $\lambda^{\mathsf{T}} c_{u, u}(u, z) v$ |
| `[Mv] = c_uz_Apply(v, u, z, lambda)`  | Compute the product $\lambda^{\mathsf{T}} c_{u, z}(u, z) v$ |
| `[Mv] = c_zu_Apply(v, u, z, lambda)`  | Compute the product $\lambda^{\mathsf{T}} c_{z, u}(u, z) v$ |
| `[Mv] = c_zz_Apply(v, u, z, lambda)`  | Compute the product $\lambda^{\mathsf{T}} c_{z, z}(u, z) v$ |
:::

:::{danger}
Because of how MATLAB's `fminunc()` function is designed, the `c_x_XXX()` methods (e.g., `c_z_Apply()`) must be implemented in a _vectorized_ fashion, i.e., assuming that `v` is a matrix where each column is a test direction.
The `c_xx_Apply()` methods may assume that `v` is a column vector.
:::

(sabl:constraint-template)=
#### Constraint Template

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

(sabl:optimizer-class)=
### Optimizer Class

The `Reduced_Space_Optimization` class combines an `Objective` and a `Constraint` to represent and solve an optimization problem of the form {eq}`sabl:opt_prob`.
The private methods `Jhat()` and `Jhat_hessVec()` use the methods of the objective and constraint to implement {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_hessvec`.
Unlike the previous classes, the user does not need to subclass `Reduced_Space_Optimization`---it is ready to be used as is.

:::{table} Methods of the `Reduced_Space_Optimization` class.
:align: center
:name: tab:optimization_methods

| Function Signature                    | Mathematical Description |
| :------------------------------------ | :----------------------- |
| `Reduced_Space_Optimization(obj, con)` | Constructor taking an objective and constraints |
| `[u, z] = Optimize(z0)` | Solve the optimization problem with an initial control guess |
| `Finite_Difference_Gradient_Check(z)` | Compare `Jhat()` to finite differences |
| `Finite_Difference_Hessian_Check(z)` | Compare `Jhat_hessVec()` to finite differences |
:::

### Time-Dependent Objective

:::{warning}
The rest of this page is under construction, please check back later.
:::

### Differential Constraints

## Bayesian Inversion
