# Constrained Optimization

:::{admonition} Summary
:class: dropdown

SABL defines the following MATLAB classes.

- [`Objective`](sabl:optimization-objective), representing $J(u,z)$.
- [`Constraint`](sabl:optimization-constraints), representing $c(u,z)$ and $S:u\mapsto z$ such that $c(S(z),z) = 0$.
- [`Reduced_Space_Optimization`](sabl:optimizer-class), representing $\min_{z}J(S(z),z)$.

The user must define subclasses of `Objective` and `Constraint`, implement their abstract methods, and use `Reduced_Space_Optimization.Optimize()` to solve an optimization problem.

In addition, SABL defines specialized classes [`Dynamic_Objective`](sabl:optimization-dynamicobjective) and [`DynamicConstraint`](sabl:optimization-dynamicconstraint) for representing a wide class of optimization problems where the state and/or control are time dependent.
:::

This page describes tools in SABL for solving [constrained optimization problems](../../problems/optimization) of the form

$$
\begin{aligned}
    \min_{u,z} ~& J(u,z)
    \qquad\Longleftrightarrow\qquad\min_{z}\hat{J}(z)=J(S(z),z),
    \\
    s.t. ~~& c(u,z) = 0,
\end{aligned}
$$ (sabl:opt_prob)

::::{margin}
:::{admonition} Abstract Classes
:class: note

An _abstract class_ is a class with _abstract methods_, functions that must be implemented before the class can be instantiated.
Abstract classes serve as a template for classes that inherit from them.
See [MATLAB's page on abstract classes](https://www.mathworks.com/help/matlab/matlab_oop/abstract-classes-and-interfaces.html) for details.
:::
::::

where $u$ is the state variable, $z$ is the control variable, $J$ is the objective function, $c$ specifies the constraints, $S(z)$ is the solution operator corresponding to the constraints (i.e., $c(S(z),z) = 0$ for all $z$), and $\hat{J}$ is the objective function in reduced space.

To implement {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_hessvec` and solve {eq}`sabl:opt_prob`, SABL adopts an object-oriented design in MATLAB by defining abstract classes [`Objective`](sabl:optimization-objective) for the optimization objective $J(u, z)$ and [`Constraint`](sabl:optimization-constraints) for the constraints $c(u, z)$.
Each of these classes defines abstract methods that must be implemented in user-defined subclasses to solve specific instances of {eq}`sabl:opt_prob`.
The [`Reduced_Space_Optimization`](sabl:optimizer-class) class couples `Objective` and `Constraint` objects and solves the corresponding constrained optimization problem.

(sabl:optimization-objective)=
## Optimization Objective

The abstract `Objective` class is a template for representing an objective functional $J(u, z)$ and its derivatives.

### Abstract Objective Methods

To represent a particular objective function with the `Objective` class, define a new class that inherits from `Objective` and implements the abstract methods listed in [Table 1](tab:objective_abstract).

:::{table} Abstract methods of the `Objective` class.
:align: center
:name: tab:objective_abstract

| Function Signature                | Mathematical Description         |
| :-------------------------------- | :------------------------------- |
| `[val, grad_u, grad_z] = J(u, z)` | Evaluate $J(u,z)$, $\nabla_{\!u} J(u,z)$, and $\nabla_{\!z} J(u,z)$ |
| `[Mv] = J_uu_Apply(v, u, z)`      | Compute $\nabla_{\!u,u} J(u, z) v$ |
| `[Mv] = J_uz_Apply(v, u, z)`      | Compute $\nabla_{\!u,z} J(u, z) v$ |
| `[Mv] = J_zu_Apply(v, u, z)`      | Compute $\nabla_{\!z,u} J(u, z) v$ |
| `[Mv] = J_zz_Apply(v, u, z)`      | Compute $\nabla_{\!z,z} J(u, z) v$ |
:::

:::{danger}
Because of how MATLAB's [`fminunc()`](https://www.mathworks.com/help/optim/ug/fminunc.html) is designed, the `J_xx_Apply()` functions must be implemented in a _vectorized_ fashion, i.e., assuming that `v` is a matrix where each column is a test direction.
This is demonstrated in the examples.
:::

(sabl:optimization-dynamicobjective)=
### Time-Dependent Objectives

The `Dynamic_Objective` class inherits from `Objective` and can be used in the special case in which the optimization state $u$ is a function of time and the objective functional can be written as

$$
\begin{align*}
    J(u, z)
    = \int_{0}^{T}g(u(t),t)dt + R(z).
\end{align*}
$$

The integral in the objective is discretized using the trapezoidal rule, i.e.,

$$
\begin{align*}
    J(\mathbf{u},\mathbf{z})
    = \sum_{j=1}^{N}w_{j}g(\mathbf{y}_{j}, t_{j}) + R(\mathbf{z}),
    \qquad
    \mathbf{u}
    = \left(\begin{array}{c}
        \mathbf{y}_{1} \\ \vdots \\ \mathbf{y}_{N}
    \end{array}\right),
\end{align*}
$$

where $w_{j} = \delta t = \frac{T}{N - 1}$ for $j=2,\ldots,N-1$ and $w_{1} = w_{N} = \frac{1}{2}\delta t$, with $N$ being the number of points in the temporal discretization.
Here, $\mathbf{y}_{j}$ is the state at time $t = t_{j}$.

Because the derivatives of $J$ can be described in terms of the functions $g$ and $R$, specific instances of `Dynamic_Objective` must implement methods describing $g$, $R$, and their derivatives instead of the methods in [Table 1](tab:objective_abstract).

:::{table} Abstract methods of the `Dynamic_Objective` class.
:align: center
:name: tab:dynamicobjective_abstract

| Function Signature                                 | Mathematical Description                        |
| :------------------------------------------------- | :---------------------------------------------- |
| `[val, grad_y] = Time_Instance_Objective(y, t)`    | Evaluate $g(y(t), t)$ and $\nabla_{\!y} g(y(t), t)$ |
| `[val, grad_z] = Regularization_Objective(z)`      | Evaluate $R(z)$ and $\nabla_{\!z} R(z)$             |
| `[Mv] = Time_Instance_Objective_yy_Apply(v, y, t)` | Compute the product $\nabla_{\!y,y} g(y(t), t) v$ |
| `[Mv] = Regularization_Objective_zz_Apply(v, z)`   | Compute the product $\nabla_{\!z,z} R(z) v$       |
:::

### Objective Constructor

The `Objective` class takes no arguments in the constructor.
However, the `Dynamic_Objective` class needs user-provided information to define the discretization of the time integral. [Table 3](tab:dynamicobjective_constructor) lists the arguments for the constructor.

:::{table} Constructor arguments of the `Dynamic_Objective` class.
:align: center
:name: tab:dynamicobjective_constructor

| Argument | Description                                                     |
| :------- | :-------------------------------------------------------------- |
| `m`      | Dimension of the state variable at each time step               |
| `n`      | Dimension of the control variable at each time step             |
| `T`      | Final time $T$, the upper limit of the integral                 |
| `N`      | Number of points in the temporal discretization of the integral |
:::

### Objective Verification

The `Objective` class (and any class that inherits from it) is equipped with the following methods to verify the consistency between the objective function itself and its gradient and Hessian functions.

:::{table} Verification functions in the `Objective` class.
:align: center
:name: tab:objective_checkers

| Function Signature                       | Mathematical Description             |
| :--------------------------------------- | :----------------------------------- |
| `Finite_Difference_Gradient_Check(u, z)` | Check gradients of $J$               |
| `Finite_Difference_Hessian_Check(u, z)`  | Check Hessian-vector products of $J$ |
:::

(sabl:optimization-template)=
### Objective Templates

The following template can be used to start a new `Objective` subclass.

```matlab
classdef My_Objective < Objective

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

The following template can be used for a new `DynamicObjective` subclass.

```matlab
classdef My_Dynamic_Objective < Dynamic_Objective

    methods (Access = public)

        function [val, grad_y] = Time_Instance_Objective(this, y, t)
            error('Time_Instance_Objective() not implemented');
        end

        function [val, grad_z] = Regularization_Objective(this, z)
            error('Regularization_Objective() not implemented');
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this, v, y, t)
            error('Time_Instance_Objective_yy_Apply() not implemented');
        end

        function [Mv] = Regularization_Objective_zz_Apply(this, v, z)
            error('Regularization_Objective_zz_Apply() not implemented');
        end

    end
end
```

(sabl:optimization-constraints)=
## Optimization Constraints

The `Constraint` class encodes the constraint function $c(u, z) = 0$, its derivatives, and the solution operator $S:z\mapsto u$ satisfying $c(S(z), z) = 0$ for all $z$.

### Abstract Constraint Methods

To represent a particular set of constraints with the `Constraint` class, define a new class that inherits from `Constraint` and implements the abstract methods listed in [Table 5](tab:constraint_abstract).

:::{table} Abstract methods of the `Constraint` class.
:align: center
:name: tab:constraint_abstract

| Function Signature                    | Mathematical Description             |
| :------------------------------------ | :----------------------------------- |
| `[u] = State_Solve(z)`                | Given $z$, solve $c(u, z)=0$ for $u$ |
| `[Mv] = c_u_Transpose_Inverse_Apply(v, u, z)` | Compute $c_u(u, z)^{-\mathsf{T}} v$ |
| `[Mv] = c_z_Transpose_Apply(v, u, z)` | Compute $c_z(u, z)\trp v$   |
| `[Mv] = c_u_Inverse_Apply(v, u, z)`   | Compute $c_u(u, z)^{-1} v$           |
| `[Mv] = c_z_Apply(v, u, z)`           | Compute $c_z(u, z) v$                |
:::

:::{danger}
Because of how MATLAB's [`fminunc()`](https://www.mathworks.com/help/optim/ug/fminunc.html) is designed, the `c_x_XXX()` methods (e.g., `c_z_Apply()`) must be implemented in a _vectorized_ fashion, i.e., assuming that `v` is a matrix where each column is a test direction.
:::

The following methods are also required **unless** `Gauss_Newton_Hess` is set to `true` in the `Reduced_Space_Optimization`, which uses {prf:ref}`alg:adjoint_gaussnewton` in place of {prf:ref}`alg:adjoint_hessvec`.

:::{table} Methods of the `Constraint` class required for Gauss--Newton minimization.
:align: center
:name: tab:constraint_fullhessian

| Function Signature                   | Mathematical Description                        |
| :----------------------------------- | :---------------------------------------------- |
| `[Mv] = c_uu_Apply(v, u, z, lambda)` | Compute $\lambda\trp c_{u, u}(u, z) v$ |
| `[Mv] = c_uz_Apply(v, u, z, lambda)` | Compute $\lambda\trp c_{u, z}(u, z) v$ |
| `[Mv] = c_zu_Apply(v, u, z, lambda)` | Compute $\lambda\trp c_{z, u}(u, z) v$ |
| `[Mv] = c_zz_Apply(v, u, z, lambda)` | Compute $\lambda\trp c_{z, z}(u, z) v$ |
:::

(sabl:optimization-dynamicconstraint)=
### Time-Dependent Differential Constraints

The `Dynamic_Constraint` class inherits from `Constraint` and is designed to encode constraints represented as a system of ordinary differential equations (possibly a spatially discretized partial differential equation), i.e.,

$$
\begin{align*}
    \frac{\textup{d}}{\textup{d}t}y(t) = f(y(t), z, t),
    \qquad
    y(0) = h(z).
\end{align*}
$$

Because the derivatives of the constraint can be described in terms of the functions $f$ and $h$, specific instances of `Dynamic_Objective` must implement methods describing $f$, $h$, and their derivatives instead of the methods in [Table 5](tab:constraint_abstract) and [Table 6](tab:constraint_fullhessian).

:::{table} Abstract methods of the `DynamicConstraint` class.
:align: center
:name: tab:dynamicconstraint_abstract

| Function Signature                               | Mathematical Description             |
| :----------------------------------------------- | :----------------------------------- |
| `[f, f_y, f_z] = Time_Instance_RHS(y, z, t)`     | Evaluate $f(y(t),z,t)$, $\nabla_{\!y}f$, and $\nabla_{\!z}f$ |
| `[h, h_z] = Initial_Condition(z)`                | Evaluate $h(z)$ and $\nabla_{\!z}h(z)$ |
| `Time_Instance_RHS_yy_Apply(v, y, z, t, lambda)` | Compute $\lambda\trp \nabla_{\!yy}f(y(t), z) v$ |
| `Time_Instance_RHS_yz_Apply(v, y, z, t, lambda)` | Compute $\lambda\trp \nabla_{\!yz}f(y(t), z) v$ |
| `Time_Instance_RHS_zy_Apply(v, y, z, t, lambda)` | Compute $\lambda\trp \nabla_{\!zy}f(y(t), z) v$ |
| `Time_Instance_RHS_zz_Apply(v, y, z, t, lambda)` | Compute $\lambda\trp \nabla_{\!zz}f(y(t), z) v$ |
| `Initial_Condition_zz_Apply(v, z, lambda)`       | Compute $\lambda\trp\nabla_{\!z}h(z) v$ |
:::

### Constraint Constructor

The `Constraint` class takes no arguments in the constructor.
However, the `Dynamic_Constraint` class needs user-provided information to define the time integration scheme. [Table 8](tab:dynamicconstraint_constructor) lists the arguments for the constructor.

:::{table} Constructor arguments of the `Dynamic_Objective` class.
:align: center
:name: tab:dynamicconstraint_constructor

| Argument | Description                                                     |
| :------- | :-------------------------------------------------------------- |
| `m`      | Dimension of the state variable at each time step               |
| `n`      | Dimension of the control variable at each time step             |
| `T`      | Final time $T$, the upper limit of the integral                 |
| `N`      | Number of points in the temporal discretization of the integral |
:::

### Constraint Verification

The `Constraint` class (and any class that inherits from it) is equipped with the following methods to verify the consistency between the objective function itself and its gradient and Hessian functions.
However, these methods require implementing the method `[c] = c(this, u, z)` to explicitly construct the constraints.

:::{table} Verification functions in the `Constraint` class.
:align: center
:name: tab:constraint_checkers

| Function Signature                       | Mathematical Description             |
| :--------------------------------------- | :----------------------------------- |
| `Finite_Difference_Constraint_Check(u, z)` | Check gradients of $c$             |
:::

:::{table} Verification functions in the `Dynamic_Constraint` class.
:align: left
:name: tab:dynamicconstraint_checkers

| Function Signature                            | Mathematical Description |
| :-------------------------------------------- | :----------------------- |
| `Time_Instance_RHS_Jacobian_y_Check(y, z, t)` | Check $\nabla_{\!y}f$    |
| `Time_Instance_RHS_Jacobian_z_Check(y, z, t)` | Check $\nabla_{\!z}f$    |
| `Time_Instance_RHS_Hessian_yy_Check(y, z, t)` | Check $\nabla_{\!yy}f$   |
| `Time_Instance_RHS_Hessian_yz_Check(y, z, t)` | Check $\nabla_{\!yy}f$   |
| `Time_Instance_RHS_Hessian_zy_Check(y, z, t)` | Check $\nabla_{\!yy}f$   |
| `Time_Instance_RHS_Hessian_zz_Check(y, z, t)` | Check $\nabla_{\!yy}f$   |
:::

(sabl:constraint-template)=
### Constraint Templates

The following template can be used to start a new `Constraint` subclass.

```matlab
classdef My_Constraint < Constraint

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

        % The following methods are used if Gauss_Newton_Hess=false
        % in the Reduced_Space_Optimization (this is the default).
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

The following template can be used for a new `DynamicConstraint` subclass.

```matlab
classdef My_Dynamic_Constraint < Dynamic_Constraint

    methods (Abstract, Access = public)

        function [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
            error('Time_Instance_RHS() not implemented');
        end

        function [h, h_z] = Initial_Condition(this, z)
            error('Initial_ConditionApply() not implemented');
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)
            error('Time_Instance_RHS_yy_Apply() not implemented');
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
            error('Time_Instance_RHS_yz_Apply() not implemented');
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
            error('Time_Instance_RHS_zy_Apply() not implemented');
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
            error('Time_Instance_RHS_zz_Apply() not implemented');
        end

        function [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)
            error('Initial_Condition_zz_Apply() not implemented');
        end

    end
end
```

(sabl:optimizer-class)=
## Optimizer Class

The `Reduced_Space_Optimization` class combines an `Objective` and a `Constraint` to represent and solve an optimization problem of the form {eq}`sabl:opt_prob`.

### Solving the Optimization Problem

The private methods `Jhat()` and `Jhat_hessVec()` use the methods of the objective and constraint to implement {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_hessvec`/{prf:ref}`alg:adjoint_gaussnewton`.
Unlike the previous classes, the user does not need to subclass `Reduced_Space_Optimization`---it is ready to be used as is.

:::{table} Public methods of the `Reduced_Space_Optimization` class.
:align: center
:name: tab:optimization_methods

| Function Signature                     | Mathematical Description                        |
| :------------------------------------- | :---------------------------------------------- |
| `Reduced_Space_Optimization(obj, con)` | Constructor taking an objective and constraints |
| `[u, z] = Optimize(z0)`   | Solve the optimization problem with an initial control guess |
| `Finite_Difference_Gradient_Check(z)`  | Compare `Jhat()` to finite differences          |
| `Finite_Difference_Hessian_Check(z)`   | Compare `Jhat_hessVec()` to finite differences  |
:::

### Customizing the Optimization

After instantiating a `Reduced_Space_Optimization` object, the optimization routine can be customized by setting the object attributes.
The options are listed in the following table, most of which correspond to an option for MATLAB's [`fminunc()`](https://www.mathworks.com/help/optim/ug/fminunc.html) or [`fmincon()`](https://www.mathworks.com/help/optim/ug/fmincon.html).

:::{table} Attributes of the `Reduced_Space_Optimization` class for customizing optimization.
:align: center
:name: tab:optimization_attributes

| Attribute           | Description                                       | `fminunc()` argument  | Default value |
| :------------------ | :------------------------------------------------ | :-------------------- | :------------ |
| `Gauss_Newton_Hess` | Use the Gauss-Newton approximation of the Hessian ({prf:ref}`alg:adjoint_gaussnewton`) | | `false`       |
| `use_trust_region`  | Use trust region for the optimization             |                       | `true`        |
| `z_lb`              | Lower bounds for control                          |                       | `[]`          |
| `z_ub`              | Upper bounds for control                          |                       | `[]`          |
| `opt_tol`           | Optimality tolerance                              | `OptimalityTolerance` | `10^-8`       |
| `fun_tol`           | Function tolerance                                | `FunctionTolerance`   | `10^-6`       |
| `iteration_limit`   | Maximum number of iterations                      | `MaxIterations`       | `1000`        |
| `step_tol`          | Step tolerance                                    | `StepTolerance`       | `10^-6`       |
| `max_cg_iter`       | Maximum number of conjugate gradient iterations   | `MaxPCGIter`          | `50`          |
| `cg_tol`            | Conjugate gradient tolerance                      | `TolPCG`              | `10^-4`       |
| `verbose`           | Print optimization info at each iteration         | `Display`             | `true` (`'iter-detailed'`) |
:::

The following code is an example of how to do this.

```matlab
opt = Reduced_Space_Optimization(obj, con);
opt.iteration_limit = 200;
opt.Gauss_Newton_Hess = true;
opt.use_trust_region = false;
opt.verbose = false;
[u_solution, z_solution] = opt.Optimize(z0);
```
