# Constrained Optimization

:::{admonition} Summary
:class: dropdown

SABL defines the following MATLAB classes.

- {class}`Objective`, representing $J(\u,\z)$.
- {class}`Constraint`, representing $\c(\u,\z)$ and $\S:\u\mapsto\z$ such that $\c(\S(\z),\z) = \0$.
- {class}`Reduced_Space_Optimization`, representing $\min_{z}J(\S(\z),\z)$.

The user must define subclasses of {class}`Objective` and {class}`Constraint`, implement their abstract methods, and use {meth}`Reduced_Space_Optimization.Optimize()` to solve an optimization problem.

In addition, SABL defines specialized classes {class}`Dynamic_Objective` and {class}`Dynamic_Constraint` for representing a wide class of optimization problems where the state and/or control are time dependent.
:::

This page describes tools in SABL for solving [constrained optimization problems](../../problems/optimization) of the form

$$
\begin{aligned}
    \min_{\u,\z} ~& J(\u,\z)
    \qquad\Longleftrightarrow\qquad
    \min_{\z}\hat{J}(\z) = J(\S(\z),\z),
    \\
    s.t. ~~& \c(\u,\z) = \0,
\end{aligned}
$$ (sabl:opt_prob)

where $\u\in\R^{n_u}$ is the state variable, $\z\in\R^{n_z}$ is the control variable, $J(\u,\z)$ is the objective function, $\c(\u,\z)$ specifies the constraints, $\S(\z)$ is the solution operator corresponding to the constraints (i.e., $\c(\S(\z),\z) = \0$ for all $\z$), and $\hat{J}(\z)$ is the objective function in reduced space.

## Step 1: Extend Abstract Classes

::::{margin}
:::{admonition} Abstract Classes
:class: note

An _abstract class_ is a class with _abstract methods_, functions that must be implemented before the class can be instantiated.
Abstract classes serve as a template for classes that inherit from them.
See [MATLAB's page on abstract classes](https://www.mathworks.com/help/matlab/matlab_oop/abstract-classes-and-interfaces.html) for details.
:::
::::

To implement {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_hessvec` and solve {eq}`sabl:opt_prob`, SABL adopts an object-oriented design in MATLAB by defining the following abstract classes.

- {class}`Objective` represents an optimization objective $J(\u,\z)$.
- {class}`Constraint` represents constraints $\c(\u,\z)$ and the corresponding solution operator $\S(\z)$.

Each of these classes defines abstract methods that must be implemented in user-defined subclasses to define specific instances of {eq}`sabl:opt_prob`.
For instance, {meth}`Objective.J()` computes the objective value $J(\u,\z)$ and its gradients $\grad{u}J(\u,\z)$ and $\grad{z}J(\u,\z)$ for use in {prf:ref}`alg:adjoint_gradient`; {meth}`Constraint.c_u_Inverse_Apply()` solves $\c_u(\u,\z)\bfmu = \v$ for $\bfmu$ to facilitate step 6 of {prf:ref}`alg:adjoint_hessvec`.
See the [documentation](../api) on each class for a list of abstract methods and templates for writing new inherited classes.

:::{admonition} Time-dependent Problems
:class: important

SABL defines special cases of {class}`Objective` and {class}`Constraint` for [a class of constrained optimization problems](optimization:differential) where the state is a function of time.

- {class}`Dynamic_Objective` represents the objective function $J(\u,\z) = \int_{0}^{T} g(\u(t), t) dt + R(\z)$.
- {class}`Dynamic_Constraint` represents constraints expressed by ordinary differential equations, $\frac{\textrm{d}}{\textup{d}t}\y(t) = \f(\y(t),\z,t)$.

Each of these classes implements the abstract method defined by their parent classes and defines new abstract methods.
For example, the value and gradients of $J$ can be written in terms of $g$ and $R$, so {class}`Dynamic_Objective` implements {meth}`Objective.J()` and instead requires the user to implement methods to evaluate $g$, $R$, and their gradients.
:::

## Step 2: Verify Implementation

Most SABL classes are equipped with methods for checking the implementation of the abstract methods with finite differences.
These sanity checks do not guarantee that an implementation is correct, but they are good indicators for finding mistakes.

- {meth}`Objective.Finite_Difference_Gradient_Check()`
- {meth}`Objective.Finite_Difference_Hessian_Check()`
- {meth}`Constraint.Finite_Difference_Constraint_Check()`
- {meth}`Dynamic_Constraint.Time_Instance_RHS_Jacobian_Check()`
- {meth}`Dynamic_Constraint.Time_Instance_RHS_Hessian_Check()`
- {meth}`Reduced_Space_Optimization.Finite_Difference_Gradient_Check()`
- {meth}`Reduced_Space_Optimization.Finite_Difference_Hessian_Check()`

:::{warning}
{meth}`Constraint.Finite_Difference_Constraint_Check()` requires the optional method {meth}`Constraint.c()` to be implemented.
:::

## Step 3: Solve the Problem

The {class}`Reduced_Space_Optimization` class combines an {class}`Objective` and a {class}`Constraint` to represent and solve an optimization problem of the form {eq}`sabl:opt_prob`.
Specifically, {meth}`Reduced_Space_Optimization.Jhat()` and {meth}`Reduced_Space_Optimization.Jhat_hessVec()` utilize the abstract methods of the objective and the constraints to implement {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_hessvec`/{prf:ref}`alg:adjoint_gaussnewton`.
Unlike the previous classes, the user does not need to subclass {meth}`Reduced_Space_Optimization`---it is ready to be used as is.

:::{admonition} Try It Out!
:class: seealso

[Tutorial 1](../tutorials/example1) shows this workflow in action for a simple problem.
:::
