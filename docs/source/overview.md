# Project Overview

WOLF is the **W**aanders **O**uter-**L**oop **F**ramework for iteratively solving large-scale problems such as optimization problems constrained by partial differential equations.
This project includes documentation and tutorials for two separate packages:

1. The **S**andbox for **A**djoint-**B**ased outer **L**oop analysis ([SABL](./sabl/about)), a [MATLAB](https://www.mathworks.com/products/matlab.html) framework designed for prototyping outer-loop problems.
2. The **M**ulti-**r**esolution **Hy**bridizable **D**ifferential **E**quations solver ([MrHyDE](./mrhyde/about)), a powerful [C++](https://cplusplus.com/) package for solving partial differential equations at scale.

This project shows how to use SABL and MrHyDE to solve the following types of problems.

:::{important}
This documentation is aimed at advanced undergraduates and graduate students with some MATLAB experience and light familiarity with C++.
The goal is to lower the barrier to entry for large-scale high-performance computing tasks using sophisticated software libraries.
See the [Fundamentals](./fundamentals) page for the necessary mathematical background and references to materials for getting up to speed.
:::

(overview:optimization)=
## Constrained Optimization

The prototypical outer-loop problem is constrained optimization.
WOLF addresses constrained optimization problems of the form

$$
\begin{aligned}
    \min_{u,z} ~& J(u,z)
    \\
    s.t. ~~& c(u,z) = 0
\end{aligned}
$$ (overview:opt_prob)

Here, $u$ indicates the state variable, $z$ is the control variable, $J$ is the objective function, and $c$ specifies the constraints.
If $c(u,z)=0$ admits a unique solution for any given $z$, then there exists a solution operator $S(z)$ such that $c(S(z),z)=0$ for all possible controls $z$.
In this case, {eq}`overview:opt_prob` can be formulated as an equivalent _unconstrained_ optimization problem

::::{margin}
:::{note}
The minimization problem {eq}`overview:rs_opt_prob` is said to be in _reduced space_ because the optimization is over the space of possible controls only, not the joint space of possible states and controls.
:::
::::

$$
    \min_{z} \hat{J}(z)=J(S(z),z).
$$ (overview:rs_opt_prob)

Note that {eq}`overview:opt_prob` and {eq}`overview:rs_opt_prob` may be posed as infinite-dimensional optimizations over function spaces, but the state and control must be discretized to do any computation.
After discretization, we have finite-dimensional state and control vectors

$$
    u\in\mathbb{R}^{n_u},
    \qquad
    z\in\mathbb{R}^{n_z},
$$

where $n_{u}\in\mathbb{R}$ is the state dimension and $n_{z}\in\mathbb{R}$ is the control dimension.

With these dimensions, we have

$$
\begin{aligned}
    &J:\mathbb{R}^{n_{u}}\times\mathbb{R}^{n_{z}}\to\mathbb{R},
    &
    &c:\mathbb{R}^{n_{u}}\times\mathbb{R}^{n_{z}}\to\mathbb{R}^{n_{c}},
    \\
    &S:\mathbb{R}^{n_{z}}\to\mathbb{R}^{n_{u}},
    &
    &\hat{J}:\mathbb{R}^{n_{z}}\to\mathbb{R},
\end{aligned}
$$

where $n_{c}$ is the number of equations needed to express the constraints (typically $n_{c} = n_{u}$).
In problems where the constraints are described by a partial differential equation, the state dimension $n_{u}$ and/or the control dimension $n_{z}$ can be very large.

### Adjoint-based Optimization

The unconstrained problem {eq}`overview:rs_opt_prob` can be solved with off-the-shelf minimizers as long as the gradient $\nabla_z \hat{J}(z)\in\mathbb{R}^{n_{z}}$ and Hessian-vector products $v\mapsto \nabla_{z,z} \hat{J}(z) v$ can be computed efficiently.
To that end, we utilize adjoint-based derivative formulas.
{prf:ref}`alg:adjoint_gradient` shows how to efficiently calculate the gradient, and {prf:ref}`alg:adjoint_hessvec` shows how to compute Hessian-vector products using incremental state and incremental adjoint equations, which avoids explicitly forming the (very large) Hessian matrix $\nabla_{z,z}\hat{J}(z)\in\mathbb{R}^{n_{z}\times n_{z}}$.

In the algorithms and throughout this document, we use $c_u(u,z)\in\mathbb{R}^{n_{c}\times n_{u}}$ and $c_z(u,z)\in\mathbb{R}^{n_{c}\times n_{z}}$ to denote the Jacobians of $c$ with respect to $u$ and $z$, respectively, $\nabla_u J(u, z)\in\mathbb{R}^{n_{u}}$ and $\nabla_z J(u, z)\in\mathbb{R}^{n_{z}}$ to denote the gradients (where $u$ and $z$ are independent of one another) of $J$ with respect to $u$ and $z$, respectively, and a superscript $\mathsf{T}$ to denote a matrix transpose.
We also use $\lambda^{\mathsf{T}} c_{u,u}(u,z) \in \mathbb{R}^{n_{u}\times n_{u}}$ and $\nabla_{u,u} J(u,z)\in\mathbb{R}^{n_{u}\times n_{u}}$ to denote the $(u,u)$ Hessian of the scalar functions $\lambda^{\mathsf{T}} c$ and $J$, respectively, and similar expressions for the $(u,z)$, $(z,u)$, and $(z,z)$ Hessians.

:::{prf:algorithm} Adjoint-based gradient calculation
:label: alg:adjoint_gradient

**Input:** $\overline{z}\in\mathbb{R}^{n_{z}}$

1. Solve the state equation to determine $\overline{u}\in\mathbb{R}^{n_{u}}$ such that $c(\overline{u},\overline{z}) = 0$
2. Solve the adjoint equation $c_u(\overline{u},\overline{z})^{\mathsf{T}} \overline{\lambda} = - \nabla_u J(\overline{u},\overline{z})$ for $\overline{\lambda}\in\mathbb{R}^{n_{c}}$ (a linear solve)
3. Compute $\nabla_z \hat{J}(\overline{z}) = c_z(\overline{u},\overline{z})^{\mathsf{T}} \overline{\lambda} + \nabla_z J(\overline{u},\overline{z})$

**Return**: $\nabla_z \hat{J}(\overline{z}) \in \mathbb{R}^{n_{z}}$
:::

:::{prf:algorithm} Adjoint-based Hessian-vector product calculation
:label: alg:adjoint_hessvec

**Inputs:** $v\in\mathbb{R}^{n_{z}}$, $\overline{u}\in\mathbb{R}^{n_{u}}$, $\overline{z}\in\mathbb{R}^{n_{z}}$, $\overline{\lambda}\in\mathbb{R}^{n_{c}}$

1. Compute $w = c_z(\overline{u},\overline{z}) v \in\mathbb{R}^{n_{c}}$
2. Solve the incremental state equation $c_u(\overline{u},\overline{z}) \overline{\mu} = - w$ for $\overline{\mu} \in \mathbb{R}^{n_{u}}$ (a linear solve)
3. Compute $y_J = \nabla_{u,u} J(\overline{u},\overline{z}) \overline{\mu} + \nabla_{u,z} J(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{u}}$
4. Compute $y_c = \overline{\lambda}^{\mathsf{T}} c_{u,u}(\overline{u},\overline{z}) \overline{\mu} +  \overline{\lambda}^{\mathsf{T}} c_{u,z}(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{u}}$
5. Solve the incremental adjoint equation $c_u(\overline{u},\overline{z})^{\mathsf{T}} \overline{\gamma} = -(y_J + y_c)$ for $\overline{\gamma} \in \mathbb{R}^{n_{c}}$ (a linear solve)
6. Compute $x_J= \nabla_{z,u}J(\overline{u},\overline{z}) \overline{\mu} +  \nabla_{z,z}J(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{z}}$
7. Compute $x_c = c_z(\overline{u},\overline{z})^{\mathsf{T}} \overline{\gamma} + \overline{\lambda}^{\mathsf{T}} c_{z,u}(\overline{u},\overline{z}) \overline{\mu} + \overline{\lambda}^{\mathsf{T}} c_{z,z}(\overline{u},\overline{z})v \in\mathbb{R}^{n_{z}}$
8. Compute $\nabla_{z,z} \hat{J}(\overline{z}) v = x_J + x_c$

**Return**: $\nabla_{z,z} \hat{J}(\overline{z}) v \in \mathbb{R}^{n_{z}}$
:::

These algorithms are implemented in [SABL](./sabl/anatomy) and [MrHyDE](./mrhyde/anatomy).
Examples **TODO** deal with constrained optimization.

### Time-dependent Problems

:::{warning}
The rest of this page is under construction, please check back later.
:::

We are particularly interested in problems where the state $u$ represents the state of an ordinary or partial differential equation and described by the constraints $c(u, z) = 0$.

## Bayesian Inversion

## Optimal Experimental Design (OED)

## Hyper-differential Sensitivity Analysis (HDSA)
