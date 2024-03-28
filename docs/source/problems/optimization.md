# Constrained Optimization

## Problem Statement

The prototypical outer-loop problem is constrained optimization where solving the constraints incurs a significant computational cost.
WOLF addresses constrained optimization problems that can be written as follows.

$$
\begin{aligned}
    \min_{u,z} ~& \J(u,z)
    \\
    s.t. ~~& c(u,z) = 0.
\end{aligned}
$$ (optimization:problem)

::::{margin}
:::{note}
The minimization problem {eq}`optimization:problem_reduced` is said to be in _reduced space_ because the optimization is over the space of possible controls only, not the joint space of possible states and controls.
:::
::::

Here, $u$ indicates the state variable, $z$ is the control variable, $\J$ is the objective function, and $c$ specifies the constraints.
If $c(u,z)=0$ admits a unique solution for any given $z$, then there exists a solution operator $S(z)$ such that $c(S(z),z)=0$ for all possible controls $z$.
In this case, {eq}`optimization:problem` can be formulated as an equivalent _unconstrained_ optimization problem

$$
    \min_{z} \hat{\J}(z) = \J(S(z),z).
$$ (optimization:problem_reduced)

The state and control may be finite-dimensional vectors, or they may be functions of space and/or time.
In the latter case, {eq}`optimization:problem` and {eq}`optimization:problem_reduced` may be posed as infinite-dimensional optimizations over function spaces.
However, the state and control must be discretized (represented with finitely many numbers) to do any computation.
After discretization, we have finite-dimensional state and control vectors

$$
    \u\in\R^{n_u},
    \qquad
    \z\in\R^{n_z},
$$

where $n_u\in\mathbb{N}$ is the state dimension and $n_z\in\mathbb{N}$ is the control dimension.
Discretization also introduces finite-dimensional versions of the functions $\J$, $c$, $S$, and $\hat{\J}$.

:::{table}
:align: center

| Original function | Discretized function | Dimensions                             |
| ----------------: | :------------------- | :------------------------------------- |
| $\J(u, z)$        | $J(\u,\z)$           | $J:\R^{n_u}\times\R^{n_z}\to\R$        |
| $c(u, z)$         | $\c(\u,\z)$          | $\c:\R^{n_u}\times\R^{n_z}\to\R^{n_u}$ |
| $S(u, z)$         | $\S(\z)$             | $\S:\R^{n_z}\to\R^{n_u}$               |
| $\hat{\J}(u, z) = \J(S(z), z)$ | $\hat{J}(\z) = J(\S(\z), \z)$          | $\hat{J}:\R^{n_z}\to\R$ |
:::

The discretized version of {eq}`optimization:problem`--{eq}`optimization:problem_reduced` is given by

$$
\begin{aligned}
    \min_{\u,\z} ~& J(\u,\z)
    \qquad\Longleftrightarrow\qquad
    \min_{\z}\hat{J}(\z) = J(\S(\z),\z),
    \\
    s.t. ~~& \c(\u,\z) = \0,
\end{aligned}
$$ (optimization:problem_discrete)

where $\S$ satisfies $\c(\S(\z), \z) = \0$ for all admissible controls $\z$.

:::{important}
In many settings, such as problems where the constraints are described by a partial differential equation, the state dimension $n_u$ and/or the control dimension $n_z$ can be very large.
This makes solving complex problems of this type at scale computationally challenging.
:::

## Adjoint-based Optimization

The minimization problem {eq}`optimization:problem_discrete` can be solved with off-the-shelf minimizers as long as the gradient $\grad{z}\hat{J}(\z)\in\R^{n_z}$ and, for any vector $\v\in\R^{n_z}$, Hessian-vector products $\v\mapsto \grad{z,z}\hat{J}(\z) \v \in \R^{n_z}$ can be computed efficiently.
{prf:ref}`alg:adjoint_gradient` shows how to efficiently calculate the gradient by utilizing adjoint-based derivative formulas.
Similarly, {prf:ref}`alg:adjoint_hessvec` shows how to compute Hessian-vector products using incremental state and incremental adjoint equations, which avoids explicitly forming the (very large) Hessian matrix $\grad{z,z}\hat{J}(\z)\in\R^{n_z \times n_z}$.
These algorithms stem from analysis of the original (possible infinite-dimensional) problem {eq}`optimization:problem`, then applying the discretization.

:::{admonition} Notation
:class: note

In the algorithms and throughout this project, we use the following notation.

- A superscript $\mathsf{T}$ is the matrix transpose.
- $\grad{u}J(\u,\z)\in\R^{n_u}$ and $\grad{z}J(\u, \z)\in\R^{n_z}$ denote the gradients of $J$ with respect to $\u$ and $\z$, respectively.
- $\c_u(\u,\z)\in\R^{n_u \times n_u}$ and $\c_z(\u,\z)\in\R^{n_u \times n_z}$ denote the Jacobians of $\c$ with respect to $\u$ and $\z$, respectively.
- $\bflambda\trp\c_{u,u}(\u,\z) \in \R^{n_u \times n_u}$ and $\grad{u,u}J(\u,\z)\in\R^{n_u \times n_u}$ denote the $(u,u)$ Hessian of the scalar functions $\bflambda\trp\c(\u,\z)$ and $J(\u,\z)$, respectively. Similar expressions are used for the $(u,z)$, $(z,u)$, and $(z,z)$ Hessians.

See [Notation](../appendix/notation) for a comprehensive list of notation and [Fundamentals](../appendix/fundamentals) for the definitions of these derivatives.
:::

:::{prf:algorithm} Adjoint-based gradient calculation
:label: alg:adjoint_gradient

**Input:** Control $\z\in\R^{n_z}$

1. Solve the state equation to determine $\u\in\R^{n_u}$ such that $\c(\u,\z) = \0$ (i.e., $\u = \S(\z)$)
2. Solve the (linear) adjoint equation $\c_u(\u,\z)\trp \bflambda = - \grad{u}J(\u,\z)$ for $\bflambda\in\R^{n_u}$
3. Set $\grad{z}\hat{J}(\z) = \c_z(\u,\z)\trp \bflambda + \grad{z}J(\u,\z)$

**Return**: $\grad{z}\hat{J}(\z) \in \R^{n_z}$
:::

:::{prf:algorithm} Adjoint-based Hessian-vector product calculation
:label: alg:adjoint_hessvec

**Inputs:** State $\u\in\R^{n_u}$, control $\z\in\R^{n_z}$, adjoint $\bflambda\in\R^{n_u}$, search direction $\v\in\R^{n_z}$

1. Compute $\w = \c_z(\u,\z) \v \in\R^{n_u}$
2. Solve the (linear) incremental state equation $\c_u(\u,\z) \bfmu = -\w$ for $\bfmu \in \R^{n_u}$
3. Compute $\x_J= \grad{z,u}J(\u,\z) \bfmu +  \grad{z,z}J(\u,\z) \v \in \R^{n_z}$
4. Compute $\y_J = \grad{u,u}J(\u,\z) \bfmu + \grad{u,z}J(\u,\z) \v \in \R^{n_u}$
5. Compute $\y_c = \bflambda\trp \c_{u,u}(\u,\z) \bfmu +  \bflambda\trp \c_{u,z}(\u,\z) \v \in \R^{n_u}$
6. Solve the (linear) incremental adjoint equation $\c_u(\u,\z)\trp \bfgamma = -(\y_J + \y_c)$ for $\bfgamma \in \R^{n_u}$
7. Compute $\x_c = \c_z(\u,\z)\trp \bfgamma + \bflambda\trp \c_{z,u}(\u,\z) \bfmu + \bflambda\trp \c_{z,z}(\u,\z)\v \in\R^{n_z}$
8. Set $\grad{z,z}\hat{J}(\z) \v = \x_J + \x_c$

**Return**: $\grad{z,z}\hat{J}(\z) \v \in \R^{n_z}$
:::

The Gauss--Newton iteration is an alternative to {prf:ref}`alg:adjoint_hessvec` that requires only first derivatives of $\c$.
The steps marked with an asterisk are where this approach differs from {prf:ref}`alg:adjoint_hessvec`.

:::{prf:algorithm} Gauss-Newton Hessian-vector product calculation
:label: alg:adjoint_gaussnewton

**Inputs:** State $\u\in\R^{n_u}$, control $\z\in\R^{n_z}$, adjoint $\bflambda\in\R^{n_u}$, search direction $\v\in\R^{n_z}$

1. Compute $\w = \c_z(\u,\z) \v \in\R^{n_u}$
2. Solve the (linear) incremental state equation $\c_u(\u,\z) \bfmu = -\w$ for $\bfmu \in \R^{n_u}$
3. Compute $\x_J= \grad{z,u}J(\u,\z) \bfmu +  \grad{z,z}J(\u,\z) \v \in \R^{n_z}$
4. Compute $\y_J = \grad{u,u}J(\u,\z) \bfmu + \grad{u,z}J(\u,\z) \v \in \R^{n_u}$
5. ${}^\ast$Solve the (linear) incremental adjoint equation $\c_u(\u,\z)\trp \bfgamma = -\y_J$ for $\bfgamma \in \R^{n_u}$
6. ${}^\ast$Compute $\x_c = \c_z(\u,\z)\trp \bfgamma$
7. Set $\grad{z,z}\hat{J}(\z) \v = \x_J + \x_c$

**Return**: $\grad{z,z}\hat{J}(\z) \v \in \R^{n_z}$
:::

These algorithms are implemented in [SABL](../sabl/guides/optimization) and [MrHyDE](../mrhyde/guides/optimization).
Tutorials 1--**TODO** deal with constrained optimization.

(optimization:differential)=
## Differential Equation Constraints

Of particular interest are problems where the state $u$ is the solution of an ordinary or partial differential equation.
In this case, the constraints $\c(\u, \z) = \0$ describe the differential equation and the solution mapping $\S:\z\mapsto \u$ means solving the differential equation with the specified control.
For steady-state problems, the discrete state $\u$ represents a spatial discretization of the state $u$;
for time-dependent problems, $\u$ includes the spatial discretization of the $u$ _at each time_ in the temporal discretization.
In the latter setting, we use $\y(t)$ to denote the state at fixed time $t$.
For example, suppose the state $u(x,t)$ is to be computed at spatial points $x_{1},\ldots,x_{n_y}$ and temporal points $t_{1},\ldots,t_{n_t}$.
Then the discrete state $\u$ is given by

$$
\begin{align*}
    \u =
    \left(\begin{array}{c}
        \y_1 \\ \y_2 \\ \vdots \\ \y_{n_t}
    \end{array}\right)
    = \left(\begin{array}{c}
        y_{1,1} \\ y_{1,2} \\ \vdots \\ y_{1,n_t} \\
        y_{2,1} \\ y_{2,2} \\ \vdots \\ y_{2,n_t} \\
        \vdots \\
        y_{n_y,1} \\ y_{n_y,2} \\ \vdots \\ y_{n_y,n_t} \\
    \end{array}\right)
    \approx
    \left(\begin{array}{c}
        u(x_{1},t_{1}) \\ u(x_{2},t_{1}) \\ \vdots \\ u(x_{n_y},t_{1}) \\
        u(x_{1},t_{2}) \\ u(x_{2},t_{2}) \\ \vdots \\ u(x_{n_y},t_{2}) \\
        \vdots \\
        u(x_{1},t_{n_t}) \\ u(x_{2},t_{n_t}) \\ \vdots \\ u(x_{n_y},t_{n_t})
    \end{array}\right)
    \in \R^{n_u},
\end{align*}
$$

where $\y_j = (y_{j,1},y_{j,2},\ldots,y_{j,n_t})\trp\in\R^{n_t}$ with $y_{j,k} \approx u(x_k, t_j)$.
The discrete state dimension is $n_u = n_y n_t$.

The following ODE-constrained optimization problem is a common occurrence of {eq}`optimization:problem_discrete` for scientific problems.

$$
\begin{align*}
    \min_{\y(t),\z} ~& \int_{0}^{T} g(\y(t),t) dt + R(\z)
    \\
    s.t. ~~& \frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t), \z, t), ~~ \y(0) = \h(\z),
\end{align*}
$$

where $T>0$ and

$$
\begin{align*}
    \y &: [0,T] \to \R^{n_y},
    &
    g &: \R^{n_y} \times [0, T] \to \R,
    &
    R &: \R^{n_z} \to \R,
    \\
    \f &: \R^{n_y} \times \R^{n_z} \times [0,T]  \to \R^m,
    &
    \h &: \R^{n_z} \to \R^{n_y}.
\end{align*}
$$

In this case, the derivatives of $J$ and $\c$ needed to implement {prf:ref}`alg:adjoint_gradient` and {prf:ref}`alg:adjoint_hessvec` can be expressed in terms of the derivatives of $g$, $R$, $\f$, and $\h$.
