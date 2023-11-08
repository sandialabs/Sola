# Constrained Optimization

The prototypical outer-loop problem is constrained optimization where solving the constraints incurs a significant computational cost.
WOLF addresses constrained optimization problems that can be written as follows.

$$
\begin{aligned}
    \min_{u,z} ~& J(u,z)
    \\
    s.t. ~~& c(u,z) = 0.
\end{aligned}
$$ (optimization:problem)

Here, $u$ indicates the state variable, $z$ is the control variable, $J$ is the objective function, and $c$ specifies the constraints.
If $c(u,z)=0$ admits a unique solution for any given $z$, then there exists a solution operator $S(z)$ such that $c(S(z),z)=0$ for all possible controls $z$.
In this case, {eq}`optimization:problem` can be formulated as an equivalent _unconstrained_ optimization problem

::::{margin}
:::{note}
The minimization problem {eq}`optimization:problem_reduced` is said to be in _reduced space_ because the optimization is over the space of possible controls only, not the joint space of possible states and controls.
:::
::::

$$
    \min_{z} \hat{J}(z)=J(S(z),z).
$$ (optimization:problem_reduced)

The state and control may be finite-dimensional vectors, or they may be functions of space and/or time.
In the latter case, {eq}`optimization:problem` and {eq}`optimization:problem_reduced` may be posed as infinite-dimensional optimizations over function spaces.
However, the state and control must be discretized (represented with finitely many numbers) to do any computation.
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

## Adjoint-based Optimization

The unconstrained problem {eq}`optimization:problem_reduced` can be solved with off-the-shelf minimizers as long as the gradient $\nabla_z \hat{J}(z)\in\mathbb{R}^{n_{z}}$ and Hessian-vector products $v\mapsto \nabla_{z,z} \hat{J}(z) v$ can be computed efficiently.
{prf:ref}`alg:adjoint_gradient` shows how to efficiently calculate the gradient by utilizing adjoint-based derivative formulas.
Similarly, {prf:ref}`alg:adjoint_hessvec` shows how to compute Hessian-vector products using incremental state and incremental adjoint equations, which avoids explicitly forming the (very large) Hessian matrix $\nabla_{z,z}\hat{J}(z)\in\mathbb{R}^{n_{z}\times n_{z}}$.

:::{admonition} Notation
In the algorithms and throughout this project, we use the following notation.

- a superscript $\mathsf{T}$ is the matrix transpose.
- $c_u(u,z)\in\mathbb{R}^{n_{c}\times n_{u}}$ and $c_z(u,z)\in\mathbb{R}^{n_{c}\times n_{z}}$ denote the Jacobians of $c$ with respect to $u$ and $z$, respectively.
- $\nabla_u J(u, z)\in\mathbb{R}^{n_{u}}$ and $\nabla_z J(u, z)\in\mathbb{R}^{n_{z}}$ denote the gradients (where $u$ and $z$ are independent of one another) of $J$ with respect to $u$ and $z$, respectively.
- $\lambda^{\mathsf{T}} c_{u,u}(u,z) \in \mathbb{R}^{n_{u}\times n_{u}}$ and $\nabla_{u,u} J(u,z)\in\mathbb{R}^{n_{u}\times n_{u}}$ denote the $(u,u)$ Hessian of the scalar functions $\lambda^{\mathsf{T}} c$ and $J$, respectively. Similar expressions are for the $(u,z)$, $(z,u)$, and $(z,z)$ Hessians.
:::

:::{prf:algorithm} Adjoint-based gradient calculation
:label: alg:adjoint_gradient

**Input:** $\overline{z}\in\mathbb{R}^{n_{z}}$

1. Solve the state equation to determine $\overline{u}\in\mathbb{R}^{n_{u}}$ such that $c(\overline{u},\overline{z}) = 0$
2. Solve the adjoint equation $c_u(\overline{u},\overline{z})^{\mathsf{T}} \overline{\lambda} = - \nabla_u J(\overline{u},\overline{z})$ for $\overline{\lambda}\in\mathbb{R}^{n_{c}}$ (a linear solve)
3. Set $\nabla_z \hat{J}(\overline{z}) = c_z(\overline{u},\overline{z})^{\mathsf{T}} \overline{\lambda} + \nabla_z J(\overline{u},\overline{z})$

**Return**: $\nabla_z \hat{J}(\overline{z}) \in \mathbb{R}^{n_{z}}$
:::

:::{prf:algorithm} Adjoint-based Hessian-vector product calculation
:label: alg:adjoint_hessvec

**Inputs:** $v\in\mathbb{R}^{n_{z}}$, $\overline{u}\in\mathbb{R}^{n_{u}}$, $\overline{z}\in\mathbb{R}^{n_{z}}$, $\overline{\lambda}\in\mathbb{R}^{n_{c}}$

1. Compute $w = c_z(\overline{u},\overline{z}) v \in\mathbb{R}^{n_{c}}$
2. Solve the incremental state equation $c_u(\overline{u},\overline{z}) \overline{\mu} = - w$ for $\overline{\mu} \in \mathbb{R}^{n_{u}}$ (a linear solve)
3. Compute $y_J = \nabla_{u,u} J(\overline{u},\overline{z}) \overline{\mu} + \nabla_{u,z} J(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{u}}$
4. Compute $x_J= \nabla_{z,u}J(\overline{u},\overline{z}) \overline{\mu} +  \nabla_{z,z}J(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{z}}$
5. Compute $y_c = \overline{\lambda}^{\mathsf{T}} c_{u,u}(\overline{u},\overline{z}) \overline{\mu} +  \overline{\lambda}^{\mathsf{T}} c_{u,z}(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{u}}$
6. Solve the incremental adjoint equation $c_u(\overline{u},\overline{z})^{\mathsf{T}} \overline{\gamma} = -(y_J + y_c)$ for $\overline{\gamma} \in \mathbb{R}^{n_{c}}$ (a linear solve)
7. Compute $x_c = c_z(\overline{u},\overline{z})^{\mathsf{T}} \overline{\gamma} + \overline{\lambda}^{\mathsf{T}} c_{z,u}(\overline{u},\overline{z}) \overline{\mu} + \overline{\lambda}^{\mathsf{T}} c_{z,z}(\overline{u},\overline{z})v \in\mathbb{R}^{n_{z}}$
8. Set $\nabla_{z,z} \hat{J}(\overline{z}) v = x_J + x_c$

**Return**: $\nabla_{z,z} \hat{J}(\overline{z}) v \in \mathbb{R}^{n_{z}}$
:::

The Gauss--Newton iteration is an alternative to {prf:ref}`alg:adjoint_hessvec` that requires only first derivatives of $c$.
The steps marked with an asterisk are where this approach differs from {prf:ref}`alg:adjoint_hessvec`.

:::{prf:algorithm} Gauss-Newton Hessian-vector product calculation
:label: alg:adjoint_gaussnewton

**Inputs:** $v\in\mathbb{R}^{n_{z}}$, $\overline{u}\in\mathbb{R}^{n_{u}}$, $\overline{z}\in\mathbb{R}^{n_{z}}$, $\overline{\lambda}\in\mathbb{R}^{n_{c}}$

1. Compute $w = c_z(\overline{u},\overline{z}) v \in\mathbb{R}^{n_{c}}$
2. Solve the incremental state equation $c_u(\overline{u},\overline{z}) \overline{\mu} = - w$ for $\overline{\mu} \in \mathbb{R}^{n_{u}}$ (a linear solve)
3. Compute $y_J = \nabla_{u,u} J(\overline{u},\overline{z}) \overline{\mu} + \nabla_{u,z} J(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{u}}$
4. Compute $x_J= \nabla_{z,u}J(\overline{u},\overline{z}) \overline{\mu} +  \nabla_{z,z}J(\overline{u},\overline{z}) v \in \mathbb{R}^{n_{z}}$
5. ${}^\ast$Solve the incremental adjoint equation $c_u(\overline{u},\overline{z})^{\mathsf{T}} \overline{\gamma} = -y_J$ for $\overline{\gamma} \in \mathbb{R}^{n_{c}}$ (a linear solve)
6. ${}^\ast$Compute $x_c = c_z(\overline{u},\overline{z})^{\mathsf{T}} \overline{\gamma}$
7. Set $\nabla_{z,z} \hat{J}(\overline{z}) v = x_J + x_c$

**Return**: $\nabla_{z,z} \hat{J}(\overline{z}) v \in \mathbb{R}^{n_{z}}$
:::

These algorithms are implemented in [SABL](../sabl/anatomy.md) and [MrHyDE](../mrhyde/anatomy).
Examples 1--**TODO** deal with constrained optimization.

## Differential Equation Constraints

We are particularly interested in problems where the state $u$ is the solution of an ordinary or partial differential equation.
In this case, the constraints $c(u, z) = 0$ describe the differential equation and the solution mapping $S:z\mapsto u$ means solving the differential equation with the specified control profile.
In steady-state problems, the computational state consists of the spatial discretization of the analytical state $u$;
in time-dependent problems, the computational state includes the spatial discretization of the analytical state _at each time_ in the temporal discretization.
For example, suppose the state $u(x,t)$ is to be computed at spatial points $x_{1},\ldots,x_{m}$ and temporal points $t_{1},\ldots,t_{N}$.
Then the full computational state approximates the vector

$$
\begin{align*}
    \left(\begin{array}{c}
        u(x_{1},t_{1}) \\ u(x_{2},t_{1}) \\ \vdots \\ u(x_{m},t_{1}) \\
        u(x_{1},t_{2}) \\ u(x_{2},t_{2}) \\ \vdots \\ u(x_{m},t_{2}) \\
        \vdots \\
        u(x_{1},t_{N}) \\ u(x_{2},t_{N}) \\ \vdots \\ u(x_{m},t_{N})
    \end{array}\right)
    \in \mathbb{R}^{n_{u}}
\end{align*}
$$

and the computational state dimension is $n_{u} = mN$.
Depending on the problem, the control may be similarly represented as vector with values for space and time.

The following ODE-constrained optimization problem is a common occurrence of {eq}`optimization:problem` for scientific problems.

$$
\begin{align*}
    \min_{y,z} ~& \int_{0}^{T} g(y(t),t) dt + R(z)
    \\
    s.t. ~~& \frac{dy}{dt} = f(y(t), z, t), ~~ y(0) = h(z)
\end{align*}
$$

where $T>0$ is the final time, $y$ is the differential equation state, and

$$
\begin{align*}
    y &: [0,T] \to \mathbb{R}^m,
    &
    f &: \mathbb{R}^m \times \mathbb{R}^{n_{z}} \times [0,T]  \to \mathbb{R}^m,
    &
    G &: \mathbb{R}^{m} \times [0, T] \to \mathbb{R}
    &
    h &: \mathbb{R}^{n_{z}} \to \mathbb{R}^m.
\end{align*}
$$

In this setting, the full optimization state $u$ consists of the ODE state $y(t)$ evaluated at all monitored times $t$.

:::{admonition} TODO
Need to clarify the difference between the infinite-dimensional (function analytic) setting and the finite-dimensional (computational) setting.
Perhaps use $u$ for the _state_ $\mathbf{u}$ for the _computational state_? If so, we can write $u = y$ but $\mathbf{u} = (\mathbf{y}_{1}^{\mathsf{T}}~\cdots~\mathbf{y}_{N}^{\mathsf{T}})^{\mathsf{T}}$ where $\mathbf{y}_{j}$ approximates $(u(x_{1}, t_{j})~\cdots~u(x_{m}, t_{j}))^{\mathsf{T}}$.
:::

:::{warning}
The rest of this page is under construction, please check back later.
:::
