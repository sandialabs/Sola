# Transient Advection-Diffusion in 1D

This directory contains examples for advection-diffusion problems with one spatial dimension.

## Contents

**Optimization Classes**
- [`Adv_Diff.m`](./Adv_Diff.m): Defines an optimization constrained by advection-diffusion dynamics where the control is a source term.
- [`Adv_Diff_Gaussian_Source.m`](./Adv_Diff_Gaussian_Source.m): Defines an optimization problem constrained by advection-diffusion dynamics where the controls trigger Gaussian source terms.

**Drivers**
- [`Driver_Data_Generation.m`](./Driver_Data_Generation.m): Solve the `Adv_Diff_Gaussian_Source` state equation for several random controls.
- [`Driver_MMS.m`](./Driver_MMS.m): Test the `Adv_Diff` state equation solver with a specific control.
- [`Driver_Opt_Test.m`](./Driver_Opt_Test.m): Run finite difference checks for `Adv_Diff_Gaussian_Source`.
- [`Driver_Opt.m`](./Driver_Opt.m): Solve the `Adv_Diff_Gaussian_Source` optimization problem and compare the optimal state to the target.

We detail the advection-diffusion optimization classes after a short overview of the code architecture.

## Refresher: Code Structure

This library defines two main classes for optimization which must be extended in order to define specific problems.
The base class `Constrained_Optimization` is used for steady-state problems of the form

$$
    \min_{u,z} J(u,z)
    \qquad s.t. \qquad
    c(u, z) = 0.
$$

The documentation lists the functions that must be implemented in order to define $J$, $c$, and their derivatives.
For time-dependent problems, we use the subclass `Constrained_ODE_Optimization`, which considers problems with ODE constraints:

$$
    \min_{y,z} \int_{0}^{T}g(y(t),t)dt + R(z)
    \qquad s.t. \qquad
    \frac{dy}{dt} = f(y,z,t), \quad y(0) = h(z).
$$

To use this class, we must implement functions defining $g$, $R$, $f$, $h$, and their gradients and Hessian actions.
The base class translates the integral equation into an objective function $J(u,z)$ by using the trapezoidal rule for integration and defining $u$ to consist of $y(t)$ evaluated over the discretized temporal domain.
The ODE constraints are discretized using the first-order implicit Euler scheme; more sophisticated time steppers would require further modifications to the `Constrained_ODE_Optimization` class (e.g., in the private functions `State_Solve()`, `State_Eq_Time_Step()`, `Nonlinear_Step()`, etc.).

## `Adv_Diff`: Target-based Minimization for an Advection-Diffusion Process

Consider the spatial domain $\Omega = (0, 1)$.
The `Adv_Diff` class implements an optimization constrained by advection-diffusion dynamics over the space-time domain $\Omega\times[0,T]$, written in PDE form as

$$
\begin{align*}
    \min_{\alpha(x,t),z(x,t)} \int_{0}^{T}\int_{\Omega}\frac{1}{2}|\alpha(x,t) - \mu(x,t)|^{2} \:dx\:dt
    \\
    \text{s.t.}\quad
    \alpha_{t} = \alpha_{xx} - \text{Pe}\,\alpha_{x} + z,
    \\ \alpha_{x}(0,t) = \alpha_{x}(1,t) = 0,
\end{align*}
$$

where $\mu(x, t) = t^{2}e^{-50\left(x - \frac{1}{2}\right)^2}$ is the _target function_.
We want the PDE state $\alpha(x,t)$ to match the target function $\mu(x,t)$ as much as possible while still obeying the dynamics.

To put this problem in the ODE form described above, discretize the state in space over $m$ uniformly spaced nodes $x_{1} = 0, \ldots, x_{m} = 1$ with spacing $\delta x = x_{i+1} - x_{i}$.
The spatial integral is handled with a mass matrix $M\in\mathbb{R}^{m}$, i.e.,

$$
    \|\alpha\|_{L^{2}(\Omega)}^{2}
    = \int_\Omega |a(x)|^{2}\:dx
    \approx a(\mathbf{x})^{\mathsf{T}} M a(\mathbf{x}),
$$

where $\mathbf{x} = [~x_{1}~~\cdots~~x_{m}~]^{\mathsf{T}}\in\mathbb{R}^{m}$ collects the spatial nodes.
Writing $y(t) = \alpha(\mathbf{x},t)\in\mathbb{R}^{m}$, the constraint becomes

$$
    M\frac{dy}{dt}
    = -Ay(t) + z(t),
$$

where $A = S + \text{Pe}\,V$ consists of the discrete diffusion ($S$) and advection ($V$) operators and where $\text{Pe} > 0$ is the Peclet number ($\frac{\text{advection}}{\text{diffusion}}$; for now, $\text{Pe} = 1$).
We therefore have the ODE-constrained minimization

$$
\begin{align*}
    \min_{y(t),z(t)} \int_{0}^{T}\frac{1}{2}(y(t) - \mu(\mathbf{x},t))^{\mathsf{T}}M(y(t) - \mu(\mathbf{x},t))\:dt
    \\
    \text{s.t.}\quad
    M\frac{dy}{dt} = -Ay(t) + z(t), \quad y(0) = 0.
\end{align*}
$$

In the generic notation of `Constrained_ODE_Optimization`, we have

$$
\begin{align*}
    g(y(t),t)
    &= \frac{1}{2}(y(t) - \mu(\mathbf{x},t))^{\mathsf{T}}M(y(t) - \mu(\mathbf{x},t)),
    \\
    \nabla_{y}g(y(t),t)
    &= M(y(t) - \mu(\mathbf{x}, t)),
    \\
    R(z) &= 0, \qquad \nabla R(z) = 0,
    \\
    f(y, z, t)
    &= M^{-1}\left(-Ay + Mz\right),
    \\
    \nabla_{y}f(y, z, t)
    &= -M^{-1}A,
    \qquad
    \nabla_{z}f(y, z, t)
    = I,
    \\
    h(z)
    &= 0.
\end{align*}
$$

In addition, all Hessian actions of $g$, $f$, and $h$ are zero except for $\nabla_{y,y} g(y, z, t)v = Mv$.

### Spatial Discretization Details

We use a Galerkin finite element approach.
The matrices $M$, $S$, and $V$ are given by

$$
    M
    = \frac{\delta x}{6}\left(\begin{array}{ccccccc}
    2 & 1 &        &        &        &   &   \\
    1 & 4 &      1 &        &        &   &   \\
      & 1 &      4 &      1 &        &   &   \\
      &   & \ddots & \ddots & \ddots &   &   \\
      &   &        &      1 &      4 & 1 &   \\
      &   &        &        &      1 & 4 & 1 \\
      &   &        &        &        & 1 & 2 \\
    \end{array}\right)
    \in\mathbb{R}^{m\times m},
$$

$$
    S
    = \frac{1}{\delta x}\left(\begin{array}{ccccccc}
         1 & -1 &        &        &        &    &    \\
        -1 &  2 &     -1 &        &        &    &    \\
           & -1 &      2 &     -1 &        &    &    \\
           &    & \ddots & \ddots & \ddots &    &    \\
           &    &        &     -1 &      2 & -1 &    \\
           &    &        &        &     -1 &  2 & -1 \\
           &    &        &        &        & -1 &  1 \\
    \end{array}\right)
    \in \mathbb{R}^{m\times m},
$$

and

$$
    V
    = \frac{1}{2}\left(\begin{array}{ccccccc}
        -1 &   1 &        &        &        &    &    \\
        -1 &   0 &      1 &        &        &    &    \\
           &  -1 &      0 &      1 &        &    &    \\
           &     & \ddots & \ddots & \ddots &    &    \\
           &     &        &     -1 &      0 &  1 &    \\
           &     &        &        &     -1 &  0 &  1 \\
           &     &        &        &        & -1 &  1 \\
    \end{array}\right)
    \in \mathbb{R}^{m\times m}
$$


## `Adv_Diff_Gaussian_Source`: Advection-diffusion with "Bubble" Sources

The `Adv_Diff_Gaussian_Source` class implements a constrained optimization similar to `Adv_Diff`, but with a slightly different control profile.
Consider the one-dimensional radial basis functions indexed by a centering point$x = c$:

$$
    \phi_{c}(x)
    = e^{-200(x - c)^{2}}.
$$

Given a number $n_{c}$ of control nodes, we define the source term $B_{r}\zeta(t)\in\mathbb{R}^{m}$ via the time-dependent vector $\zeta(t)\in\mathbb{R}^{n_{c}}$ and the matrix

$$
    B_{r}
    = \left(\begin{array}{c|c|c}
        & & \\
        \phi_{c_{1}}(\mathbf{x})
        & \cdots &
        \phi_{c_{n_{c}}}(\mathbf{x})
        \\ & &
    \end{array}\right)
    \in \mathbb{R}^{m\times n_{c}},
$$

where the centers $c_{1} < c_{2} < \ldots < c_{n_{c}}$ are equally spaced points in $\Omega$.
The control space thus has dimension $n_{z} = n_{c}N$, i.e., the control profile $z$ dictates the coefficient vector $\xi(t)\in\mathbb{R}^{n_{c}}$ at each of the $N$ time points.

We select the target function $\mu(x,t) = \frac{t^{2}}{5}e^{-10\left(x - \frac{1}{2}\right)^{2}}$.
As before, we have homogeneous Neumann boundary conditions and trivial initial conditions.
The regularization function $R$ in the objective function (which was zero in the previous case) is defined to be

$$
    \frac{\beta}{2}\int_{0}^{T} \|B_{r}\xi(t)\|_{L^{2}(\Omega)}^{2} \:dt
    ~\approx~
    \frac{\beta}{2} \sum_{k=1}^{N}w_{k}\xi(t_{k})^\mathsf{T}B_{r}^\mathsf{T} M B_{r} \xi(t_{k}),
$$

with trapezoidal weights $w_{1} = w_{N} = \frac{\delta t}{2}$, $w_{k} = \delta t$ for $k=2,\ldots, N-1$, with $\delta t = t_{k+1} - t_{k} = \frac{T}{N - 1}$.
