## Nonlinear heat/reaction equation

Let $x\in \Omega=(0,1)$, with homogeneous Dirichlet boundary conditions

$$
y(0)=y(1)=0.
$$

Interpret $y(x)$ as a temperature perturbation and $u(x)$ as a distributed heat source/control.

---

## High-fidelity state model

Use a nonlinear diffusion coefficient representing temperature-dependent conductivity:

$$
-\frac{d}{dx}\left(k_\gamma(y)\frac{dy}{dx}\right) + \mu y^3 = u(x),
\qquad x\in(0,1),
$$

where

$$
k_\gamma(y)=1+\gamma y^2.
$$

Thus the high-fidelity model is

$$
-\frac{d}{dx}\left((1+\gamma y^2)y_x\right)+\mu y^3=u(x).
$$

Here:

- $\mu>0$ controls the strength of the nonlinear reaction/radiation term.
- $\gamma\ge 0$ controls the discrepancy between the high- and low-fidelity models.
- When $\gamma=0$, the high- and low-fidelity models coincide.

A reasonable default is

$$
\mu=1,
\qquad
\gamma\in\{0,0.1,0.25,0.5,1.0\}.
$$

---

## Low-fidelity state model

The low-fidelity model assumes constant thermal conductivity, i.e., it neglects the temperature dependence of $k(y)$, but keeps the nonlinear reaction term:

$$
-y_{xx}+\mu y^3=u(x),
\qquad x\in(0,1),
$$

with

$$
y(0)=y(1)=0.
$$

This is a physically plausible approximation: in many heat-transfer settings one linearizes or freezes material properties while retaining some dominant nonlinear source/sink physics.

---

## Model discrepancy

At the continuous operator level, the high-fidelity residual differs from the low-fidelity residual by

$$
\delta_\gamma(y)
=
-\frac{d}{dx}\left((1+\gamma y^2)y_x\right)
+
y_{xx}.
$$

Since

$$
-\frac{d}{dx}\left((1+\gamma y^2)y_x\right)
=
-y_{xx}
-\gamma \frac{d}{dx}\left(y^2 y_x\right),
$$

the discrepancy is

$$
\delta_\gamma(y)
=
-\gamma \frac{d}{dx}\left(y^2 y_x\right).
$$

Expanding,

$$
\delta_\gamma(y)
=
-\gamma\left(2y y_x^2+y^2y_{xx}\right).
$$

So the discrepancy is nonlinear in the state and is controlled directly by $\gamma$.

# Optimal control problem

Define an optimal control problem. Let $u(x)$ be the control, and define a desired target state $y_d(x)$.
Then solve

$$
\min_{u}
\frac{1}{2}
\int_0^1
\left(y(x)-y_d(x)\right)^2\,dx
+
\frac{\alpha}{2}
\int_0^1 u(x)^2\,dx,
$$

subject to either the high-fidelity or low-fidelity state equation.

---

# Remarks
This benchmark gives you:

- nonlinear high-fidelity PDE,
- nonlinear low-fidelity PDE,
- physically motivated simplification,
- discrepancy controlled by a single knob $\gamma$,
- easy 1D implementation,
- compatibility with both inverse problems and optimal control.

A useful extension, if you want a slightly harder discrepancy, is to replace

$$
k_\gamma(y)=1+\gamma y^2
$$

with

$$
k_\gamma(y)=\exp(\gamma y).
$$

Then the low-fidelity model corresponds to the small-$\gamma y$ or constant-conductivity approximation, but the discrepancy becomes asymmetric in positive and negative $y$. For a first test problem, though, $k_\gamma(y)=1+\gamma y^2$ is simpler and more robust.