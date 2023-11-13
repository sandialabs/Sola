# Notation

This page lists the mathematical notation used throughout this project.

## Dimensions

| Symbol  | Description                                                                          |
| :------ | :----------------------------------------------------------------------------------- |
| $n_{u}$ | Dimension of the state after discretization                                          |
| $n_{z}$ | Dimension of the control after discretization                                        |
| $n_{y}$ | ODE state dimension, i.e., number of degrees of freedom in the state at a fixed time |
| $n_{t}$ | Number of time steps in a temporal discretization for time-dependent problems        |

## Variables

| Symbol   | Description         | Dimension        |
| :------- | :------------------ | :--------------- |
| $u$      | State               |                  |
| $z$      | Control             |                  |
| | | |
| $\u$     | Discretized state   | $n_{u}$          |
| $\z$     | Discretized control | $n_{z}$          |
| $x$      | Space               | $1$, $2$, or $3$ |
| $t$      | Time                | $1$              |
| $T$      | Final time          | $1$              |
| $\y$     | ODE state           | $n_y$            |
| $\bflambda$ | Adjoint of state (discrete) | $n_y$ |
| $\bfmu$ | Incremental adjoint (discrete)  | TODO  |

## Functions

| Symbol           | Description                           | Dimensions                        |
| :--------------- | :------------------------------------ | :-------------------------------- |
| $\mathcal{J}(u,z)$ | Objective                           |                                   |
| $\hat{J}(z) = J(S(z),z)$ | Reduced-space objective       |                                   |
| $c(u,z)$         | Constraint                            |                                   |
| $S(z)$           | Solution operator for the constraints |                                   |
|                  |                                       |                                   |
| $J(\u,\z)$       | Discretized objective                 | $n_u \times n_z \to 1$            |
| $c(u,z)$         | Discretized constraint                | $n_u \times n_z \to n_u$          |
| $S(z)$           | Discretized solution operator         | $n_z \to n_u$                     |
| $g(\y(t),t)$     | Integrand of objective for time-dependent problems | $n_y \times 1 \to 1$ |
| $R(\z)$          | Control-dependent portion of objective for time-dependent problems | $n_{z} \to 1$ |
| $\f(\y(t),\z,t)$ | Differential equation constraint      | $n_y \times n_z \times 1 \to n_y$ |
| $\h(\z)$         | Initial condition of ODE constraint   | $n_{z} \to n_y$                   |

## Derivatives

For a scalar-valued function $J(\u,\z) \to \R$ where $\u = (u_1,\ldots,u_{n_u})\trp\in\R^{n_u}$ and $\z = (z_1,\ldots,z_{n_u})\trp\in\R^{n_z}$, we denote the derivatives of $J$ as follows.

| Symbol             | Description                             | Dimensions       | $i$-th entry      |
| :----------------- | :-------------------------------------- | :--------------- | :---------------- |
| $\grad{u}J(\u,\z)$ | Gradient of $J$ w.r.t. the state $\u$   | $n_u \times 1$   | $\frac{\partial}{\partial u_{i}}J(\u,\z)$ |
| $\grad{z}J(\u,\z)$ | Gradient of $J$ w.r.t. the control $\z$ | $n_z \times 1$   | $\frac{\partial}{\partial z_{i}}J(\u,\z)$ |
|                    |                                         |                  | $ij$**-th entry** |
| $\grad{u,u}J(\u,\z)$ | Hessian of $J$ w.r.t. $\u$ then $\u$  | $n_u \times n_u$ | $\frac{\partial^2}{\partial u_i \partial u_j}J(\u,\z)$ |
| $\grad{u,z}J(\u,\z)$ | Hessian of $J$ w.r.t. $\z$ then $\u$  | $n_u \times n_z$ | $\frac{\partial^2}{\partial u_i \partial z_j}J(\u,\z)$ |
| $\grad{z,u}J(\u,\z)$ | Hessian of $J$ w.r.t. $\u$ then $\z$  | $n_z \times n_u$ | $\frac{\partial^2}{\partial z_i \partial u_j}J(\u,\z)$ |
| $\grad{z,z}J(\u,\z)$ | Hessian of $J$ w.r.t. $\z$ then $\z$  | $n_z \times n_z$ | $\frac{\partial^2}{\partial z_i \partial z_j}J(\u,\z)$ |

For a vector-valued function $\c(\u, \z) \to \R^{n}$ with entries $\c(\u,\z) = (c_1(\u, \z),\ldots,c_n(\u,\z))\trp$, we denote the derivatives of $\c$ as follows.

| Symbol            | Description                              | Dimensions                  | $ij$-th entry      |
| :---------------- | :--------------------------------------- | :-------------------------- | :----------------- |
| $\c_u(\u,\z)$     | Jacobian of $\c$ w.r.t. the state $\u$   | $n_u \times n_u$            | $\frac{\partial}{\partial u_{i}}J(\u,\z)$ |
| $\c_z(\u,\z)$     | Jacobian of $\c$ w.r.t. the control $\z$ | $n_u \times n_z$            | $\frac{\partial}{\partial z_{i}}J(\u,\z)$ |
|                   |                                          |                             | $ijk$**-th entry** |
| $\c_{u,u}(\u,\z)$ | Hessian of $\c$ w.r.t. $\u$ then $\u$    | $n_u \times n_u \times n_u$ | $\frac{\partial^2}{\partial u_j \partial u_k}c_i(\u,\z)$ |
| $\c_{u,z}(\u,\z)$ | Hessian of $\c$ w.r.t. $\z$ then $\u$    | $n_u \times n_u \times n_z$ | $\frac{\partial^2}{\partial u_j \partial z_k}c_i(\u,\z)$ |
| $\c_{z,u}(\u,\z)$ | Hessian of $\c$ w.r.t. $\u$ then $\z$    | $n_u \times n_z \times n_u$ | $\frac{\partial^2}{\partial z_j \partial u_k}c_i(\u,\z)$ |
| $\c_{z,z}(\u,\z)$ | Hessian of $\c$ w.r.t. $\z$ then $\z$    | $n_u \times n_z \times n_z$ | $\frac{\partial^2}{\partial z_j \partial z_k}c_i(\u,\z)$ |
