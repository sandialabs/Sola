# Notation

This page lists the mathematical notation used throughout this documentation.

## Dimensions

| Symbol  | Description           |
| :------ | :-------------------- |
| $n_{u}$ | State dimension       |
| $n_{z}$ | Control dimension     |
| $n_{c}$ | Number of constraints |
| $m$     | ODE state dimension (PDE spatial discretization)  |
| $N$     | Number of time steps  |

## Variables

| Symbol   | Description         | Dimension        |
| :------- | :------------------ | :--------------- |
| $u$      | State variable      | $n_{u}$          |
| $z$      | Control variable    | $n_{z}$          |
| $x$      | Space               | $1$, $2$, or $3$ |
| $t$      | Time                | $1$              |
| $y(t)$   | ODE state variable  | $m$              |
| $T$      | Final time          | $1$              |

## Functions

| Symbol        | Description         | Dimensions |
| :------------ | :------------------ | :--------- |
| $J(u,z)$      | Objective function  | $n_{u}\times n_{z}\to 1$ |
| $c(u,z)$      | Constraint function | $n_{u}\times n_{z}\to n_{c}$ |
| $S(z)$        | Solution operator defined by the constraints | $n_{z} \to n_{u}$ |
| $\hat{J}(z) = J(S(z),z)$ | Objective function in reduced space | $n_{z} \to 1$ |
| $f(y(t),z,t)$ | Differential equation function for $y$ | $m \times n_{z} \to m$ |
| $g(y(t),t)$   | Integrand of cost functional | $m\times 1 \to 1$ |
| $h(z)$        | Initial condition of differential equation | $n_{z} \to m$ |
| $R(z)$        | Control-dependent portion of cost functional | $n_{z} \to 1$ |

<!-- ## Discretized Variables -->
