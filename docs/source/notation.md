# Notation

This page lists the mathematical notation used throughout this documentation.

## Dimensions

| Symbol  | Description           |
| :------ | :-------------------- |
| $n_{u}$ | State dimension       |
| $n_{z}$ | Control dimension     |
| $n_{c}$ | Number of constraints |
| $m$     | ODE State dimension   |

## Variables

| Symbol   | Description         | Dimension |
| :------- | :------------------ | :-------- |
| $u$      | State variable      | $n_{u}$   |
| $z$      | Control variable    | $n_{z}$   |
| $t$      | Time                | $1$       |
| $y(t)$   | ODE state variable  | $m$       |

## Functions

| Symbol       | Description         | Dimensions |
| :----------- | :------------------ | :--------- |
| $J(u,z)$     | Objective function  | $n_{u}\times n_{z}\to 1$ |
| $c(u,z)$     | Constraint function | $n_{u}\times n_{z}\to n_{c}$ |
| $S(z)$       | Solution operator defined by the constraints | $n_{z} \to n_{u}$ |
| $\hat{J}(z) = J(S(z),z)$ | Objective function in reduced space | $n_{z} \to 1$ |
| $g(t)$   | TODO                | TODO      |
