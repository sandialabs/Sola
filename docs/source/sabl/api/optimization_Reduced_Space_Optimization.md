# `optimization.Reduced_Space_Optimization`

## Class Definition

```{eval-rst}
.. currentmodule:: optimization

.. autoclass:: Reduced_Space_Optimization
   :show-inheritance:
   :members:
```

## Optimization Customization

{meth}`Reduced_Space_Optimization.Optimize()` can be customized by setting the object properties, most of which correspond to an option for MATLAB's [`fminunc()`](https://www.mathworks.com/help/optim/ug/fminunc.html) or [`fmincon()`](https://www.mathworks.com/help/optim/ug/fmincon.html).

:::{table}
:align: center
:name: tab-optimization-customization

| Property            | Description                                       | `fminunc()` argument    | Default value |
| :------------------ | :------------------------------------------------ | :---------------------- | :------------ |
| `Gauss_Newton_Hess` | Use the Gauss-Newton approximation of the Hessian ({prf:ref}`alg:adjoint_gaussnewton`) | | `false` |
| `use_trust_region`  | Use trust region for the optimization             |                         | `true`        |
| `z_lb`              | Lower bounds for control                          |                         | `[]`          |
| `z_ub`              | Upper bounds for control                          |                         | `[]`          |
| `opt_tol`           | Optimality tolerance                              | `'OptimalityTolerance'` | `10^-8`       |
| `fun_tol`           | Function tolerance                                | `'FunctionTolerance'`   | `10^-6`       |
| `iteration_limit`   | Maximum number of iterations                      | `'MaxIterations'`       | `1000`        |
| `step_tol`          | Step tolerance                                    | `'StepTolerance'`       | `10^-6`       |
| `max_cg_iter`       | Maximum number of conjugate gradient iterations   | `'MaxPCGIter'`          | `50`          |
| `cg_tol`            | Conjugate gradient tolerance                      | `'TolPCG'`              | `10^-4`       |
| `verbose`           | Print optimization info at each iteration         | `'Display'`             | `true` (`'iter-detailed'`) |
:::

For example, suppose we have an objective object `obj` and a Constraint object `con`.
The following code customizes an optimization to use the Gauss-Newton Hessian approximation without a trust region algorithm.

```matlab
% Define the optimization problem.
opt = Reduced_Space_Optimization(obj, con);

% Set custom optimization options.
opt.iteration_limit = 200;
opt.Gauss_Newton_Hess = true;
opt.use_trust_region = false;
opt.verbose = false;

% Solve the optimization problem.
[u_solution, z_solution] = opt.Optimize(z0);
```
