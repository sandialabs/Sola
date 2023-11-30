# `Constraint`

:::{admonition} Abstract Methods
:class: abstract

{class}`Constraint` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Constraint.State_Solve()`
- {meth}`Constraint.c_u_Transpose_Inverse_Apply()`
- {meth}`Constraint.c_z_Transpose_Apply()`
- {meth}`Constraint.c_u_Inverse_Apply()`
- {meth}`Constraint.c_z_Apply()`

:::{danger}
Because of how MATLAB's [`fminunc()`](https://www.mathworks.com/help/optim/ug/fminunc.html) is designed, the `c_x_XXX()` methods (e.g., `c_z_Apply()`) must be implemented in a _vectorized_ fashion, i.e., assuming that `v` is a matrix where each column is a test direction.
:::

The following equations are not abstract, but they must be implemented to use {meth}`Reduced_Space_Optimization.Optimize()` with the default options ({prf:ref}`alg:adjoint_hessvec`).
We label these _semi-abstract_.
Set `Reduced_Space_Optimization.Gauss_Newton_Hess = true` to use {prf:ref}`alg:adjoint_gaussnewton` instead, which does not rely on these methods.

- {meth}`Constraint.c_uu_Apply()`
- {meth}`Constraint.c_uz_Apply()`
- {meth}`Constraint.c_zu_Apply()`
- {meth}`Constraint.c_zz_Apply()`

Finally, {meth}`Constraint.c()` is not abstract, but it must be implemented in order to use the finite difference check {meth}`Constraint.Finite_Difference_Constraint_Check()`.

See the [Inheritance Template](optimization.Constraint.template) to start a new subclass of {class}`Constraint`.
:::

## Class Definition

```{eval-rst}
.. currentmodule:: optimization

.. autoclass:: Constraint
   :show-inheritance:
   :members:
```

(optimization.Constraint.template)=
## Inheritance Template

```matlab
classdef My_Constraint < Constraint

    methods (Access = public)

        function [u] = State_Solve(this, z)
            error('StateSolve() not implemented');
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            error('c_u_Transpose_Inverse_Apply() not implemented');
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            error('c_z_Transpose_Apply() not implemented');
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            error('c_u_Inverse_Apply() not implemented');
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            error('c_z_Apply() not implemented');
        end

        % The following methods are not required if
        % Reduced_Space_Optimization.Gauss_Newton_Hess=true.
        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            error('c_uu_Apply() not implemented');
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            error('c_uz_Apply() not implemented');
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            error('c_zu_Apply() not implemented');
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            error('c_zz_Apply() not implemented');
        end

        % This method is required for finite difference checks.
        function [c] = c(this, u, v)
            error('c() not implemented');
        end

    end
end
```
