# `Dynamic_Constraint`

::::{admonition} Abstract Methods
:class: abstract

{class}`Dynamic_Constraint` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Dynamic_Constraint.f()`
- {meth}`Dynamic_Constraint.h()`

The following methods are not abstract, but they must be implemented to use {meth}`Reduced_Space_Optimization.Optimize()` with the default options ({prf:ref}`alg:adjoint_hessvec`).
We label these _semi-abstract_.
Set `Reduced_Space_Optimization.Gauss_Newton_Hess = true` to use {prf:ref}`alg:adjoint_gaussnewton` instead, which does not rely on these methods.

- {meth}`Dynamic_Constraint.f_yy_Apply()`
- {meth}`Dynamic_Constraint.f_yz_Apply()`
- {meth}`Dynamic_Constraint.f_zy_Apply()`
- {meth}`Dynamic_Constraint.f_zz_Apply()`
- {meth}`Dynamic_Constraint.h_zz_Apply()`

:::{danger}
Because of how MATLAB's [`fminunc()`](https://www.mathworks.com/help/optim/ug/fminunc.html) is designed, the above methods must be implemented in a _vectorized_ fashion, i.e., assuming that `y_in`/`z_in` is a matrix where each column is a test direction.
:::

{class}`Dynamic_Constraint` implements the following abstract methods from its parent class, i.e., inherited classes **should not** implement these methods.

- {meth}`Constraint.State_Solve()`
- {meth}`Constraint.c_u_Transpose_Inverse_Apply()`
- {meth}`Constraint.c_z_Transpose_Apply()`
- {meth}`Constraint.c_u_Inverse_Apply()`
- {meth}`Constraint.c_z_Apply()`
- {meth}`Constraint.c_uu_Apply()`
- {meth}`Constraint.c_uz_Apply()`
- {meth}`Constraint.c_zu_Apply()`
- {meth}`Constraint.c_zz_Apply()`

See the [Inheritance Template](optimization.Dynamic_Constraint.template) to start a new subclass of {class}`Dynamic_Constraint`.
::::

## Class Definition

```{eval-rst}
.. currentmodule:: optimization

.. autoclass:: Dynamic_Constraint
   :show-inheritance:
   :members:
```

(optimization.Dynamic_Constraint.template)=
## Inheritance Template

```matlab
classdef My_Dynamic_Constraint < Dynamic_Constraint

    methods (Access = public)

        function [f, f_y, f_z] = f(this, y, z, t)
            error('f() not implemented');
        end

        function [h, h_z] = h(this, z)
            error('h() not implemented');
        end

        function [y_out] = f_yy_Apply(this, y_in, y, z, t, lambda)
            error('f_yy_Apply() not implemented');
        end

        function [y_out] = f_yz_Apply(this, z_in, y, z, t, lambda)
            error('f_yz_Apply() not implemented');
        end

        function [z_out] = f_zy_Apply(this, y_in, y, z, t, lambda)
            error('f_zy_Apply() not implemented');
        end

        function [z_out] = f_zz_Apply(this, z_in, y, z, t, lambda)
            error('f_zz_Apply() not implemented');
        end

        function [z_out] = h_zz_Apply(this, z_in, z, lambda)
            error('h_zz_Apply() not implemented');
        end

    end
end
```
