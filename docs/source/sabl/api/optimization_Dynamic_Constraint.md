# `Dynamic_Constraint`

:::{admonition} Abstract Methods
:class: abstract

{class}`Dynamic_Constraint` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Dynamic_Constraint.f()`
- {meth}`Dynamic_Constraint.h()`
- {meth}`Dynamic_Constraint.f_yy_Apply()`
- {meth}`Dynamic_Constraint.f_yz_Apply()`
- {meth}`Dynamic_Constraint.f_zy_Apply()`
- {meth}`Dynamic_Constraint.f_zz_Apply()`
- {meth}`Dynamic_Constraint.h_zz_Apply()`

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
:::

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
