# `Dynamic_Constraint`

:::{admonition} Abstract Methods
:class: abstract

{class}`Dynamic_Constraint` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Dynamic_Constraint.Time_Instance_RHS()`
- {meth}`Dynamic_Constraint.Initial_Condition()`
- {meth}`Dynamic_Constraint.Time_Instance_RHS_yy_Apply()`
- {meth}`Dynamic_Constraint.Time_Instance_RHS_yz_Apply()`
- {meth}`Dynamic_Constraint.Time_Instance_RHS_zy_Apply()`
- {meth}`Dynamic_Constraint.Time_Instance_RHS_zz_Apply()`
- {meth}`Dynamic_Constraint.Initial_Condition_zz_Apply()`

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

    methods (Abstract, Access = public)

        function [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
            error('Time_Instance_RHS() not implemented');
        end

        function [h, h_z] = Initial_Condition(this, z)
            error('Initial_Condition() not implemented');
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)
            error('Time_Instance_RHS_yy_Apply() not implemented');
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
            error('Time_Instance_RHS_yz_Apply() not implemented');
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
            error('Time_Instance_RHS_zy_Apply() not implemented');
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
            error('Time_Instance_RHS_zz_Apply() not implemented');
        end

        function [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)
            error('Initial_Condition_zz_Apply() not implemented');
        end

    end
end
```
