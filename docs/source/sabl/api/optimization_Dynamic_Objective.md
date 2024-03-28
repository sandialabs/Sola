# `Dynamic_Objective`

:::{admonition} Abstract Methods
:class: abstract

{class}`Dynamic_Objective` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Dynamic_Objective.g()`
- {meth}`Dynamic_Objective.R()`
- {meth}`Dynamic_Objective.g_yy_Apply()`
- {meth}`Dynamic_Objective.R_zz_Apply()`

{class}`Dynamic_Objective` implements the following abstract methods from its parent class, i.e., inherited classes **should not** implement these methods.

- {meth}`Objective.J()`
- {meth}`Objective.J_uu_Apply()`
- {meth}`Objective.J_uz_Apply()`
- {meth}`Objective.J_zu_Apply()`
- {meth}`Objective.J_zz_Apply()`

See the [Inheritance Template](optimization.Dynamic_Objective.template) to start a new subclass of {class}`Dynamic_Objective`.
:::

## Class Definition

```{eval-rst}
.. currentmodule:: optimization

.. autoclass:: Dynamic_Objective
   :show-inheritance:
   :members:
```

(optimization.Dynamic_Objective.template)=
## Inheritance Template

```matlab
classdef My_Dynamic_Objective < Dynamic_Objective

    methods (Access = public)

        function [val, grad_y] = g(this, y, t)
            error('g() not implemented');
        end

        function [val, grad_z] = R(this, z)
            error('R() not implemented');
        end

        function [y_out] = g_yy_Apply(this, y_in, y, t)
            error('g_yy_Apply() not implemented');
        end

        function [z_out] = R_zz_Apply(this, z_in, z)
            error('R_zz_Apply() not implemented');
        end

    end
end
```
