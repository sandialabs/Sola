# `optimization.Dynamic_Objective`

:::{admonition} Abstract Methods
:class: abstract

{class}`Dynamic_Objective` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Dynamic_Objective.Time_Instance_Objective()`
- {meth}`Dynamic_Objective.Regularization_Objective()`
- {meth}`Dynamic_Objective.Time_Instance_Objective_yy_Apply()`
- {meth}`Dynamic_Objective.Regularization_Objective_zz_Apply()`

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

        function [val, grad_y] = Time_Instance_Objective(this, y, t)
            error('Time_Instance_Objective() not implemented');
        end

        function [val, grad_z] = Regularization_Objective(this, z)
            error('Regularization_Objective() not implemented');
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this, v, y, t)
            error('Time_Instance_Objective_yy_Apply() not implemented');
        end

        function [Mv] = Regularization_Objective_zz_Apply(this, v, z)
            error('Regularization_Objective_zz_Apply() not implemented');
        end

    end
end
```
