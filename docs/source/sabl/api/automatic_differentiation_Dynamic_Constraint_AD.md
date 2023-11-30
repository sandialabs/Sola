# `Dynamic_Constraint_AD`

:::{admonition} Abstract Methods
:class: abstract

{class}`Dynamic_Constraint_AD` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Dynamic_Constraint_AD.Time_Instance_RHS_AD()`
- {meth}`Dynamic_Constraint_AD.Initial_Condition_AD()`

Abstract methods of the {class}`Dynamic_Constraint` class are implemented by automatic differentiation of the methods listed above.

See the [Inheritance Template](automatic_differentiation.Dynamic_Constraint_AD.template) to start a new subclass of {class}`Dynamic_Constraint_AD`.
:::

## Class Definition

```{eval-rst}
.. currentmodule:: automatic_differentiation

.. autoclass:: Dynamic_Constraint_AD
   :show-inheritance:
   :members:
```

(automatic_differentiation.Dynamic_Constraint_AD.template)=
## Inheritance Template

```matlab
classdef My_Constraint < Dynamic_Constraint_AD

    methods (Access = public)

        function [f] = Time_Instance_RHS_AD(this, y, z, t)
            error('Time_Instance_RHS_AD() not implemented');
        end

        function [h] = Initial_Condition_AD(this, z)
            error('Initial_Condition_AD() not implemented');
        end

    end
end
```
