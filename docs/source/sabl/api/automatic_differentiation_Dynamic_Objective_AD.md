# `Dynamic_Objective_AD`

:::{admonition} Abstract Methods
:class: abstract

{class}`Dynamic_Objective_AD` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Dynamic_Objective_AD.g_AD()`
- {meth}`Dynamic_Objective_AD.R_AD()`

Abstract methods of the {class}`Dynamic_Objective` class are implemented by automatic differentiation of the methods listed above.

See the [Inheritance Template](automatic_differentiation.Dynamic_Objective_AD.template) to start a new subclass of {class}`Dynamic_Objective_AD`.
:::

## Class Definition

```{eval-rst}
.. currentmodule:: automatic_differentiation

.. autoclass:: Dynamic_Objective_AD
   :show-inheritance:
   :members:
```

(automatic_differentiation.Dynamic_Objective_AD.template)=
## Inheritance Template

```matlab
classdef My_Objective < Dynamic_Objective_AD

    methods (Access = public)

        function [val] = g_AD(this, y, t)
            error('g_AD() not implemented');
        end

        function [val] = R_AD(this, z)
            error('R_AD() not implemented');
        end

    end
end
```
