# `Objective_AD`

:::{admonition} Abstract Methods
:class: abstract

{class}`Objective_AD` is an abstract class.
Classes that inherit from it must implement the {meth}`Objective_AD.J_AD()` method.
Abstract methods of the {class}`Objective` class are implemented by automatic differentiation of {meth}`Objective_AD.J_AD()`.

See the [Inheritance Template](automatic_differentiation.Objective_AD.template) to start a new subclass of {class}`Objective_AD`.
:::

## Class Definition

```{eval-rst}
.. currentmodule:: automatic_differentiation

.. autoclass:: Objective_AD
   :show-inheritance:
   :members:
```

(automatic_differentiation.Objective_AD.template)=
## Inheritance Template

```matlab
classdef My_Objective < Objective_AD

    methods (Access = public)

        function [val] = J_AD(this, u, z)
            error('J_AD() not implemented');
        end

    end
end
```
