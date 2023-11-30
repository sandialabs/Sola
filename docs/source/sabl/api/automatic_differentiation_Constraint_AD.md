# `Constraint_AD`

:::{admonition} Abstract Methods
:class: abstract

{class}`Constraint_AD` is an abstract class.
Classes that inherit from it must implement the {meth}`Constraint_AD.c_AD()` method.
Abstract methods of the {class}`Constraint` class are implemented by automatic differentiation of {meth}`Constraint_AD.c_AD()`.

See the [Inheritance Template](automatic_differentiation.Constraint_AD.template) to start a new subclass of {class}`Constraint_AD`.
:::

## Class Definition

```{eval-rst}
.. currentmodule:: automatic_differentiation

.. autoclass:: Constraint_AD
   :show-inheritance:
   :members:
```

(automatic_differentiation.Constraint_AD.template)=
## Inheritance Template

```matlab
classdef My_Constraint < Constraint_AD

    methods (Access = public)

        function [c] = c_AD(this, u, z)
            error('c_AD() not implemented');
        end

    end
end
```
