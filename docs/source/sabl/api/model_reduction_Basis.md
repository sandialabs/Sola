# `Basis`

:::{admonition} Abstract Methods
:class: abstract

{class}`Basis` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Basis.Compress()`
- {meth}`Basis.Decompress()`

See the [Inheritance Template](model_reduction.Basis.template) to start a new subclass of {class}`Basis`.
:::

## Class Definition

```{eval-rst}
.. currentmodule:: model_reduction

.. autoclass:: Basis
   :show-inheritance:
   :members:
```

(model_reduction.Basis.template)=
## Inheritance Template

```matlab
classdef My_Basis < Basis

    properties
        n_y
        r
    end

    methods (Access = public)

        function [states_compressed] = Compress(this, states)
            error('Compress() not implemented');
        end

        function [states] = Decompress(this, states_compressed)
            error('Decompress() not implemented');
        end

    end
end
```
