# `OpInf_Operator`

:::{admonition} Abstract Methods
:class: abstract

{class}`OpInf_Operator` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`OpInf_Operator.Apply()`
- {meth}`OpInf_Operator.Jacobian()`
- {meth}`OpInf_Operator.Galerkin()`
- {meth}`OpInf_Operator.Column_Dimension()`
- {meth}`OpInf_Operator.Datablock()`

See the [Inheritance Template](model_reduction.operators.OpInf_Operator.template) to start a new subclass of {class}`OpInf_Operator`.
:::

## Class Definition

```{eval-rst}
.. currentmodule:: model_reduction.operators

.. autoclass:: OpInf_Operator
   :show-inheritance:
   :members:
```

(model_reduction.operators.OpInf_Operator.template)=
## Inheritance Template

```matlab
classdef My_Operator < OpInf_Operator

    properties
        entries
    end

    methods (Access = public)

        function [this] = My_Operator(args)
            error('Constructor not implemented');
        end

        function [out] = Evaluate(this, y, q)
            error('Evaluate() not implemented');
        end

        function [jac] = Jacobian(this, y, q)
            error('Jacobian() not implemented');
        end

        function [reduced] = Galerkin(this, Vr, Wr)
            error('Galerkin() not implemented');
        end

    end

    methods (Static, Access = public)

        function [d] = Column_Dimension(n_y, n_t)
            error('Column_Dimension() not implemented');
        end

        function [block] = Datablock(Y, Q)
            error('Datablock() not implemented');
        end

    end
end
```
