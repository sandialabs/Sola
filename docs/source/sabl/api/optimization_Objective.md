# `Objective`

::::{admonition} Abstract Methods
:class: abstract

{class}`Objective` is an abstract class.
Classes that inherit from it must implement the following methods.

- {meth}`Objective.J()`
- {meth}`Objective.J_uu_Apply()`
- {meth}`Objective.J_uz_Apply()`
- {meth}`Objective.J_zu_Apply()`
- {meth}`Objective.J_zz_Apply()`

:::{danger}
Because of how MATLAB's [`fminunc()`](https://www.mathworks.com/help/optim/ug/fminunc.html) is designed, the `J_xx_Apply()` functions must be implemented in a _vectorized_ fashion, i.e., assuming that `v` is a matrix where each column is a test direction.
:::

See the [Inheritance Template](optimization.Objective.template) to start a new subclass of {class}`Objective`.
::::

## Class Definition

```{eval-rst}
.. currentmodule:: optimization

.. autoclass:: Objective
   :show-inheritance:
   :members:
```

(optimization.Objective.template)=
## Inheritance Template

```matlab
classdef My_Objective < Objective

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            error('J() not implemented');
        end

        function [u_out] = J_uu_Apply(this, u_in, u, z)
            error('J_uu_Apply() not implemented');
        end

        function [u_out] = J_uz_Apply(this, z_in, u, z)
            error('J_uz_Apply() not implemented');
        end

        function [z_out] = J_zu_Apply(this, u_in, u, z)
            error('J_zu_Apply() not implemented');
        end

        function [z_out] = J_zz_Apply(this, z_in, u, z)
            error('J_zz_Apply() not implemented');
        end

    end
end
```
