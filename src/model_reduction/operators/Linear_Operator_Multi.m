classdef Linear_Operator_Multi < OpInf_Operator_Multi
    % Linear state operator for multilithic states.
    %
    % Writing the state vector :math:`\y\in\R^{n_y}` as
    %
    % .. math::
    %    \y = \left[\begin{array}{c}
    %    \y_1 \\ \vdots \\ \y_L
    %    \end{array}\right]
    %
    % where each :math:`\y_\ell\in\R^{n_\ell}` with
    % :math:`\sum_{\ell=1}^{L}n_\ell = n_y`, this class represents a linear
    % transformation of on one of the :math:`\y_\ell`, i.e.,
    % :math:`\mathcal{A}_{k,\ell}(\y,\q)=\A\y_\ell,~~\A\in\R^{n_k \times n_\ell}`.

    properties (SetAccess = protected)
        index
    end

    methods (Access = public)

        function this = Linear_Operator_Multi(index, state_dimensions, A)
            arguments
                index
                state_dimensions
                A (:, :) {mustBeNumeric} = []
            end
            this = this@OpInf_Operator_Multi(state_dimensions, A);
            this.index = index;
        end

        function [out] = Apply(this, y, ~)
            out = this.entries * this.Get_Substate(this.index, y);
        end

        function [d] = Column_Dimension(this)
            d = this.state_dimensions(this.index);
        end

        function [block] = Datablock(this, Y, ~)
            block = this.Get_Substate(this.index, Y);
        end

        function [jac] = Jacobian_y(this, ~, ~)
            n_ell = this.state_dimensions(this.index);
            jac = zeros(n_ell, this.n_y);
            first = this.state_indices(this.index);
            last = this.state_indices(this.index + 1) - 1;
            jac(:, first:last) = this.entries;
        end

    end
end
