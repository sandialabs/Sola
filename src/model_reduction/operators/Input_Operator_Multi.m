classdef Input_Operator_Multi < OpInf_Operator_Multi
    % Linear input operator for multilithic states.
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
    % transformation :math:`\mathcal{B}(\y,\q)=\B\q,~~\B\in\R^{n_\ell \times n_q}`
    % for some :math:`n_\ell`.

    properties (SetAccess = protected)
        n_q
    end

    methods (Access = public)

        function this = Input_Operator_Multi(n_q, state_dimensions, B)
            arguments
                n_q
                state_dimensions
                B (:, :) {mustBeNumeric} = []
            end
            this = this@OpInf_Operator_Multi(state_dimensions, B);
            this.n_q = n_q;
        end

        function [out] = Apply(this, ~, q)
            out = this.entries * q;
        end

        function [d] = Column_Dimension(this)
            d = this.n_q;
        end

        function [block] = Datablock(this, Y, Q)
            block = this.Get_Substate(this.index, Y);
        end

        function [jac] = Jacobian_y(this, ~, ~)
            jac = zeros(size(this.entries, 1), this.n_y);
        end

        function [jac] = Jacobian_q(this, ~, ~)
            jac = this.entries;
        end

    end
end
