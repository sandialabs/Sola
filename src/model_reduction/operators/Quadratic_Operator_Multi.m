classdef Quadratic_Operator_Multi < OpInf_Operator_Multi
    % Quadratic state operator for multilithic states.
    %
    % Writing the state vector :math:`\y\in\R^{n_y}` as
    %
    % .. math::
    %    \y = \left[\begin{array}{c}
    %    \y_1 \\ \vdots \\ \y_L
    %    \end{array}\right]
    %
    % where each :math:`\y_\ell\in\R^{n_\ell}` with
    % :math:`\sum_{\ell=1}^{L}n_\ell = n_y`, this class represents a quadratic
    % function of two different :math:`\y_\ell`, i.e.,
    % :math:`\mathcal{H}(\y,\q)=\H[\y_k\otimes\y_\ell]`,
    % :math:`\H\in\R^{n_\ell \times n_j n_k}`.

    properties (SetAccess = protected)
        index1
        index2
    end

    methods (Access = public)

        function this = Quadratic_Operator_Multi(index1, index2, rs, H)
            arguments
                index1
                index2
                rs
                H (:, :) {mustBeNumeric} = []
            end
            this = this@OpInf_Operator_Multi(rs, H);
            if index1 == index2
                error('index2 must not equal index1');
            end
            this.index1 = index1;
            this.index2 = index2;
        end

        function [out] = Apply(this, y, ~)
            y1 = this.Get_Substate(this.index1, y);
            y2 = this.Get_Substate(this.index2, y);
            out = this.entries * kron(y1, y2);
        end

        function [d] = Column_Dimension(this)
            n1 = this.state_dimensions(this.index1);
            n2 = this.state_dimensions(this.index2);
            d = n1 * n2;
        end

        function [block] = Datablock(this, Y, ~)
            Y1 = this.Get_Substate(this.index1, Y);
            Y2 = this.Get_Substate(this.index2, Y);
            k = size(Y, 2);
            block = zeros(size(Y1, 1) * size(Y2, 1), k);
            for j = 1:k
                block(:, k) = kron(Y1(:, k), Y2(:, k));
            end
        end

        function [jac] = Jacobian_y(this, y, ~)
            n_ell = size(this.entries, 1);
            n1 = this.state_dimensions(this.index1);
            n2 = this.state_dimensions(this.index2);
            first1 = this.state_indices(this.index1);
            last1 = this.state_indices(this.index1 + 1) - 1;
            first2 = this.state_indices(this.index2);
            last2 = this.state_indices(this.index2 + 1) - 1;
            y1 = this.Get_Substate(this.index1, y);
            y2 = this.Get_Substate(this.index2, y);
            jac = zeros(n_ell, this.n_y);
            jac(:, first1:last1) = this.entries * kron(eye(n1), y2);
            jac(:, first2:last2) = this.entries * kron(y1, eye(n2));
        end

        function [Mv] = Hessian_yy_Apply(this, v, ~, ~, lambda)
            Mv = zeros(size(this.entries, 1), size(v, 2));
            for j = 1:num_vecs
                Mv(:, j) = this.Jacobian_y(v(:, j))' * lambda;
                % TODO: verify the transpose.
            end
        end

    end
end
