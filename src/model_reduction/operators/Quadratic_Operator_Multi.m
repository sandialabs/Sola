%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % .. math::
    %    \mathcal{H}(\y,\q)=\H[\y_j\otimes\y_k],
    %    ~~
    %    \H\in\R^{n_\ell \times n_j n_k},~ j != k.

    properties (SetAccess = protected)
        in_index1
        in_index2
    end

    methods (Access = public)

        function this = Quadratic_Operator_Multi(out_index, in_index1, in_index2, rs, H)
            % Initialize the operator dimensions.
            %
            % Parameters
            % ----------
            % out_index
            %   Integer index :math:`\ell` of the substate that this operator
            %   maps to. For example, for the equation
            %   :math:`\ddt\y_2(t) = \H[\y_1(t) \otimes \y_3(t)]`,
            %   the operation :math:`(\y1, \y_3) \mapsto \H[\y_1\otimes\y_2]`
            %   can be represented with this class with ``out_index = 2``,
            %   ``in_index1 = 1``, and ``in_index2 = 3``.
            % in_index1
            %   Integer index :math:`j` of the first substate that this
            %   operator acts on.
            % in_index2
            %   Integer index :math:`k` of the second substate that this
            %   operator acts on.
            % state_dimensions
            %   Dimensions :math:`n_1,\ldots,n_L` of the substates.
            %   That is, ``state_dimensions(i)`` is the size of :math:`\y_i`.
            % H
            %   (Optional) Entries of the matrix :math:`\H`.
            arguments
                out_index
                in_index1
                in_index2
                rs
                H (:, :) {mustBeNumeric} = []
            end

            this = this@OpInf_Operator_Multi(out_index, rs, H);

            ndims = size(this.state_dimensions, 1);
            if (in_index1 < 1) || (in_index1 > ndims)
                error('in_index1 not aligned with state_dimensions');
            end
            this.in_index1 = in_index1;

            if (in_index2 < 1) || (in_index2 > ndims)
                error('in_index2 not aligned with state_dimensions');
            end
            this.in_index2 = in_index2;
        end

        function Set_Entries(this, entries)
            % Set the entries of the operator.
            if size(entries, 1) ~= this.state_dimensions(this.out_index)
                error('invalid entries size (dimension 1)');
            end
            r2 = this.state_dimensions(this.in_index1) * this.state_dimensions(this.in_index2);
            if size(entries, 2) ~= r2
                error('invalid entries size (dimension 2)');
            end
            this.entries = entries;
        end

        function [out] = Apply(this, y, ~)
            y1 = this.Get_Substate(this.in_index1, y);
            y2 = this.Get_Substate(this.in_index2, y);
            k = size(y, 2);
            out = zeros(this.state_dimensions(this.out_index), k);
            for j = 1:k
                out(:, j) = this.entries * kron(y1(:, j), y2(:, j));
            end
        end

        function [d] = Column_Dimension(this)
            n1 = this.state_dimensions(this.in_index1);
            n2 = this.state_dimensions(this.in_index2);
            d = n1 * n2;
        end

        function [block] = Datablock(this, Y, ~)
            Y1 = this.Get_Substate(this.in_index1, Y);
            Y2 = this.Get_Substate(this.in_index2, Y);
            k = size(Y, 2);
            block = zeros(size(Y1, 1) * size(Y2, 1), k);
            for j = 1:k
                block(:, j) = kron(Y1(:, j), Y2(:, j));
            end
        end

        function [jac] = Jacobian_y(this, y, ~)
            n1 = this.state_dimensions(this.in_index1);
            n2 = this.state_dimensions(this.in_index2);
            first1 = this.state_indices(this.in_index1);
            last1 = this.state_indices(this.in_index1 + 1) - 1;
            first2 = this.state_indices(this.in_index2);
            last2 = this.state_indices(this.in_index2 + 1) - 1;
            y1 = this.Get_Substate(this.in_index1, y);
            y2 = this.Get_Substate(this.in_index2, y);
            jac = zeros(this.state_dimensions(this.out_index), this.n_y);
            jac(:, first1:last1) = this.entries * kron(eye(n1), y2);
            jac(:, first2:last2) = this.entries * kron(y1, eye(n2));
        end

        function [y_out] = Hessian_yy_Apply(this, y_in, ~, ~, lambda)
            num_vecs = size(y_in, 2);
            y_out = zeros(this.n_y, num_vecs);
            subadjoint = this.Get_Substate(this.out_index, lambda);
            for j = 1:num_vecs
                y_out(:, j) = this.Jacobian_y(y_in(:, j))' * subadjoint;
            end
        end

    end
end
