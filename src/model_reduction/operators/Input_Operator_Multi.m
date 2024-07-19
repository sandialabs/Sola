classdef Input_Operator_Multi < OpInf_Operator_Multi
    % Linear input operator for multilithic states.
    %
    % Write the state vector :math:`\y\in\R^{n_y}` as
    %
    % .. math::
    %    \y = \left[\begin{array}{c}
    %    \y_1 \\ \vdots \\ \y_L
    %    \end{array}\right],
    %
    % where each :math:`\y_\ell\in\R^{n_\ell}` with
    % :math:`\sum_{\ell=1}^{L}n_\ell = n_y`.
    % This class represents a linear transformation of the input to the space
    % of one of the substates, i.e.,
    % .. math:
    %    \mathcal{B}(\y,\q) = \B\q,
    %    ~~
    %    \B\in\R^{n_\ell \times n_q}.

    properties (SetAccess = protected)
        n_q
    end

    methods (Access = public)

        function this = Input_Operator_Multi(out_index, n_q, state_dimensions, B)
            % Initialize the operator dimensions.
            %
            % Parameters
            % ----------
            % out_index
            %   Integer index :math:`\ell` of the substate that this operator
            %   maps to. For example, for :math:`\ddt\y_1(t) = \A\y_2(t)`, the
            %   operation :math:`\y_2 \mapsto \A\y_2` can be represented
            %   with this class with ``out_index = 1`` and ``in_index = 2``.
            % in_index
            %   Integer index of the substate that this operator acts on.
            % state_dimensions
            %   Dimensions :math:`n_1,\ldots,n_L` of the substates.
            %   That is, ``state_dimensions(i)`` is the size of :math:`\y_i`.
            % A
            %   (Optional) Entries of the matrix :math:`\A_{\ell, k}`.
            arguments
                out_index
                n_q
                state_dimensions
                B (:, :) {mustBeNumeric} = []
            end
            this = this@OpInf_Operator_Multi(out_index, state_dimensions, B);
            this.n_q = n_q;
        end

        function Set_Entries(this, entries)
            % Set the entries of the operator.
            if size(entries, 1) ~= this.state_dimensions(this.out_index)
                error('invalid entries size (dimension 1)');
            elseif size(entries, 2) ~= this.n_q
                error('invalid entries size (dimension 2)');
            end
            this.entries = entries;
        end

        function [out] = Apply(this, ~, q)
            out = this.entries * q;
        end

        function [d] = Column_Dimension(this)
            d = this.n_q;
        end

        function [block] = Datablock(~, ~, Q)
            block = Q;
        end

        function [jac] = Jacobian_q(this, ~, ~)
            jac = this.entries;
        end

    end
end
