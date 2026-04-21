%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Linear_Operator_Multi < OpInf_Operator_Multi
    % Linear state operator for multilithic states.
    %
    % Write the state vector :math:`\y\in\R^{n_y}` as substates
    %
    % .. math::
    %    \y = \left[\begin{array}{c}
    %    \y_1 \\ \vdots \\ \y_L
    %    \end{array}\right],
    %
    % where each :math:`\y_\ell\in\R^{n_\ell}` with
    % :math:`\sum_{\ell=1}^{L}n_\ell = n_y`.
    % This class represents a linear transformation from the dimensions of one
    % substate on another, i.e.,
    % .. math::
    %    \mathcal{A}_{\ell,k}(\y,\q)=\A\y_k,
    %    ~~
    %    \A\in\R^{n_\ell \times n_k}.

    properties (SetAccess = protected)
        in_index
    end

    methods (Access = public)

        function this = Linear_Operator_Multi(out_index, in_index, state_dimensions, A)
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
                in_index
                state_dimensions
                A (:, :) {mustBeNumeric} = []
            end

            this = this@OpInf_Operator_Multi(out_index, state_dimensions, A);

            if (in_index < 1) || (in_index > size(this.state_dimensions, 1))
                error('in_index not aligned with state_dimensions');
            end
            this.in_index = in_index;
        end

        function Set_Entries(this, entries)
            % Set the entries of the operator.
            if size(entries, 1) ~= this.state_dimensions(this.out_index)
                error('invalid entries size (dimension 1)');
            elseif size(entries, 2) ~= this.state_dimensions(this.in_index)
                error('invalid entries size (dimension 2)');
            end
            this.entries = entries;
        end

        function [out] = Apply(this, y, ~)
            out = this.entries * this.Get_Substate(this.in_index, y);
        end

        function [d] = Column_Dimension(this)
            d = this.state_dimensions(this.in_index);
        end

        function [block] = Datablock(this, Y, ~)
            block = this.Get_Substate(this.in_index, Y);
        end

        function [jac] = Jacobian_y(this, ~, ~)
            jac = zeros(this.state_dimensions(this.out_index), this.n_y);
            first = this.state_indices(this.in_index);
            last = this.state_indices(this.in_index + 1) - 1;
            jac(:, first:last) = this.entries;
        end

    end
end
