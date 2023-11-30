classdef Constant_Operator < OpInf_Operator
    % Constant operator
    % :math:`\mathcal{C}(\y,\q) = \C \in\R^{n_y}`.

    methods (Access = public)

        function this = Constant_Operator(C)
            arguments
                C (:, :) {mustBeNumeric} = []
            end
            this = this@OpInf_Operator(C);
        end

        function Set_Entries(this, C)
            % Set the entries of the operator.
            %
            % Parameters
            % ----------
            % C
            %   Vector representation of the operator, :math:`\C\in\R^{n_y}`.
            if size(C, 2) ~= 1
                error('Constant_Operator entries must be a single column vector');
            end
            Set_Entries@OpInf_Operator(this, C);
        end

        function [out] = Apply(this, y, ~)
            % Apply the operator to the given state(s) and input(s):
            % :math:`\mathcal{C}(\y,\q) = \C` or
            % :math:`\mathcal{C}(\Y,\Q) = [~\C~~\cdots~~\C~]`.
            %
            % Parameters
            % ----------
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}`,
            %   or an :math:`n_y \times n_t` collection of these.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`
            %   or an :math:`n_q \times n_t` collection of these.
            %
            % Returns
            % -------
            % out : :math:`n_y`-vector or :math:`n_y\times n_t` matrix
            %   Application of the operator, :math:`\C` or :math:`\C\1\trp`.
            out = this.entries(:, ones(1, size(y, 2)));
        end

        function [reduced] = Galerkin(this, Vr, ~)
            % Compute the Galerkin projection of this operator :math:`\C\in\R^{n_y}`
            % with respect to a trial basis :math:`\V_r\in\R^{n_y \times n_y'}`,
            % i.e., :math:`\hat{\C} = \V_r\trp\C \in \R^{n_y'}`.
            %
            % Parameters
            % ----------
            % Vr
            %   Basis matrix :math:`\V_r\in\R^{n_y \times n_y'}` for the trial space.
            %
            % Returns
            % -------
            % reduced : Constant_Operator
            %   Galerkin projection of this operator (a new object).
            reduced = Constant_Operator();
            reduced.Set_Entries(Vr' * this.entries);
        end

    end

    methods (Static, Access = public)

        function [d] = Column_Dimension(~, ~)
            % Column dimension of the operator entries.
            %
            % Parameters
            % ----------
            % n_y
            %   Dimension :math:`n_y` of the ODE state :math:`\y(t)`.
            % n_q
            %   Dimension :math:`n_q` of the input :math:`\q(t)`.
            %
            % Returns
            % -------
            % d : uint8
            %   Operator column dimension :math:`d = 1`.
            d = 1;
        end

        function [block] = Datablock(Y, ~)
            % Construct the data matrix block corresponding to the operator.
            %
            % Since
            %
            % .. math::
            %  \min_{\C}\sum_{j=1}^{n_t}\left\| \C(\y_j, \q_j) \right\|_{2}^{2}
            %  = \min_{\C}\left\| [~\C~~\cdots~~\C~]\right\|_{F}^{2}
            %  = \min_{\C}\left\| \C\1\trp\right\|_{F}^{2},
            %
            % the data block is :math:`\D = \1\trp\in\R^{1 \times n_t}`.
            %
            % Parameters
            % ----------
            % Y
            %   State data :math:`\Y\in\R^{n_y \times n_t}`.
            %   Each column is a single state vector :math:`\y(t_j)`.
            % Q
            %   Input data :math:`\Q\in\R^{n_q \times n_t}`.
            %   Each column is a single input vector :math:`\q(t_j)`.
            %
            % Returns
            % -------
            % block : :math:`1 \times n_t` matrix
            %     Data matrix block of ones, :math:`\1\trp\in\R^{1 \times n_t}`
            block = ones(1, size(Y, 2));
        end

    end

end
