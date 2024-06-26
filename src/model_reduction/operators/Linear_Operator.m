classdef Linear_Operator < OpInf_Operator
    % Linear state operator
    % :math:`\mathcal{A}(\y,\q)=\A\y,~~\A\in\R^{n_y \times n_y}`.

    methods (Access = public)

        function this = Linear_Operator(A)
            arguments
                A (:, :) {mustBeNumeric} = []
            end
            this = this@OpInf_Operator(A);
        end

        function Set_Entries(this, A)
            % Set the entries of the operator.
            %
            % Parameters
            % ----------
            % A
            %   Matrix representation of the operator,
            %   :math:`\A\in\R^{n_y \times n_y}`.
            if size(A, 1) ~= size(A, 2)
                error('Linear_Operator entries must be a square matrix');
            end
            Set_Entries@OpInf_Operator(this, A);
        end

        function [out] = Apply(this, y, ~)
            % Apply the operator to the given state(s):
            % :math:`\mathcal{A}(\y,\q) = \A\y` or
            % :math:`\mathcal{A}(\Y,\Q) = \A\Y`.
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
            %   Application of the operator, :math:`\A\y` or :math:`\A\Y`.
            out = this.entries * y;
        end

        function [jac] = Jacobian_y(this, ~, ~)
            % Construct the state Jacobian of the operator.
            %
            % Parameters
            % ----------
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}` at time :math:`t`.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`.
            %
            % Returns
            % -------
            % jac : :math:`n_y \times n_y` matrix
            %   State Jacobian
            %   :math:`\mathcal{A}_y(\y,\q) = \A\in\R^{n_y \times n_y}`.
            jac = this.entries;
        end

    end

    methods (Static, Access = public)

        function [d] = Column_Dimension(n_y, ~)
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
            %   Operator column dimension :math:`d = n_y`.
            d = n_y;
        end

        function [block] = Datablock(Y, ~)
            % Construct the data matrix block corresponding to the operator.
            %
            % Since
            %
            % .. math::
            %  \min_{\A}\sum_{j=1}^{n_t}\left\| \A(\y_j, \q_j) \right\|_{2}^{2}
            %  = \min_{\A}\left\| \A[~\y_1~~\cdots~~\y_{n_t}~]\right\|_{F}^{2}
            %  = \min_{\A}\left\| \A\Y\right\|_{F}^{2},
            %
            % the data block is
            % :math:`\D = \Y = [~\y_1~~\cdots~~\y_{n_t}~] \in \R^{n_y \times n_t}`.
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
            % block : :math:`n_y \times n_t` matrix
            %     Data matrix block :math:`\D = \Y \in\R^{n_y \times n_t}`.
            block = Y;
        end

    end
end
