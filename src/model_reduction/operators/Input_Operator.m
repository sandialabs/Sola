classdef Input_Operator < OpInf_Operator
    % Linear input operator
    % :math:`\mathcal{B}(\y,\q) = \B\q,~~\B\in\R^{n_y \times n_q}`.

    properties (Dependent)
        n_q             % Dimension :math:`n_q` of the input :math:`\q(t)`.
    end

    methods

        function [n_q] = get.n_q(this)
            n_q = size(this.entries, 2);
        end

    end

    methods (Access = public)

        function this = Input_Operator(B)
            arguments
                B (:, :) {mustBeNumeric} = []
            end
            this = this@OpInf_Operator(B);
        end

        function Set_Entries(this, B)
            % Set the entries of the operator.
            %
            % Parameters
            % ----------
            % B
            %   Matrix representation of the operator,
            %   :math:`\B\in\R^{n_y \times n_q}`.
            Set_Entries@OpInf_Operator(this, B);
        end

        function [out] = Apply(this, ~, q)
            % Apply the operator to the given state(s) and input(s):
            % :math:`\mathcal{B}(\y,\q) = \B\q` or
            % :math:`\mathcal{B}(\Y,\Q) = \B\Q`.
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
            %   Application of the operator, :math:`\B\q` or :math:`\B\Q`.
            out = this.entries * q;
        end

        function [jac] = Jacobian_q(this, ~, ~)
            % Construct the input Jacobian of the operator.
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
            % jac : :math:`n_y \times n_q` matrix
            %   Input Jacobian
            %   :math:`\mathcal{B}_{q}(\y,\q) = \B\in\R^{n_y \times n_q}`.
            jac = this.entries;
        end

        function Finite_Difference_Check(this, n_t)
            % Check consistency between :meth:`Apply()` and :meth:`Datablock()`
            % and do finite difference checks for :meth:`Jacobian_q()`.
            arguments
                this
                n_t uint8 = 30
            end

            % Ensure consistency with Apply() and Datablock().
            Y = randn(this.n_y, n_t);
            Q = randn(this.n_q, n_t);
            D = this.Datablock(Y, Q);
            assert(size(D, 1) == this.Column_Dimension(this.n_y, this.n_q));
            appliedtoQ = this.Apply(Y, Q);
            fromDblock = this.entries * D;
            assert(size(appliedtoQ, 1) == size(fromDblock, 1));
            assert(size(appliedtoQ, 2) == n_t);
            assert(size(fromDblock, 2) == n_t);
            assert(norm(appliedtoQ - fromDblock) < 1e-12);

            % TODO: finite difference check for Jacobian_q().

            % Ensure Galerkin() gives a dimension reduction.
            if this.n_y > 1
                r = floor(this.n_y / 2);
                Vr = randn(this.n_y, r);
                op = this.Galerkin(Vr, Vr);
                assert(size(op.entries, 1) == r);
                assert(size(op.entries, 2) == this.n_q);
                assert(isa(op, class(this)));
            end
        end

    end

    methods (Static, Access = public)

        function [d] = Column_Dimension(~, n_q)
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
            %   Operator column dimension :math:`d = n_q`.
            d = n_q;
        end

        function [block] = Datablock(~, Q)
            % Construct the data matrix block corresponding to the operator.
            %
            % Since
            %
            % .. math::
            %  \min_{\B}\sum_{j=1}^{n_t}\left\| \B(\y_j, \q_j) \right\|_{2}^{2}
            %  = \min_{\B}\left\| \B[~\q_1~~\cdots~~\q_{n_t}~]\right\|_{F}^{2}
            %  = \min_{\B}\left\| \B\Q\right\|_{F}^{2},
            %
            % the data block is
            % :math:`\D = \Q = [~\q_1~~\cdots~~\q_{n_t}~] \in \R^{n_q \times n_t}`.
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
            % block : :math:`n_q \times n_t` matrix
            %     Data matrix block :math:`\D = \Q \in\R^{n_q \times n_t}`.
            block = Q;
        end

    end
end
