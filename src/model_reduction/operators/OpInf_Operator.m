classdef OpInf_Operator < handle
    % Base class for operators to be used in operator inference ROMs.

    properties
        entries         % Matrix representation of the operator.
    end

    properties (Dependent)
        n_y             % Dimension :math:`n_y` of the ODE state :math:`\y(t)`.
    end

    methods

        function [n_y] = get.n_y(this)
            n_y = size(this.entries, 1);
        end

    end

    methods (Abstract, Access = public)

        [out] = Apply(this, y, q)
        % Apply the operator to the given state and input.
        %
        % Parameters
        % ----------
        % y
        %   Differential equation state :math:`\y(t)\in\R^{n_y}`,
        %   or an :math:`n_y \times n_t` collection of these.
        % q
        %   Input :math:`\q(t)\in\R^{n_q}`,
        %   or an :math:`n_q \times n_t` collection of these.
        %
        % Returns
        % -------
        % out : vector or matrix
        %   Application of the operator to ``(y,q)``,
        %   either an :math:`n_y`-vector or an :math:`n_y\times n_t` matrix.

    end

    methods (Access = public)

        %% Constructor and Initializer.

        function this = OpInf_Operator(entries)
            % Initialize the operator and (optionally) set its entries.
            %
            % Parameters
            % ----------
            % entries
            %   (Optional) Operator entries.
            arguments
                entries (:, :) {mustBeNumeric} = []
            end

            if size(entries, 1) > 0
                this.Set_Entries(entries);
            end
        end

        function Set_Entries(this, entries)
            % Set the entries of the operator.
            this.entries = entries;
        end

        %% Operator derivatives.
        % If these are not implemented by child classes,
        % they are assumed to be zero.

        function [jac] = Jacobian_y(this, y, q)
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
            %   State Jacobian :math:`\mathcal{F}_{y}(\y,\q)\in\R^{n_y \times n_y}`.
            jac = 0;
        end

        function [jac] = Jacobian_q(this, y, q)
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
            %   Input Jacobian :math:`\mathcal{F}_{q}(\y,\q)\in\R^{n_y \times n_q}`.
            jac = 0;
        end

        function [Mv] = Hessian_yy_Apply(this, v, y, q, lambda)
            % Compute the action of the :math:`y,y` Hessian of the operator.
            %
            % Parameters
            % ----------
            % v
            %   Search direction :math:`\v\in\R^{n_y}`.
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}` at time :math:`t`.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`.
            % lambda
            %   Adjoint state :math:`\lambda(t)\in\R^{n_y}` at time :math:`t`.
            %
            % Returns
            % -------
            % Mv : :math:`n_y`-vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp\mathcal{F}_{y,y}(\y,\q)\v\in\R^{n_y}`.
            Mv = 0;
        end

        function [Mv] = Hessian_yq_Apply(this, v, y, q, lambda)
            % Compute the action of the :math:`y,q` Hessian of the operator.
            %
            % Parameters
            % ----------
            % v
            %   Search direction :math:`\v\in\R^{n_q}`.
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}` at time :math:`t`.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`.
            % lambda
            %   Adjoint state :math:`\lambda(t)\in\R^{n_y}` at time :math:`t`.
            %
            % Returns
            % -------
            % Mv : :math:`n_y`-vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp\mathcal{F}_{y,q}(\y,\q)\v\in\R^{n_y}`.
            Mv = 0;
        end

        function [Mv] = Hessian_qy_Apply(this, v, y, q, lambda)
            % Compute the action of the :math:`q,y` Hessian of the operator.
            %
            % Parameters
            % ----------
            % v
            %   Search direction :math:`\v\in\R^{n_y}`.
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}` at time :math:`t`.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`.
            % lambda
            %   Adjoint state :math:`\lambda(t)\in\R^{n_y}` at time :math:`t`.
            %
            % Returns
            % -------
            % Mv : :math:`n_q`-vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp\mathcal{F}_{q,y}(\y,\q)\v\in\R^{n_q}`.
            Mv = 0;
        end

        function [Mv] = Hessian_qq_Apply(this, v, y, q, lambda)
            % Compute the action of the :math:`q,q` Hessian of the operator.
            %
            % Parameters
            % ----------
            % v
            %   Search direction :math:`\v\in\R^{n_q}`.
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}` at time :math:`t`.
            % q
            %   Input :math:`\q(t)\in\R^{n_q}` at time :math:`t`.
            % lambda
            %   Adjoint state :math:`\lambda(t)\in\R^{n_y}` at time :math:`t`.
            %
            % Returns
            % -------
            % Mv : :math:`n_q`-vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp\mathcal{F}_{q,q}(\y,\q)\v\in\R^{n_q}`.
            Mv = 0;
        end

        %% Galerkin projection (dimensionality reduction).

        function [reduced] = Galerkin(this, Vr, Wr)
            % Compute the Galerkin projection of this operator with respect to
            % :math:`n_y'`-dimensional trial and test bases.
            %
            % Parameters
            % ----------
            % Vr
            %   Basis matrix :math:`\V_r\in\R^{n_y \times n_y'}` for the trial space.
            % Wr
            %   Basis matrix :math:`\W_r\in\R^{n_y \times n_y'}` for the test space.
            %
            % Returns
            % -------
            % reduced : OpInf_Operator
            %   Galerkin projection of this operator (a new object).
            error('Galerkin() not implemented');
        end

        %% Verification.

        function Finite_Difference_Check(this, n_t)
            % Check consistency between :meth:`Apply()` and :meth:`Datablock()`
            % and do a finite difference checks for :meth:`Jacobian_y()`
            % and :meth:`Hessian_yy()`.
            arguments
                this
                n_t uint8 = 30
            end

            % Ensure consistency with Apply() and Datablock().
            Y = randn(this.n_y, n_t);
            D = this.Datablock(Y);
            assert(size(D, 1) == this.Column_Dimension(this.n_y));
            appliedtoY = this.Apply(Y);
            fromDblock = this.entries * D;
            assert(size(appliedtoY, 1) == size(fromDblock, 1));
            assert(size(appliedtoY, 2) == n_t);
            assert(size(fromDblock, 2) == n_t);
            assert(norm(appliedtoY - fromDblock) < 1e-12);

            % Finite difference check for Jacobian_y().
            y = Y(:, 1);
            v = randn(this.n_y, 1);
            v = v / norm(v);
            h = 10.^(-2:-1:-9);
            p = length(h);
            out = this.Apply(y);
            jac = this.Jacobian_y(y, y);
            jac_v = jac * v;
            if norm(jac_v) > 0
                disp(' ');
                disp(['Finite difference check for ', class(this), '.Jacobian_y()']);
                for k = 1:p
                    out_k = this.Apply(y + h(k) * v);
                    diff = norm(jac_v - (out_k - out) / h(k)) / norm(jac_v);
                    disp(['h = ', num2str(h(k)), ' error = ', num2str(diff)]);
                end
                disp(' ');
            end

            % Finite difference check for Hessian_yy().
            lambda = randn(this.n_y, 1);
            Mv = this.Hessian_yy_Apply(v, y, y, lambda);
            diffs_yy = zeros(p, 1);
            if norm(Mv) > 0
                disp(' ');
                disp(['Finite difference check for ', class(this), '.Hessian_yy_Apply()']);
                for k = 1:p
                    out_k = this.Jacobian_y(y + h(k) * v);
                    fd_Mv = (out_k' * lambda - jac' * lambda) / h(k);
                    diff = norm(fd_Mv - Mv) / norm(Mv);
                    disp(['h = ', num2str(h(k)), ' error = ', num2str(diff)]);
                end
                disp(' ');
            end

            % Ensure Galerkin() gives a dimension reduction.
            if this.n_y > 1
                r = floor(this.n_y / 2);
                Vr = randn(this.n_y, r);
                op = this.Galerkin(Vr, Vr);
                assert(op.n_y == r);
                assert(isa(op, class(this)));
            end
        end

    end

    methods (Static, Access = public)

        function [d] = Column_Dimension(n_y, n_q)
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
            %   Column dimension of the operator entries.
            error('Column_Dimension() not implemented');
        end

        function [block] = Datablock(Y, Q)
            % Construct the data matrix block corresponding to the operator.
            %
            % Let :math:`\A(\y,\q)` represent the operator acting on a pair of
            % state and input vectors. The data matrix block is the matrix
            % :math:`\D` containing the state and input data such that
            %
            % .. math::
            %  \min_{\A}\sum_{j=1}^{n_t}\left\| \A(\y_j, \q_j)\right\|_{2}^{2},
            %
            % can be written equivalently as
            %
            % .. math:: \min_{\X}\left\| \X\D \right\|_{F}^{2}
            %
            % where :math:`\X` are the operator entries.
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
            % block : :math:`d \times n_t` matrix
            %     Data matrix block :math:`\D\in\R^{d \times n_t}`.
            error('Datablock() not implemented');
        end

    end
end
