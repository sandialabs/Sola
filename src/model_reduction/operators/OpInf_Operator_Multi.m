classdef OpInf_Operator_Multi < handle
    % Base class for operators acting on multilithic states.
    %
    % Writing the state vector :math:`\y\in\R^{n_y}` as
    %
    % .. math::
    %    \y = \left[\begin{array}{c}
    %    \y_1 \\ \vdots \\ \y_L
    %    \end{array}\right]
    %
    % where each :math:`\y_\ell\in\R^{n_\ell}` with
    % :math:`\sum_{\ell=1}^{L}n_\ell = n_y`, this class represents operators
    % :math:`\mathcal{F}:\R^{n_y}\times\R^{n_q}\to\R^{n_\ell}`
    % for some :math:`\ell`.

    properties (SetAccess = protected)
        state_dimensions  % Dimensions :math:`n_1,\ldots,n_L` of the substates.
        n_y               % Dimension :math:`n_y` of the ODE state :math:`\y(t)`.
        entries           % Matrix representation of the operator.
        out_index         % Index :math:`\ell` of the substate that this operator maps to.
    end

    properties (Access = protected)
        state_indices
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
        %   Application of the operator to ``(y,q)``, either an
        %   :math:`n_\ell`-vector or an :math:`n_\ell\times n_t` matrix,
        %   where :math:`n_\ell` is the dimension of one of the substates.

        [d] = Column_Dimension(this)
        % Column dimension of the operator entries.
        %
        % Returns
        % -------
        % d : uint8
        %   Column dimension of the operator entries.

        [block] = Datablock(this, Y, Q)
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

    end

    methods (Access = public)

        %% Constructor and Initializer.

        function this = OpInf_Operator_Multi(out_index, state_dimensions, entries)
            % Initialize the operator and (optionally) set its entries.
            %
            % Parameters
            % ----------
            % out_index
            %   Integer index :math:`\ell` of the substate that this operator
            %   maps to.
            % state_dimensions
            %   Dimensions :math:`n_1,\ldots,n_L` of the substates.
            % entries
            %   (Optional) Operator entries.
            arguments
                out_index
                state_dimensions
                entries (:, :) {mustBeNumeric} = []
            end

            this.state_dimensions = reshape(state_dimensions, [], 1);
            this.state_indices = cumsum([1; this.state_dimensions]);
            this.n_y = sum(this.state_dimensions);
            this.entries = [];
            nstates = size(this.state_dimensions, 1);

            if (out_index < 1) || (out_index > nstates)
                error('out_index not aligned with state_dimensions');
            end
            this.out_index = out_index;

            if size(entries, 1) > 0
                this.Set_Entries(entries);
            end
        end

        function Set_Entries(this, entries)
            % Set the entries of the operator.
            this.entries = entries;
        end

        function [substate] = Get_Substate(this, index, state)
            % Extract one of the substates.
            %
            % Parameters
            % ----------
            % index
            %   Which substate to extract.
            % state
            %   Differential equation state :math:`\y(t)\in\R^{n_y}` to extract from.
            %
            % Returns
            % -------
            % substate : :math:`n_i` vector
            first = this.state_indices(index);
            last = this.state_indices(index + 1) - 1;
            substate = state(first:last, :);
        end

        %% Operator derivatives.
        % If these are not implemented by child classes,
        % they are assumed to be zero.

        function [jac] = Jacobian_y(this, y, q)
            % Construct the partial state Jacobian of the operator.
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
            % jac : :math:`n_\ell \times n_y` matrix
            %   Partial state Jacobian
            %   :math:`\mathcal{F}_{y}(\y,\q)\in\R^{n_\ell \times n_y}`.
            jac = 0;
        end

        function [jac] = Jacobian_q(this, y, q)
            % Construct the partial input Jacobian of the operator.
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
            % jac : :math:`n_\ell \times n_q` matrix
            %   Partial input Jacobian
            %   :math:`\mathcal{F}_{q}(\y,\q)\in\R^{n_\ell \times n_q}`.
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

        %% Verification.

        function Finite_Difference_Check(this, n_t, hasinputs)
            % Check consistency between :meth:`Apply()` and :meth:`Datablock()`
            % and do a finite difference checks for :meth:`Jacobian_y()`
            % and :meth:`Hessian_yy()`.
            arguments
                this
                n_t uint8 = 30
                hasinputs = false
            end

            % Ensure consistency with Apply() and Datablock().
            Y = randn(this.n_y, n_t);
            Q = randn(1, n_t);
            if hasinputs
                Q = randn(this.n_q, n_t);
            end
            D = this.Datablock(Y, Q);
            assert(size(D, 1) == this.Column_Dimension());
            assert(size(D, 2) == n_t);
            appliedtoY = this.Apply(Y, Q);
            fromDblock = this.entries * D;
            assert(size(appliedtoY, 1) == size(fromDblock, 1));
            assert(size(appliedtoY, 2) == n_t);
            assert(size(fromDblock, 2) == n_t);
            assert(norm(appliedtoY - fromDblock) < 1e-12);

            % Finite difference check for Jacobian_y().
            y = Y(:, 1);
            q = Q(:, 1);
            v = randn(this.n_y, 1);
            v = v / norm(v);
            h = 10.^(-2:-1:-9);
            p = length(h);
            out = this.Apply(y, q);
            jac = this.Jacobian_y(y, q);
            jac_v = jac * v;
            if norm(jac_v) > 0
                disp(' ');
                disp(['Finite difference check for ', class(this), '.Jacobian_y()']);
                for k = 1:p
                    out_k = this.Apply(y + h(k) * v);
                    diff = norm(jac_v - (out_k - out) / h(k)) / norm(jac_v);
                    fprintf('h = %.2e error = %.3e\n', h(k), diff);
                end
                disp(' ');
            end

            % Finite difference check for Jacobian_q().
            if hasinputs
                vq = randn(this.n_q, 1);
                vq = vq / norm(vq);
                jacq = this.Jacobian_q(y, q);
                jacq_v = jacq * vq;
                if norm(jacq_v) > 0
                    disp(' ');
                    disp(['Finite difference check for ', class(this), '.Jacobian_q()']);
                    for k = 1:p
                        out_k = this.Apply(y, q + h(k) * vq);
                        diff = norm(jacq_v - (out_k - out) / h(k)) / norm(jacq_v);
                        fprintf('h = %.2e error = %.3e\n', h(k), diff);
                    end
                end
            end

            % Finite difference check for Hessian_yy_Apply().
            lambda = randn(this.n_y, 1);
            Mv = this.Hessian_yy_Apply(v, y, q, lambda);
            if norm(Mv) > 0
                disp(' ');
                disp(['Finite difference check for ', class(this), '.Hessian_yy_Apply()']);
                first = this.state_indices(this.out_index);
                last = this.state_indices(this.out_index + 1) - 1;
                fulljac = zeros(this.n_y, this.n_y);
                fulljac(first:last, :) = jac;
                for k = 1:p
                    out_k = zeros(this.n_y, this.n_y);
                    out_k(first:last, :) = this.Jacobian_y(y + h(k) * v);
                    fd_Mv = (out_k' * lambda - fulljac' * lambda) / h(k);
                    diff = norm(fd_Mv - Mv) / norm(Mv);
                    fprintf('h = %.2e error = %.3e\n', h(k), diff);
                end
                disp(' ');
            end

            if ~hasinputs
                return
            end

            % Finite difference check for Hessian_qq_Apply().
            Mv = this.Hessian_qq_Apply(vq, y, q, lambda);
            if norm(Mv) > 0
                disp(' ');
                disp(['Finite difference check for ', class(this), '.Hessian_qq_Apply()']);
                first = this.state_indices(this.out_index);
                last = this.state_indices(this.out_index + 1) - 1;
                fulljac = zeros(this.n_y, this.n_q);
                fulljac(first:last, :) = jacq;
                for k = 1:p
                    out_k = zeros(this.n_y, this.n_q);
                    out_k(first:last, :) = this.Jacobian_q(y, q + h(k) * vq);
                    fd_Mv = (out_k' * lambda - fulljac' * lambda) / h(k);
                    diff = norm(fd_Mv - Mv) / norm(Mv);
                    fprintf('h = %.2e error = %.3e\n', h(k), diff);
                end
                disp(' ');
            end

        end

    end
end
