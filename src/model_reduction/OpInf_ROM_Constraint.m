classdef OpInf_ROM_Constraint < Dynamic_Constraint
    % Operator Inference Reduced-order Model.
    %
    % This class represents systems of ordinary differential equations (ODEs)
    % with the polynomial
    %
    % .. math:: \ddt\y(t) = \c + \A\y(t) + \H[\y(t)\otimes\y(t)] + \B\q(t),
    %
    % where
    %
    % * :math:`\y(t)\in\R^{n_y}` is the ODE state,
    % * :math:`\q(t)\in\R^{n_q}` is the time-dependent (state-independent) input, and
    % * :math:`\c,\A,\H,\B` are the system "operators".
    %
    % The system operators are determined by solving a regression problem,
    %
    % .. math::
    %  \min_{\c,\A,\H,\B}\sum_{j=2}^{n_t}\left\|
    %  \c + \A\y_j + \H[\y_j\otimes\y_j] + \B\q_j - \dot{\y}_j
    %  \right\|_{2}^{2}
    %  + \left\|\bfGamma[~\c~~\A~~\H~~\B~]\trp\right\|_{F}^{2},
    %
    % where
    %
    % * :math:`\y_j\in\R^{n_y}` is a snapshot :math:`\y(t_j)` of the ODE state at time :math:`t_j`,
    % * :math:`\q_j\in\R^{n_q}` is the input :math:`\q(t_j)` at time :math:`t_j`,
    % * :math:`\dot{\y}_j\in\R^{n_y}` is the time derivative of the ODE state at time :math:`t_j`,
    %   i.e., :math:`\ddt\y(t)\big|_{t = t_j}`,
    % * :math:`n_t \in \N` is the number of snapshots, and
    % * :math:`\bfGamma\in\R^{d \times d}` is a Tikhonov regularization matrix.
    %
    % This class inherits from :class:`Dynamic_Constraint`, hence the model is integrated in time
    % using the implicit Euler method. In this context, the optimization state and control are
    % given by
    %
    % .. math::
    %  \u = \left(\begin{array}{c}\y_1 \\ \vdots \\ \y_{n_t}\end{array}\right)\in\R^{n_u},
    %  \qquad
    %  \z = \left(\begin{array}{c}\q_2 \\ \vdots \\ \q_{n_t}\end{array}\right)\in\R^{n_z},
    %
    % with :math:`n_u = n_y n_t` and :math:`n_z = n_q (n_t - 1)`.
    % Note that the state and the control share a time mesh, which is important for the definition
    % of the regression problem. However, the control is not measured at the initial condition
    % because the system is solved with a one-step implicit time integration scheme.
    %
    % The user indicates which terms appear in the system via the ``operators`` argument.
    %
    % Example
    % -------
    % .. code-block:: matlab
    %
    %  % Define dimensions, final time, and initial conditions.
    %  n_y = 20;            % ODE state dimension (at a fixed time).
    %  n_q = 4;             % ODE control dimension (at a fixed time).
    %  T = 1;               % Final simulation time.
    %  n_t = 400;           % Number of time steps.
    %  y0 = randn(n_y, 1);  % Initial condition.
    %
    %  % Define a list of operator terms and construct the constraint.
    %  operators = {Constant_Operator(), Linear_Operator(), ...
    %               Quadratic_Operator(), Input_Operator()};
    %  constraint = OpInf_ROM_Constraint(n_y, 4, T, 400, y0, operators);

    properties
        regularizer         % Tikhonov regularization matrix.
        operators           % Operator objects defining the differential system.
        y0                  % Initial condition :math:`\y(0)\in\R^{n_y}` for the ODE.
    end

    properties (SetAccess = protected)
        n_q                 % Dimension :math:`n_q` of the input :math:`\q(t)`.
        ds
        d                   % Number of operator entries for each system mode.
        inferred_operators  % Track which operators need to be learned.
    end

    methods

        %% Property setters.

        function set.operators(this, operators)
            % Set the dynamical system operators. These may be uninitialized,
            % meaning Set_Operators() has not been called.
            num_operators = length(operators);
            this.ds = zeros(num_operators, 1);
            this.inferred_operators = zeros(num_operators, 1);
            for i = 1:num_operators
                op = operators{i};
                this.ds(i) = op.Column_Dimension(this.n_y, this.n_q);
                if size(op.entries, 1) == 0
                    this.inferred_operators(i) = 1;
                end
            end
            this.d = sum(this.ds);
            this.operators = operators;
        end

        function set.regularizer(this, reg)
            % Set the Tikhonov regularizer.
            %
            % Parameters
            % ----------
            % reg
            %   Tikhonov regularization, one of the following:
            %
            %   * scalar: regularizer = reg * I
            %   * vector of length ``length(this.operators)``:
            %     regularizer = diag(reg(1)I_1, reg(2)I_2, ...) where each I_i
            %     is the column dimension of the ith operator.
            %   * vector of length d: regularizer = diag(reg)
            %   * matrix of size (d, d): regularizer = reg.

            % Scalar regularizer.
            if isscalar(reg)
                % Scalar regularizer --> Gamma = regularizer * I
                reg = reg * eye(this.d);
            end

            % Vector regularizer.
            num_operators = length(this.operators);
            if min(size(reg)) == 1
                % Vector regularizer
                if length(reg) == num_operators
                    % Scalar regularizer for each operator
                    newreg = [];
                    for i = 1:num_operators
                        newreg = horzcat(newreg, reg(i) * ones(1, this.ds(i)));
                    end
                    reg = newreg;
                end
                reg = diag(reg);
            end

            % Matrix regularizer.
            if size(reg, 1) ~= this.d || size(reg, 2) ~= this.d
                error('invalid regularizer');
            end

            this.regularizer = reg;
        end

        function set.y0(this, init)
            % Set the initial condition.
            if length(init) ~= this.n_y
                error('y0 must have n_y entries');
            end
            this.y0 = reshape(init, this.n_y, 1);
        end

    end

    methods (Access = public)

        %% Constructor.

        function this = OpInf_ROM_Constraint(n_y, n_q, T, n_t, y0, operators)
            % Parameters
            % ----------
            % n_y
            %   Dimension :math:`n_y` of the state :math:`\y_{j}` at each time.
            %   Note that if this is a reduced-order model, this is the REDUCED
            %   state dimension, denoted elsewhere as :math:`n_y'`.
            % n_q
            %   Dimension :math:`n_q` of the input :math:`\q`.
            % T
            %   Final time :math:`T > 0`.
            % n_t
            %   Number of nodes :math:`N` in the time mesh.
            % y0
            %   Initial condition :math:`\y(0)\in\R^{n_y}` for the ODE.
            % operators
            %   ODE system operators. These may be uninitialized,
            %   i.e., with null entries.

            this@Dynamic_Constraint(n_y, n_q * (n_t - 1), T, n_t);
            % Note: no control at the initial time step b/c bkwd Euler.
            this.n_q = n_q;
            this.y0 = y0;
            this.operators = operators;
            this.regularizer = 0;
        end

        %% Training.

        function [Y] = State_Solve2(this, Q)
            % Solve the constraint equation :math:`\c(\u,\z) = \0` by
            % integrating the ordinary differential equation
            % :math:`\ddt\y(t) = \f(\y(t),\z,t),~~\y(0)=\y_0`
            % in time using the first-order implicit Euler method.
            %
            % This method calls :meth:`State_Solve` with some convenience reshaping.
            %
            % Parameters
            % ----------
            % Q
            %   Control profile :math:`\Q\in\R^{n_q\times(n_t - 1)}` representing
            %   the :math:`n_q` control nodes at all time points except for the
            %   initial time (because of the Backward Euler scheme).
            %
            % Returns
            % -------
            % Y : :math:`n_y \times n_t` matrix
            %   State solution :math:`[~\y_1~~\cdots~~\y_{n_t}~]\in\R^{n_y\times n_t}`
            %   where :math:`\y_j` is the ODE state at time :math:`t_j`.
            Z = reshape(Q, this.n_z, 1);
            U = this.State_Solve(Z);
            Y = reshape(U, this.n_y, this.n_t);
        end

        function [dYdt] = Estimate_State_ddts(this, Y)
            % Use first-order backward differences to estimate the time
            % derivatives of the training states.
            %
            % Parameters
            % ----------
            % Y
            %   State snapshots :math:`\Y \in \R^{n_y \times n_t}`.
            %   This includes the initial condition.
            %
            % Returns
            % -------
            % dYdt
            %   Time derivative corresponding to the snapshots
            %   except at the initial condition,
            %   :math:`\dot{\Y} \in \R^{n_y \times (n_t - 1)}`,

            dt = this.t_mesh(2) - this.t_mesh(1);
            dYdt = (Y(:, 2:end) - Y(:, 1:(end - 1))) / dt;
        end

        function [dYdt] = Estimate_State_ddts_2ndOrder(this, Y)
            % Use second-order differences to estimate the time
            % derivatives of the training states.
            %
            % Parameters
            % ----------
            % Y
            %   State snapshots :math:`\Y \in \R^{n_y \times n_t}`.
            %   This includes the initial condition.
            %
            % Returns
            % -------
            % dYdt
            %   Time derivative corresponding to the snapshots
            %   except at the initial condition,
            %   :math:`\dot{\Y} \in \R^{n_y \times (n_t - 1)}`,

            dt = this.t_mesh(2) - this.t_mesh(1);
            dYdt = zeros(size(Y));
            % Forward difference for the first point.
            % dYdt(:, 1) = -3 * Y(:, 1) + 4 * Y(:, 2) - Y(:, 3);
            % Central difference for the interior points.
            dYdt(:, 2:end - 1) = -Y(:, 1:end - 2) + Y(:, 3:end);
            % Backward difference for the last point.
            dYdt(:, end) = 3 * Y(:, end - 1) - 4 * Y(:, end - 2) + Y(:, end - 3);
            dYdt = dYdt(:, 2:end) / (2 * dt);
        end

        function [dYdt] = Estimate_State_ddts_4thOrder(this, Y)
            % Use fourth-order differences to estimate the time
            % derivatives of the training states.
            %
            % Parameters
            % ----------
            % Y
            %   State snapshots :math:`\Y \in \R^{n_y \times n_t}`.
            %   This includes the initial condition.
            %
            % Returns
            % -------
            % dYdt
            %   Time derivative corresponding to the snapshots
            %   except at the initial condition,
            %   :math:`\dot{\Y} \in \R^{n_y \times (n_t - 1)}`,
            dt = this.t_mesh(2) - this.t_mesh(1);
            dYdt = zeros(size(Y));
            % Forward difference for the first two points.
            for j = 2:2
                dYdt(:, j) = -25 * Y(:, j) + 48 * Y(:, j + 1) - 36 * Y(:, j + 2) + 16 * Y(:, j + 3) - 3 * Y(:, j + 4);
            end
            % Central difference for the interior points.
            dYdt(:, 3:end - 2) = Y(:, 1:end - 4) - 8 * Y(:, 2:end - 3) + 8 * Y(:, 4:end - 1) - Y(:, 5:end);
            % Backward difference for the last three points.
            for j = 0:1
                dYdt(:, end - j) = 25 * Y(:, end - j) - 48 * Y(:, end - j - 1) + 36 * Y(:, end - j - 2) - 16 * Y(:, end - j - 3) + 3 * Y(:, end - j - 4);
            end
            dYdt = dYdt(:, 2:end) / (12 * dt);
            % Yk = Y(:, 3:end - 2);
            % Qk = Q(:, 2:end - 2);
            % dYdt = (Y(:, 1:end - 4) - 8 * Y(:, 2:end - 3) + 8 * Y(:, 4:end - 1) - Y(:, 5:end)) / (12 * dt);
        end

        function [dYdt] = Estimate_State_ddts_6thOrder(this, Y)
            % Use sixth-order differences to estimate the time
            % derivatives of the training states.
            %
            % Parameters
            % ----------
            % Y
            %   State snapshots :math:`\Y \in \R^{n_y \times n_t}`.
            %   This includes the initial condition.
            %
            % Returns
            % -------
            % dYdt
            %   Time derivative corresponding to the snapshots
            %   except at the initial condition,
            %   :math:`\dot{\Y} \in \R^{n_y \times (n_t - 1)}`,
            dt = this.t_mesh(2) - this.t_mesh(1);
            dYdt = zeros(size(Y));
            % Forward difference for the first three points.
            for j = 2:3
                dYdt(:, j) = -147 * Y(:, j) + 360 * Y(:, j + 1) - 450 * Y(:, j + 2) + 400 * Y(:, j + 3) - 225 * Y(:, j + 4) + 72 * Y(:, j + 5) - 10 * Y(:, j + 6);
            end
            % Central difference for the interior points.
            dYdt(:, 4:end - 3) = -Y(:, 1:end - 6) + 9 * Y(:, 2:end - 5) - 45 * Y(:, 3:end - 4) + 45 * Y(:, 5:end - 2) - 9 * Y(:, 6:end - 1) + Y(:, 7:end);
            % Backward difference for the last three points.
            for j = 0:2
                dYdt(:, end - j) = 147 * Y(:, end - j) - 360 * Y(:, end - j - 1) + 450 * Y(:, end - j - 2) - 400 * Y(:, end - j - 3) + 225 * Y(:, end - j - 4) - 72 * Y(:, end - j - 5) + 10 * Y(:, end - j - 6);
            end
            dYdt = dYdt(:, 2:end) / (60 * dt);
            % Yk = Y(:, 4:end - 3);
            % Qk = Q(:, 3:end - 3);
            % dYdt = (-Y(:, 1:end - 6) + 9 * Y(:, 2:end - 5) - 45 * Y(:, 3:end - 4) + 45 * Y(:, 5:end - 2) - 9 * Y(:, 6:end - 1) + Y(:, 7:end)) / (60 * dt);
        end

        function Learn_Operators(this, Y, Q, dYdt)
            % Infer the entries of the system operators by solving the
            % regression problem
            %
            % .. math::
            %  \min_{\c,\A,\H,\B}\sum_{j=1}^{n_t}\left\|
            %  \c + \A\y_j + \H[\y_j\otimes\y_j] + \B\q_j - \dot{\y}_j
            %  \right\|_{2}^{2}
            %  + \left\|\bfGamma[~\c~~\A~~\H~~\B~]\trp\right\|_{F}^{2}
            %
            % where
            %
            % * :math:`\Y = [~\y_1~~\cdots~~\y_{n_t}~] \in \R^{n_y \times n_t}`
            %   are state snapshots,
            % * :math:`\dot{\Y} = [~\dot{\y}_1~~\cdots~~\dot{\y}_{n_t}~]
            %   \in \R^{n_y \times n_t}` collects the time derivative of the
            %   state snapshots,
            % * :math:`\Q = [~\q_1~~\cdots~~\q_{n_t}~] \in \R^{n_q \times n_t}`
            %   are the corresponding inputs, and
            % * :math:`\bfGamma\in\R^{d \times d}` is a Tikhonov regularization
            %   matrix (the ``regularizer`` property).
            %
            % Parameters
            % ----------
            % Y
            %   State snapshots :math:`\Y \in \R^{n_y \times n_t}`.
            % Q
            %   Inputs :math:`\Q \in \R^{n_q \times n_t}`.
            % dYdt
            %   Time derivative of the state snapshots,
            %   :math:`\dot{\Y} \in \R^{n_y \times n_t}`.
            arguments
                this
                Y (:, :) {mustBeNumeric}
                Q (:, :) {mustBeNumeric}
                dYdt (:, :) {mustBeNumeric}
            end

            % Check that there is something to do.
            if sum(this.inferred_operators) == 0
                error('All operator entries already populated');
            end

            % Check that Y, Q, and dYdt have appropriate sizes.
            if ~isequal(size(Y, 2), size(Q, 2))
                error('Y and Q not aligned');
            elseif ~isequal(size(Y), size(dYdt))
                error('Y and dYdt not aligned');
            end

            % Construct the data matrix.
            rhs = dYdt;
            num_operators = length(this.operators);
            D = [];
            for i = 1:num_operators
                op = this.operators{i};
                if this.inferred_operators(i) == 1
                    D = horzcat(D, op.Datablock(Y, Q)');
                else
                    rhs = rhs - op.Apply(Y, Q);
                end
            end

            % Solve the regression problem for the operator entries.
            if norm(this.regularizer) > 0
                % Tikhonov regularization.
                Ohat = lscov([D; this.regularizer], [rhs'; zeros(this.d, this.n_y)])';
            else
                % No regularization.
                Ohat = lscov(D, rhs')';
            end

            % Unpack the operator entries.
            index = 1;
            for i = 1:num_operators
                if this.inferred_operators(i) == 1
                    newindex = index + this.ds(i);
                    this.operators{i}.Set_Entries(Ohat(:, index:(newindex - 1)));
                    index = newindex;
                end
            end
        end

        function [best_reg] = Select_Regularization(this, states, controls, reg_candidates, ddt_strategy)
            % Use a grid search to select a scalar regularization hyperparameter.
            %
            % Parameters
            % ----------
            % states
            %   Compressed training states :math:`\Y_1,\ldots\Y_k\in\R^{n_y \times n_t}`
            %   where :math:`k` is the number of trajectories. This is an
            %   :math:`n_y \times n_t \times k` tensor or (optionally), if :math:`k = 1`,
            %   an :math:`n_y \times n_t` matrix.
            % controls
            %   Control profiles :math:`\Q_1,\ldots,\Q_k\in\R^{n_q \times (n_t - 1)}`
            %   corresponding to the training states, where :math:`k` is the number of
            %   trajectories. This is an :math:`n_q \times (n_t - 1) \times k` tensor or
            %   (optionally), if :math:`k = 1`, an :math:`n_q \times (n_t - 1)` matrix.
            % reg_candidates
            %   Candidate regularization values to check.
            %
            % Returns
            % -------
            % best_reg : float
            %   Best regularization hyperparameter.
            arguments
                this
                states {mustBeNumeric}
                controls {mustBeNumeric}
                reg_candidates {mustBeNumeric}
                ddt_strategy {mustBeText} = "bwd1"
            end

            if length(size(states)) == 2
                states = reshape(states, size(states, 1), size(states, 2), 1);
            end
            if length(size(controls)) == 2
                controls = reshape(controls, size(controls, 1), size(controls, 2), 1);
            end
            num_trajectories = size(states, 3);
            num_snapshots_per_trajectory = this.n_t - 1;

            if size(states, 1) ~= this.n_y
                error('states should have n_y rows');
            elseif size(states, 2) ~= this.n_t
                error('states should have n_t columns');
            elseif size(controls, 1) ~= this.n_q
                error('controls should have n_q rows');
            elseif size(controls, 2) ~= num_snapshots_per_trajectory
                error('controls should have n_t - 1 columns');
            elseif size(controls, 3) ~= num_trajectories
                error('states and controls not aligned, different number of trajectories');
            end

            % Estimate time derivatives.
            if ddt_strategy == "bwd1"
                driver = @this.Estimate_State_ddts;
            elseif ddt_strategy == "2ndOrder"
                driver = @this.Estimate_State_ddts_2ndOrder;
            elseif ddt_strategy == "4thOrder"
                driver = @this.Estimate_State_ddts_4thOrder;
            elseif ddt_strategy == "6thOrder"
                driver = @this.Estimate_State_ddts_6thOrder;
            else
                error('unexpected ddt_strategy');
            end

            Y = zeros(this.n_y, num_snapshots_per_trajectory, num_trajectories);
            dYdt = zeros(this.n_y, num_snapshots_per_trajectory, num_trajectories);
            for k = 1:num_trajectories
                dYk = driver(states(:, :, k));
                Y(:, :, k) = states(:, 2:end, k);
                dYdt(:, :, k) = dYk;
            end

            num_training_snapshots = num_snapshots_per_trajectory * num_trajectories;
            Y = reshape(Y, this.n_y, num_training_snapshots);
            dYdt = reshape(dYdt, this.n_y, num_training_snapshots);
            Q = reshape(controls, this.n_q, num_training_snapshots);

            num_candidates = size(reg_candidates, 2);
            reconstruction_errors = zeros(1, num_candidates);
            original_initial_condition = this.y0;
            disp(['Regularization selection (' num2str(num_candidates), ' candidates)']);
            for i = 1:num_candidates
                reg = reg_candidates(i);
                this.regularizer = reg;
                this.Learn_Operators(Y, Q, dYdt);

                % Solve the model for each training control profile.
                total_error = 0;
                for k = 1:num_trajectories
                    this.y0 = states(:, 1, k);
                    Yk = this.State_Solve2(controls(:, :, k));
                    local_error = norm(Yk - states(:, :, k)) / norm(states(:, :, k));
                    total_error = total_error + local_error;
                end
                reconstruction_errors(i) = total_error;

                disp(['reg = ', num2str(reg), '; error = ', num2str(total_error)]);
            end
            this.y0 = original_initial_condition;

            % Choose the best regularization.
            [err, idx] = min(reconstruction_errors);
            best_reg = reg_candidates(idx);
            disp(['winner = ', num2str(best_reg), '; error = ', num2str(err)]);

            this.regularizer = best_reg;
            this.Learn_Operators(Y, Q, dYdt);
        end

        %% Implement abstract methods from the parent class.

        function [f, f_y, f_z] = f(this, y, z, t)

            % Extract the input q(t) from the control z based on the time t.
            I = this.Input_Indices(t);
            q = z(I);

            % Allocate space for outputs.
            f = zeros(this.n_y, 1);
            f_y = zeros(this.n_y, this.n_y);
            f_z = zeros(this.n_y, this.n_z);

            % Evaluate each operator.
            for i = 1:length(this.operators)
                op = this.operators{i};
                f = f + op.Apply(y, q);
                f_y = f_y + op.Jacobian_y(y, q);
                f_z(:, I) = f_z(:, I) + op.Jacobian_q(y, q);
            end
        end

        function [h, h_z] = h(this, ~)
            h = this.y0;
            h_z = zeros(this.n_y, this.n_z);
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);
            Mv = zeros(this.n_y, size(v, 2));
            for i = 1:length(this.operators)
                Mv = Mv + this.operators{i}.Hessian_yy_Apply(v, y, q, lambda);
            end
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);
            vt = v(I, :);
            Mv = zeros(this.n_y, size(v, 2));
            for i = 1:length(this.operators)
                Mv = Mv + this.operators{i}.Hessian_yq_Apply(vt, y, q, lambda);
            end
        end

        function [Mv] = f_zy_Apply(this, v, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);
            Mv = zeros(this.n_z, size(v, 2));
            for i = 1:length(this.operators)
                Mv(I, :) = Mv(I, :) + this.operators{i}.Hessian_qy_Apply(v, y, q, lambda);
            end
        end

        function [Mv] = f_zz_Apply(this, v, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);
            vt = v(I, :);
            Mv = zeros(this.n_z, size(v, 2));
            for i = 1:length(this.operators)
                Mv(I, :) = Mv(I, :) + this.operators{i}.Hessian_yq_Apply(vt, y, q, lambda);
            end
        end

        function [Mv] = h_zz_Apply(this, v, ~, ~)
            Mv = zeros(this.n_z, size(v, 2));
        end

    end

    methods (Access = protected)

        function [mask] = Input_Indices(this, t)
            [~, t_index] = min(abs(t - this.t_mesh));
            if t_index == 1
                error('no control at initial time!');
            end
            idx = t_index - 1;
            mask = (this.n_q * (idx - 1) + 1):(this.n_q * idx);
        end

    end
end
