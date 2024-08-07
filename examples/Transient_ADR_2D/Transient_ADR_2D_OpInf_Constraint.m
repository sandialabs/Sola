classdef Transient_ADR_2D_OpInf_Constraint < Dynamic_Constraint
    % Operator Inference Reduced-order Model for a two-species
    % advection-diffusion-reaction system.
    %
    % This class represents a system of ordinary differential equations for
    % a state variable :math:`\y(t)\in\R^{n_y}` which can be partitioned into
    % two substates:
    %
    % .. math::
    %  \y(t) = \left[\begin{array}{c}
    %  \y_1(t) \\ \y_2(t)
    %  \end{array}\right]
    %
    % where each :math:`\y_\ell\in\R^{n_\ell}`, with
    % :math:`n_1 + n_2 = n_y`. The ODE system is given by
    %
    % .. math::
    %  \ddt\y_1(t) = \A_1\y_1(t) + \H_1[\y_1(t)\otimes\y_2(t)],
    %  \\
    %  \ddt\y_2(t) = \A_2\y_2(t) + \H_2[\y_1(t)\otimes\y_2(t)] + \B\q(t),
    %
    % where :math:`\A_\ell\in\R^{n_\ell\times n_\ell}`,
    % :math:`\H_\ell\in\R^{n_\ell\times n_1 n_2}`,
    % and :math:`\B\in\R^{n_2\times n_q}`.
    %
    % This class inherits from :class:`Dynamic_Constraint`, hence the model is integrated in time
    % using the implicit Euler method. In this context, the optimization state and control are
    % given by
    %
    % .. math::
    %  \u = \left(\begin{array}{c}\y(t_1) \\ \vdots \\ \y(t_{n_t})\end{array}\right)\in\R^{n_u},
    %  \qquad
    %  \z = \left(\begin{array}{c}\q(t_2) \\ \vdots \\ \q(t_{n_t})\end{array}\right)\in\R^{n_z},
    %
    % with :math:`n_u = n_y n_t` and :math:`n_z = n_q (n_t - 1)`.
    % Note that the state and the control share a time mesh, which is important for the definition
    % of the regression problem. However, the control is not measured at the initial condition
    % because the system is solved with a one-step implicit time integration scheme.

    properties
        toinfer             % Which operators to infer the entries of.
        y0                  % Initial condition :math:`\y(0)\in\R^{n_y}` for the ODE.
    end

    properties (SetAccess = protected)
        n_q                 % Dimension :math:`n_q` of the input :math:`\q(t)`.
        n_1                 % Dimension of the first substate, :math:`n_1`.
        n_2                 % Dimension of the second substate, :math:`n_2`.
        d_1                 % First operator dimension, :math:`d_1 = n_1 + n_1 n_2`.
        d_2                 % Second operator dimension, :math:`d_2 = n_2 + n_1 n_2 + n_q`.
        A_1                 % Linear term for the first variable, :math:`\A_1\y_1(t)`.
        H_1                 % Quadratic term for the first variable, :math:`\H_1\[y_1(t)\otimes\y_2(t)]`.
        A_2                 % Linear term for the second variable, :math:`\A_2\y_2(t)`.
        H_2                 % Quadratic term for the second variable, :math:`\H_2\[y_1(t)\otimes\y_2(t)]`.
        B_2                 % Input term for the second variable, :math:`\B\q(t)`.
    end

    methods

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

        function this = Transient_ADR_2D_OpInf_Constraint(n_1, n_2, n_q, T, n_t, y0)
            % Parameters
            % ----------
            % n_1
            %   Dimension :math:`n_1` of the first substate, :math:`\y_1(t)`.
            % n_2
            %   Dimension :math:`n_2` of the second substate, :math:`\y_2(t)`.
            %   The full state dimension is :math:`n_y = n_1 + n_2`.
            % n_q
            %   Dimension :math:`n_q` of the input :math:`\q(t)`.
            % T
            %   Final time :math:`T > 0`.
            % n_t
            %   Number of nodes :math:`n_t` in the time mesh.
            % y0
            %   Initial condition :math:`\y(0)\in\R^{n_y}` for the ODE.

            n_y = n_1 + n_2;
            this@Dynamic_Constraint(n_y, n_q * (n_t - 1), T, n_t);
            % Note: no control at the initial time step b/c bkwd Euler.

            % Dimensions.
            this.n_1 = n_1;
            this.n_2 = n_2;
            n1n2 = n_1 * n_2;
            this.d_1 = n_1 + n1n2;
            this.d_2 = n_2 + n1n2 + n_q;
            this.n_q = n_q;

            % Operators.
            dims = [n_1; n_2];
            this.A_1 = Linear_Operator_Multi(1, 1, dims);
            this.H_1 = Quadratic_Operator_Multi(1, 1, 2, dims);
            this.A_2 = Linear_Operator_Multi(2, 2, dims);
            this.H_2 = Quadratic_Operator_Multi(2, 1, 2, dims);
            this.B_2 = Input_Operator_Multi(2, n_q, dims);

            % Operator Inference hyperparameters.
            this.toinfer = [1; 1; 1; 1; 1];

            % Initial condition.
            this.y0 = y0;
        end

        %% Training.

        function [Y] = State_Solve2(this, Q)
            % Solve the constraint equation :math:`\c(\u,\z) = \0` by
            % integrating the model in time using the first-order implicit
            % Euler method.
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
            %   State solution
            %   :math:`[~\y(t_1)~~\cdots~~\y(t_{n_t})~]\in\R^{n_y\times n_t}`.
            z = reshape(Q, this.n_z, 1);
            u = this.State_Solve(z);
            Y = reshape(u, this.n_y, this.n_t);
        end

        function [dYdt] = Estimate_State_ddts_1stOrder(this, Y)
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

        function [Reg1, Reg2] = Assemble_Regularizer(this, ABreg, Hreg)
            % Assemble Tikhonov regularization matrices for each regression.
            %
            % Parameters
            % ----------
            % ABreg
            %   Regularization scalar for linear and input terms.
            % Hreg
            %   Regularization scalar for linear and input terms.

            arguments
                this
                ABreg {mustBeNumeric}
                Hreg {mustBeNumeric}
            end

            % [A_1 H_1]
            Reg1 = zeros(this.d_1, 1);
            Reg1(1:this.n_1) = ABreg;
            Reg1(this.n_1 + 1:end) = Hreg;
            Reg1 = diag(Reg1);

            % [A_2 H_2 B_2]
            Reg2 = zeros(this.d_2, 1);
            Reg2(1:this.n_2) = ABreg;
            Reg2(this.n_2 + 1:this.n_2 + (this.n_1 * this.n_2)) = Hreg;
            Reg2(this.n_2 + (this.n_1 * this.n_2) + 1:end) = ABreg;
            Reg2 = diag(Reg2);
        end

        function Learn_Operators(this, Y, Q, dYdt, abreg, hreg)
            % Infer the entries of the system operators by solving the
            % regression problems
            %
            % .. math::
            %  \min_{\A_1,\H_1}\sum_{j=1}^{n_t}\left\|
            %  \A_1\y_1(t_j) + \H_1[\y_1(t_j)\otimes\y_2(t_j)] - \dot{\y}_1(t_j)
            %  \right\|_{2}^{2}
            %  + \gamma_1^2\|\A_1\|_{F}^{2}
            %  + \gamma_2^2\|\H_1\|_{F}^{2},
            %  \\
            %  \min_{\A_2,\H_2,\B_2}\sum_{j=1}^{n_t}\left\|
            %  \A_2\y_2(t_j) + \H_2[\y_1(t_j)\otimes\y_2(t_j)] + \B_2\q(t_j) - \dot{\y}_2(t_j)
            %  \right\|_{2}^{2}
            %  + \gamma_1^2\|\A_2\|_{F}^{2}
            %  + \gamma_2^2\|\H_2\|_{F}^{2}
            %  + \gamma_1^2\|\B_2\|_{F}^{2}.
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
                abreg = 0
                hreg = 0
            end

            % Check that Y, Q, and dYdt have appropriate sizes.
            k = size(Y, 2);
            if ~isequal(size(Q, 2), k)
                error('Y and Q not aligned');
            elseif ~isequal(size(Y), size(dYdt))
                error('Y and dYdt not aligned');
            end

            % Construct the data matrices.
            rhs1 = dYdt(1:this.n_1, :);
            D1 = [];

            if this.toinfer(1) == 1
                D1 = vertcat(D1, this.A_1.Datablock(Y, Q));
            else
                rhs1 = rhs1 - this.A_1.Apply(Y, Q);
            end

            if this.toinfer(2) == 1
                D1 = vertcat(D1, this.H_1.Datablock(Y, Q));
            else
                rhs1 = rhs1 - this.H_1.Apply(Y, Q);
            end

            D2 = [];
            rhs2 = dYdt(this.n_1 + 1:end, :);

            if this.toinfer(3) == 1
                D2 = vertcat(D2, this.A_2.Datablock(Y, Q));
            else
                rhs2 = rhs2 - this.A_2.Apply(Y, Q);
            end

            if this.toinfer(4) == 1
                D2 = vertcat(D2, this.H_2.Datablock(Y, Q));
            else
                rhs2 = rhs2 - this.H_2.Apply(Y, Q);
            end

            if this.toinfer(5) == 1
                D2 = vertcat(D2, this.B_2.Datablock(Y, Q));
            else
                rhs2 = rhs2 - this.B_2.Apply(Y, Q);
            end

            D1 = D1';
            D2 = D2';

            [Reg1, Reg2] = this.Assemble_Regularizer(abreg, hreg);

            % Solve the regression problem for the operator entries.
            if norm(Reg1) > 0
                % Tikhonov regularization.
                Ohat1 = lscov([D1; Reg1], [rhs1'; zeros(this.d_1, this.n_1)])';
            else
                % No regularization.
                Ohat1 = lscov(D1, rhs1')';
            end

            if norm(Reg2) > 0
                % Tikhonov regularization.
                Ohat2 = lscov([D2; Reg2], [rhs2'; zeros(this.d_2, this.n_2)])';
            else
                % No regularization.
                Ohat2 = lscov(D2, rhs2')';
            end

            % Unpack the operator entries.
            index = 1;
            if this.toinfer(1) == 1
                newindex = index + this.n_1;
                this.A_1.Set_Entries(Ohat1(:, index:(newindex - 1)));
                index = newindex;
            end

            if this.toinfer(2) == 1
                newindex = index + (this.n_1 * this.n_2);
                this.H_1.Set_Entries(Ohat1(:, index:(newindex - 1)));
            end

            index = 1;
            if this.toinfer(3) == 1
                newindex = index + this.n_2;
                this.A_2.Set_Entries(Ohat2(:, index:(newindex - 1)));
                index = newindex;
            end

            if this.toinfer(4) == 1
                newindex = index + (this.n_1 * this.n_2);
                this.H_2.Set_Entries(Ohat2(:, index:(newindex - 1)));
                index = newindex;
            end

            if this.toinfer(5) == 1
                newindex = index + this.n_q;
                this.B_2.Set_Entries(Ohat2(:, index:(newindex - 1)));
            end
        end

        function [best_reg] = Select_Regularization(this, states, controls, ABreg_candidates, Hreg_candidates, ddt_strategy)
            % Use a grid search to select regularization hyperparameters.
            %
            % Parameters
            % ----------
            % states
            %   Compressed training states. This is either a cell of :math:`k`
            %   matrices of sizes :math:`\times n_y \times n_t` or (optionally),
            %   if :math:`k = 1`, a single :math:`n_y \times n_t` matrix.
            % controls
            %   Training control profiles corresponding to the training states.
            %   This is either a cell of :math:`k` matrices of sizes
            %   :math:`n_q \times (n_t - 1)` or (optionally), if :math:`k = 1`,
            %   a single :math:`n_q \times (n_t - 1)` matrix.
            % ABreg_candidates
            %   Candidate regularization values for linear model terms.
            % Hreg_candidates
            %   Candidate regularization values for quadratic model terms.
            % ddt_strategy
            %   Which finite difference strategy to use when estimating the
            %   time derivatives of the training states, one of
            %   ``"bwd1"`` (default), ``"2ndOrder"``, ``"4thOrder"``, or
            %   ``"6thOrder"``.
            %
            % Returns
            % -------
            % best_regs : [float; float]
            %   Best regularization hyperparameters [ABreg; Hreg].
            arguments
                this
                states {Transient_ADR_2D_OpInf_Constraint.mustBeNumericOrCell}
                controls {Transient_ADR_2D_OpInf_Constraint.mustBeNumericOrCell}
                ABreg_candidates {mustBeNumeric}
                Hreg_candidates {mustBeNumeric}
                ddt_strategy {mustBeText} = "bwd1"
            end

            if ~iscell(states)
                states = {states};
            end
            if ~iscell(controls)
                controls = {controls};
            end

            % Check states and controls are aligned.
            num_trajectories = length(states);
            num_snapshots_per_trajectory = this.n_t - 1;
            if length(controls) ~= num_trajectories
                error('states and controls not aligned, different number of trajectories');
            end

            % Select time derivative estimation strategy.
            if ddt_strategy == "bwd1"
                ddt_estimator = @this.Estimate_State_ddts_1stOrder;
            elseif ddt_strategy == "2ndOrder"
                ddt_estimator = @this.Estimate_State_ddts_2ndOrder;
            elseif ddt_strategy == "4thOrder"
                ddt_estimator = @this.Estimate_State_ddts_4thOrder;
            elseif ddt_strategy == "6thOrder"
                ddt_estimator = @this.Estimate_State_ddts_6thOrder;
            else
                error('invalid ddt_strategy');
            end

            % Check shapes and estimate time derivatives.
            Y = cell(num_trajectories);
            dYdt = cell(num_trajectories);
            for k = 1:num_trajectories
                Yk = states{k};
                Qk = controls{k};

                % Check shapes.
                if size(Yk, 1) ~= this.n_y
                    error('each state matrix should have n_y rows');
                elseif size(Yk, 2) ~= this.n_t
                    error('each state matrix should have n_t columns');
                elseif size(Qk, 1) ~= this.n_q
                    error('each control matrix should have n_q rows');
                elseif size(Qk, 2) ~= num_snapshots_per_trajectory
                    error('each control matrix should have n_t - 1 columns');
                end

                % Estimate time derivatives and strip off initial state.
                dYdt{k} = ddt_estimator(Yk);
                Y{k} = Yk(:, 2:end);
            end

            % Concatenate data from all trajectories.
            states_all = horzcat(Y{:});
            ddts_all = horzcat(dYdt{:});
            controls_all = horzcat(controls{:});
            original_initial_condition = this.y0;

            % Search the regularization space for the best performer.
            num_ABcandidates = numel(ABreg_candidates);
            num_Hcandidates = numel(Hreg_candidates);
            reconstruction_errors = zeros(num_ABcandidates, num_Hcandidates);
            num_experiments = num_ABcandidates * num_Hcandidates;
            disp(['Regularization selection (' num2str(num_experiments), ' candidates)']);
            for j = 1:num_Hcandidates
                Hreg = Hreg_candidates(j);
                for i = 1:num_ABcandidates
                    ABreg = ABreg_candidates(i);
                    this.Learn_Operators(states_all, controls_all, ddts_all, ABreg, Hreg);

                    total_error = 0;
                    for k = 1:num_trajectories
                        Yk_data = states{k};
                        this.y0 = Yk_data(:, 1);
                        Yk_rom = this.State_Solve2(controls{k});
                        local_error = norm(Yk_rom - Yk_data) / norm(Yk_data);
                        total_error = total_error + local_error;
                    end
                    reconstruction_errors(i, j) = total_error / num_trajectories;

                    fprintf('ABreg: %.2e, Hreg: %.2e; error = %.2e\n', ABreg, Hreg, total_error);
                end
            end
            this.y0 = original_initial_condition;

            % Choose the best regularization out of the candidates.
            [~, linearIndex] = min(reconstruction_errors(:));
            [rowIndex, colIndex] = ind2sub(size(reconstruction_errors), linearIndex);
            best_ABreg = ABreg_candidates(rowIndex);
            best_Hreg = Hreg_candidates(colIndex);
            best_err = reconstruction_errors(rowIndex, colIndex);
            fprintf('Best (ABreg, Hreg) = (%.2e, %.2e); error = %.2e\n', best_ABreg, best_Hreg, best_err);

            % Solve the problem again with the best hyperparameters.
            this.Learn_Operators(states_all, controls_all, ddts_all, best_ABreg, best_Hreg);
        end

        %% Implement abstract methods from the parent class.

        function [f, f_y, f_z] = f(this, y, z, t)

            % Extract the input q(t) from the control z based on the time t.
            I = this.Input_Indices(t);
            q = z(I);

            % Evaluate the operators for the first equation.
            f_1 = zeros(this.n_1, 1);
            f_y_1 = zeros(this.n_1, this.n_y);
            for op = {this.A_1, this.H_1}
                f_1 = f_1 + op{1}.Apply(y, q);
                f_y_1 = f_y_1 + op{1}.Jacobian_y(y, q);
            end

            % Evaluate the operators for the second equation.
            f_2 = zeros(this.n_2, 1);
            f_y_2 = zeros(this.n_2, this.n_y);
            for op = {this.A_2, this.H_2, this.B_2}
                f_2 = f_2 + op{1}.Apply(y, q);
                f_y_2 = f_y_2 + op{1}.Jacobian_y(y, q);
            end

            % Piece together the two results.
            f = [f_1; f_2];
            f_y = [f_y_1; f_y_2];
            f_z = zeros(this.n_y, this.n_z);
            f_z(this.n_1 + 1:end, I) = this.B_2.Jacobian_q(y, q);
        end

        function [Mv] = f_yy_Apply(this, v, ~, ~, ~, lambda)
            % I = this.Input_Indices(t);
            % q = z(I);
            Mv_1 = this.H_1.Hessian_yy_Apply(v, [], [], lambda);
            Mv_2 = this.H_2.Hessian_yy_Apply(v, [], [], lambda);
            Mv = Mv_1 + Mv_2;
        end

        function [Mv] = f_yz_Apply(this, v, ~, ~, ~, ~)
            Mv = zeros(this.n_y, size(v, 2));
        end

        function [Mv] = f_zy_Apply(this, v, ~, ~, ~, ~)
            Mv = zeros(this.n_z, size(v, 2));
        end

        function [Mv] = f_zz_Apply(this, v, ~, ~, ~, ~)
            Mv = zeros(this.n_z, size(v, 2));
        end

        function [h, h_z] = h(this, ~)
            h = this.y0;
            h_z = zeros(this.n_y, this.n_z);
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

    methods (Static, Access = private)

        function mustBeNumericOrCell(value)
            % Check if the value is numeric.
            if isnumeric(value)
                return
            end

            % If it's a cell, check each entry for being numeric.
            if iscell(value)
                if all(cellfun(@isnumeric, value))
                    return
                else
                    error('All cell entries must be numeric.');
                end
            end

            % If neither condition is met, throw an error
            error('Input must be either numeric or a cell array of numerics.');
        end

    end
end
