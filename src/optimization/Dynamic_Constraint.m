classdef Dynamic_Constraint < Constraint
    % Define constraint equations through an ordinary differential equation
    %
    % .. math:: \frac{\textrm{d}}{\textup{d}t}\y(t) &= \f(\y(t),\z,t),
    %  \qquad t \in [0, T], \\
    %  \y(0) &= \h(\z),
    %
    % where
    %
    % * :math:`T > 0` is the final time,
    % * :math:`\y(t) \in \R^{n_y}` is the state at time :math:`t`,
    % * :math:`\z \in \R^{n_z}` is the control,
    % * :math:`\f:\R^{n_y}\times\R^{n_z}\times\R\to\R^{n_y}` defines the ODE, and
    % * :math:`\h:\R^{n_z}\to\R^{n_y}` prescribes the initial condition.
    %
    % The ODE is integrated using the first-order implicit Euler method,
    %
    % .. math:: \frac{1}{\delta t}(\y_{j} - \y_{j-1}) = \f(\y_{j}, \z, t_{j}),
    %
    % where :math:`\y_{j}\in\R^{n_y}` approximates :math:`\y(t)`
    % at time :math:`t = t_j`. Hence, the constraints are given by
    %
    % .. math:: \c(\u,\z) = \left(\begin{array}{c}
    %  \y_1 - \h(\z) \\
    %  \y_2 - \y_1 - \delta t \f(\y_2, \z, t_2) \\
    %  \y_3 - \y_2 - \delta t \f(\y_3, \z, t_3) \\
    %  \vdots \\
    %  \y_{n_t} - \y_{N-1} - \delta t \f(\y_{n_t}, \z, t_{n_t})
    %  \end{array}\right),
    %  \qquad
    %  \u = \left(\begin{array}{c}
    %  \y_1 \\ \vdots \\ \y_{n_t}
    %  \end{array}\right).
    %
    % Instead of forming :math:`\c(\u,\z)` explicitly, MATLAB's ``fsolve()``
    % is used at each time step to solve for successive :math:`\y_{j}`.

    properties
        n_y                         % Dimension :math:`n_y` of the state :math:`\y_{j}` at each time.
        n_z                         % Dimension :math:`n_z` of the control :math:`\z`.
        n_t                         % Number of nodes :math:`n_t` in the time mesh.
        t_mesh                      % Time mesh :math:`(t_1,\ldots,t_{n_t})\trp`.
        time_step_solver_options    % Options for ``fsolve()``, used at each time step.
        verbose                     % Output verbosity.
    end

    properties (Dependent)
        T                           % Final time :math:`T = t_{n_t}`.
    end

    methods

        function finaltime = get.T(this)
            finaltime = this.t_mesh(end);
        end

    end

    methods (Abstract, Access = public)

        [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
        % *Abstract method.*
        % Evaluate the ODE function :math:`\f(\y,\z,t)` and its Jacobians
        % :math:`\f_{y}(\y,\z,t)` and :math:`\f_{z}(\y,\z,t)`.
        %
        % Parameters
        % ----------
        % y
        %   Differential equation state :math:`\y\in\R^{n_y}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % t
        %   Time :math:`t`.
        %
        % Returns
        % -------
        % f : vector
        %   Function value :math:`\f(\y,\z,t)\in\R^{n_y}`.
        % f_y : matrix
        %   Function Jacobian :math:`\f_{y}(\y,\z,t)\in\R^{n_y \times n_y}`.
        % f_z : matrix
        %   Function Jacobian :math:`\f_{z}(\y,\z,t)\in\R^{n_y \times n_z}`.

        [h, h_z] = Initial_Condition(this, z)
        % *Abstract method.*
        % Evaluate the ODE initial condition :math:`\h(\z)` and its Jacobian :math:`\h_{z}(\z)`.
        %
        % Parameters
        % ----------
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        %
        % Returns
        % -------
        % h : vector
        %   Function value :math:`\h(\z)\in\R^{n_y}`.
        % h_z : matrix
        %   Function Jacobian :math:`\h_{z}(\z)\in\R^{n_y \times n_z}`.

        [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)
        % *Abstract method.*
        % Compute the vector-Hessian-vector product :math:`\bflambda\trp \f_{y,y}(\y,\z,t)\v`.
        %
        % Parameters
        % ----------
        % v
        %   Search direction :math:`\v\in\R^{n_y}`.
        % y
        %   Differential equation state :math:`\y\in\R^{n_y}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % t
        %   Time :math:`t`.
        % lambda
        %   Adjoint (of :math:`\y`) :math:`\bflambda\in\R^{n_y}`.
        %
        % Returns
        % -------
        % Mv : vector
        %   Vector-Hessian-vector product
        % :math:`\bflambda\trp \f_{y,y}(\y,\z,t)\v\in\R^{n_y}`.

        [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
        % *Abstract method.*
        % Compute the vector-Hessian-vector product :math:`\bflambda\trp \f_{y,z}(\y,\z,t)\v`.
        %
        % Parameters
        % ----------
        % v
        %   Search direction :math:`\v\in\R^{n_z}`.
        % y
        %   Differential equation state :math:`\y\in\R^{n_y}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % t
        %   Time :math:`t`.
        % lambda
        %   Adjoint (of :math:`\y`) :math:`\bflambda\in\R^{n_y}`.
        %
        % Returns
        % -------
        % Mv : vector
        %   Vector-Hessian-vector product
        % :math:`\bflambda\trp \f_{y,z}(\y,\z,t)\v\in\R^{n_y}`.

        [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
        % *Abstract method.*
        % Compute the vector-Hessian-vector product :math:`\bflambda\trp \f_{z,y}(\y,\z,t)\v`.
        %
        % Parameters
        % ----------
        % v
        %   Search direction :math:`\v\in\R^{n_y}`.
        % y
        %   Differential equation state :math:`\y\in\R^{n_y}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % t
        %   Time :math:`t`.
        % lambda
        %   Adjoint (of :math:`\y`) :math:`\bflambda\in\R^{n_y}`.
        %
        % Returns
        % -------
        % Mv : vector
        %   Vector-Hessian-vector product
        % :math:`\bflambda\trp \f_{z,y}(\y,\z,t)\v\in\R^{n_z}`.

        [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
        % *Abstract method.*
        % Compute the vector-Hessian-vector product :math:`\bflambda\trp \f_{z,z}(\y,\z,t)\v`.
        %
        % Parameters
        % ----------
        % v
        %   Search direction :math:`\v\in\R^{n_z}`.
        % y
        %   Differential equation state :math:`\y\in\R^{n_y}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % t
        %   Time :math:`t`.
        % lambda
        %   Adjoint (of :math:`\y`) :math:`\bflambda\in\R^{n_y}`.
        %
        % Returns
        % -------
        % Mv : vector
        %   Vector-Hessian-vector product
        % :math:`\bflambda\trp \f_{z,z}(\y,\z,t)\v\in\R^{n_z}`.

        [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)
        % *Abstract method.*
        % Compute the vector-Jacobian-vector product :math:`\bflambda\trp \h_{z}(\z)\v`.
        %
        % Parameters
        % ----------
        % v
        %   Search direction :math:`\v\in\R^{n_z}`.
        % z
        %   Control :math:`\z\in\R^{n_z}`.
        % lambda
        %   Adjoint (of :math:`\y`) :math:`\bflambda\in\R^{n_y}`.
        %
        % Returns
        % -------
        % Mv : vector
        %   Vector-Jacobian-vector product
        % :math:`\bflambda\trp \h_{z}(\z)\v\in\R^{n_z}`

    end

    methods (Access = public)

        %% Implement abstract methods from the parent class

        function [u] = State_Solve(this, z)
            % Solve the constraint equation :math:`\c(\u,\z) = \0` by
            % integrating the ordinary differential equation
            % :math:`\frac{\textup{d}}{\textup{d}t}\y(t) = \f(\y(t),\z,t)`
            % in time using the first-order implicit Euler method.
            %
            % Parameters
            % ----------
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            %
            % Returns
            % -------
            % u : vector
            %   State :math:`\u = (\y_1\trp~\cdots~\y_{n_t}\trp)\trp\in\R^{n_u}`
            %   where :math:`\y_j` is the ODE state at time :math:`t_j`.

            u = zeros(this.n_y * this.n_t, 1);
            u(1:this.n_y) = this.Initial_Condition(z);
            for k = 2:this.n_t
                Im = ((k - 2) * this.n_y + 1):((k - 1) * this.n_y);      % y_{k-1} = u(Im)
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);           % y_{k} = u(I)
                u(I) = this.State_Eq_Time_Step(u(Im), z, this.t_mesh(k), this.t_mesh(k - 1));
            end
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y * this.n_t, num_vecs);
            I = ((this.n_t - 1) * this.n_y + 1):(this.n_t * this.n_y);      % y_{k} = u(I)
            dt = this.t_mesh(end) - this.t_mesh(end - 1);
            Mv(I, :) = this.Linearized_Adjoint_Time_Step_Solve(v(I, :), u(I), z, this.t_mesh(end), dt);

            for k = (this.n_t - 1):-1:2
                Im = (k * this.n_y + 1):((k + 1) * this.n_y);           % y_{k} = u(Im)
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);            % y_{k-1} = u(I)
                dt = this.t_mesh(k + 1) - this.t_mesh(k);
                Mv(I, :) = this.Linearized_Adjoint_Time_Step_Solve(v(I, :) + Mv(Im, :), u(I), z, this.t_mesh(k), dt);
            end
            I = 1:(this.n_y);
            Im = (this.n_y + 1):(2 * this.n_y);
            Mv(I, :) = Mv(Im, :) + v(I, :);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            [~, h_z] = this.Initial_Condition(z);
            Mv = -h_z' * v(1:this.n_y, :);
            for k = 2:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);            % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                [~, ~, f_z] = this.Time_Instance_RHS(u(I), z, this.t_mesh(k));
                Mv = Mv - dt * f_z' * v(I, :);
            end
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y * this.n_t, num_vecs);
            I = 1:this.n_y;
            Mv(I, :) = v(I, :);
            for k = 2:this.n_t
                Im = ((k - 2) * this.n_y + 1):((k - 1) * this.n_y);     % y_{k-1} = u(Im)
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);            % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                Mv(I, :) = this.Linearized_Time_Step_Solve(v(I, :) + Mv(Im, :), u(I), z, this.t_mesh(k), dt);
            end
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y * this.n_t, num_vecs);
            I = 1:this.n_y;
            [~, h_z] = this.Initial_Condition(z);
            Mv(I, :) = -h_z * v;
            for k = 2:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);            % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                [~, ~, f_z] = this.Time_Instance_RHS(u(I), z, this.t_mesh(k));
                Mv(I, :) = Mv(I, :) - dt * f_z * v;
            end
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y * this.n_t, num_vecs);
            for k = 2:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);            % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                f_yy = this.Time_Instance_RHS_yy_Apply(v(I, :), u(I), z, this.t_mesh(k), lambda(I));
                Mv(I, :) = -dt * f_yy;
            end
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y * this.n_t, num_vecs);
            for k = 2:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                f_yz = this.Time_Instance_RHS_yz_Apply(v, u(I), z, this.t_mesh(k), lambda(I));
                Mv(I, :) = -dt * f_yz;
            end
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_z, num_vecs);
            for k = 2:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                f_zy = this.Time_Instance_RHS_zy_Apply(v(I, :), u(I), z, this.t_mesh(k), lambda(I));
                Mv = Mv - dt * f_zy;
            end
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = -this.Initial_Condition_zz_Apply(v, z, lambda(1:this.n_y));
            for k = 2:this.n_t
                I = ((k - 1) * this.n_y + 1):(k * this.n_y);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                f_zz = this.Time_Instance_RHS_zz_Apply(v, u(I), z, this.t_mesh(k), lambda(I));
                Mv = Mv - dt * f_zz;
            end
        end

        %% Constructor

        function this = Dynamic_Constraint(n_y, n_z, T, n_t)
            % Parameters
            % ----------
            % n_y : int
            %   Dimension :math:`n_y` of the state :math:`\y_{j}` at each time.
            % n_z : int
            %   Dimension :math:`n_z` of the control :math:`\z`.
            % T : double
            %   Final time :math:`T`.
            % n_t : int
            %   Number of nodes :math:`N` in the time mesh.
            this.n_y = n_y;                         % ODE state dimension
            this.n_z = n_z;                         % control dimension
            this.n_t = n_t;                         % Number of time nodes
            this.t_mesh = linspace(0, T, n_t)';     % Discrete time domain
            this.time_step_solver_options = optimoptions('fsolve', ...
                                                         'Display', 'none', ...
                                                         'SpecifyObjectiveGradient', true);
            this.verbose = true;                    % Output verbosity
        end

        %% Finite difference tests

        function [diffs_y, diffs_z] = Time_Instance_RHS_Jacobian_Check(this, y, z, t)
            % Check the implementation of :meth:`Time_Instance_RHS()`
            % via finite differences.
            %
            % Parameters
            % ----------
            % y
            %   Differential equation state :math:`\y(t)\in\R^{n_y}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            % t
            %   Time instance :math:`t \in [0, T]`.
            %
            % Returns
            % -------
            % diffs_y : vector
            %   Finite difference errors for the Jacobian :math:`\f_{y}(\y(t),\z,t)`.
            % diffs_z : vector
            %   Finite difference errors for the Jacobian :math:`\f_{z}(\y(t),\z,t)`.

            [f, f_y, f_z] = this.Time_Instance_RHS(y, z, t);

            % Check f_y against finite differences of f.
            v = randn(this.n_y, 1);
            v = v / norm(v);
            fv = f_y * v;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_fv = zeros(this.n_y, p);
            diffs_y = zeros(p, 1);
            for k = 1:p
                [fk] = this.Time_Instance_RHS(y + h(k) * v, z, t);
                fd_fv(:, k) = (fk - f) / h(k);
                diffs_y(k) = norm(fd_fv(:, k) - fv) / norm(fv);
            end
            if this.verbose
                disp('State Jacobian finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_y(k))]);
                end
                disp(' ');
            end

            % Check f_z against finite differences of f.
            v = randn(this.n_z, 1);
            v = v / norm(v);
            fv = f_z * v;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_fv = zeros(this.n_y, p);
            diffs_z = zeros(p, 1);
            for k = 1:p
                [fk] = this.Time_Instance_RHS(y, z + h(k) * v, t);
                fd_fv(:, k) = (fk - f) / h(k);
                diffs_z(k) = norm(fd_fv(:, k) - fv) / norm(fv);
            end
            if this.verbose
                disp('Control Jacobian finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_z(k))]);
                end
                disp(' ');
            end
        end

        function [diffs_yy, diffs_yz, diffs_zy, diffs_zz] = Time_Instance_RHS_Hessian_Check(this, y, z, t)
            % Check the implementation of the following via finite differences.
            %
            % * :meth:`Time_Instance_RHS_yy_Apply()` for :math:`\bflambda\trp\f_{y,y}(\y,\z,t)\v`.
            % * :meth:`Time_Instance_RHS_yz_Apply()` for :math:`\bflambda\trp\f_{y,z}(\y,\z,t)\v`.
            % * :meth:`Time_Instance_RHS_zy_Apply()` for :math:`\bflambda\trp\f_{z,y}(\y,\z,t)\v`.
            % * :meth:`Time_Instance_RHS_zz_Apply()` for :math:`\bflambda\trp\f_{z,z}(\y,\z,t)\v`.
            %
            % Parameters
            % ----------
            % y
            %   Differential equation state :math:`\y\in\R^{n_y}`.
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            % t
            %   Time instance :math:`t \in [0, T]`.
            %
            % Returns
            % -------
            % diffs_yy : vector
            %   Finite difference errors for :math:`\bflambda\trp\f_{y,y}(\y,\z,t)\v`.
            % diffs_yz : vector
            %   Finite difference errors for :math:`\bflambda\trp\f_{y,z}(\y,\z,t)\v`.
            % diffs_zy : vector
            %   Finite difference errors for :math:`\bflambda\trp\f_{z,y}(\y,\z,t)\v`.
            % diffs_zz : vector
            %   Finite difference errors for :math:`\bflambda\trp\f_{z,z}(\y,\z,t)\v`.

            [~, f_y, f_z] = this.Time_Instance_RHS(y, z, t);
            lambda = randn(this.n_y, 1);
            h = 10.^(-2:-1:-6);
            p = length(h);

            % Check Time_Instance_RHS_yy_Apply().
            v = randn(this.n_y, 1);
            v = v / norm(v);
            Mv = this.Time_Instance_RHS_yy_Apply(v, y, z, t, lambda);
            fd_Mv = zeros(this.n_y, p);
            diffs_yy = zeros(p, 1);
            for k = 1:p
                [~, f_yk] = this.Time_Instance_RHS(y + h(k) * v, z, t);
                fd_Mv(:, k) = (f_yk' * lambda - f_y' * lambda) / h(k);
                diffs_yy(k) = norm(fd_Mv(:, k) - Mv) / norm(Mv);
            end
            if this.verbose
                disp('Hessian yy finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_yy(k))]);
                end
                disp(' ');
            end

            % Check Time_Instance_RHS_yz_Apply().
            v = randn(this.n_z, 1);
            v = v / norm(v);
            Mv = this.Time_Instance_RHS_yz_Apply(v, y, z, t, lambda);
            fd_Mv = zeros(this.n_y, p);
            diffs_yz = zeros(p, 1);
            for k = 1:p
                [~, f_yk] = this.Time_Instance_RHS(y, z + h(k) * v, t);
                fd_Mv(:, k) = (f_yk' * lambda - f_y' * lambda) / h(k);
                diffs_yz(k) = norm(fd_Mv(:, k) - Mv) / norm(Mv);
            end
            if this.verbose
                disp('Hessian yz finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_yz(k))]);
                end
                disp(' ');
            end

            % Check Time_Instance_RHS_zy_Apply().
            v = randn(this.n_y, 1);
            v = v / norm(v);
            Mv = this.Time_Instance_RHS_zy_Apply(v, y, z, t, lambda);
            fd_Mv = zeros(this.n_z, p);
            diffs_zy = zeros(p, 1);
            for k = 1:p
                [~, ~, f_zk] = this.Time_Instance_RHS(y + h(k) * v, z, t);
                fd_Mv(:, k) = (f_zk' * lambda - f_z' * lambda) / h(k);
                diffs_zy(k) = norm(fd_Mv(:, k) - Mv) / norm(Mv);
            end
            if this.verbose
                disp('Hessian zy finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_zy(k))]);
                end
                disp(' ');
            end

            % Check Time_Instance_RHS_zz_Apply().
            v = randn(this.n_z, 1);
            v = v / norm(v);
            Mv = this.Time_Instance_RHS_zz_Apply(v, y, z, t, lambda);
            fd_Mv = zeros(this.n_z, p);
            diffs_zz = zeros(p, 1);
            for k = 1:p
                [~, ~, f_zk] = this.Time_Instance_RHS(y, z + h(k) * v, t);
                fd_Mv(:, k) = (f_zk' * lambda - f_z' * lambda) / h(k);
                diffs_zz(k) = norm(fd_Mv(:, k) - Mv) / norm(Mv);
            end
            if this.verbose
                disp('Hessian zz finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs_zz(k))]);
                end
                disp(' ');
            end
        end

    end

    %% Time stepping functions

    methods (Access = protected)

        % Input:
        % ykm: the ODE state y(t_{k-1}) in R^{n_y}
        % z: the control z in R^{n_z}
        % tk: next time step t_k
        % tdm: previous time step t_{k-1}
        % Output:
        % yk: the state y(t_k) in R^{n_y}
        function [yk] = State_Eq_Time_Step(this, ykm, z, tk, tkm)
            dt = tk - tkm;
            yk = fsolve(@(y)this.Nonlinear_Step(y, ykm, z, tk, dt), ykm, this.time_step_solver_options);
        end

        % Input:
        % v: the direction v in R^{n_y}
        % y: the ODE state y in R^{n_y}
        % z: the control z in R^{n_z}
        % tk: the time in the interval [0, T]
        % dt: time step size
        % Output:
        % Mv: (I_{n_y} - dt*f_y(y, z, t_k))^{-1}v
        function [Mv] = Linearized_Time_Step_Solve(this, v, y, z, tk, dt)
            [~, f_y] = this.Time_Instance_RHS(y, z, tk);
            A = eye(this.n_y) - dt * f_y;
            Mv = linsolve(A, v);
        end

        % Input:
        % v: the direction v in R^{n_y}
        % y: the ODE state y in R^{n_y}
        % z: the control z in R^{n_z}
        % tk: the time in the interval [0, T]
        % dt: time step size
        % Output:
        % Mv: (I - dt*f_y(y, z, t_k)^T)^{-1}v
        function [Mv] = Linearized_Adjoint_Time_Step_Solve(this, v, y, z, tk, dt)
            [~, f_y] = this.Time_Instance_RHS(y, z, tk);
            A = eye(this.n_y) - dt * f_y';
            Mv = linsolve(A, v);
        end

        % Input:
        % yk: the ODE state y(t_k) in R^{n_y}
        % ykm: the ODE state y(t_{k-1}) in R^{n_y}
        % z: the control z in R^{n_z}
        % tk: the time in the interval [0, T]
        % dt: time difference t_k-t_{k-1}
        % Output:
        % f: value of the residual y_k-y_{k-1} - dt f(y_k, z, t_k) in R^{n_y}
        % Jac: state Jacobian of the residual in R^{m x m}
        function [f, Jac] = Nonlinear_Step(this, yk, ykm, z, tk, dt)
            [val, val_y] = this.Time_Instance_RHS(yk, z, tk);
            f = yk - ykm - dt * val;
            Jac = eye(this.n_y, this.n_y) - dt * val_y;
        end

    end

end
