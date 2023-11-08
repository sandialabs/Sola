% Define the constraint as an ordinary differential equation
% \frac{dy}{dt} = f(y, z)
% y(0) = h(z)
% where
% T > 0
% y(t) in R^m
% z in R^n
% c(u, z) = 0 corresponds to the discretization of the ODE system.
classdef Dynamic_Constraint < Constraint

    properties
        m                           % dimension of the ODE state y(t)
        n                           % dimension of the control z
        T                           % final time
        N                           % number of nodes in the time mesh
        t_mesh                      % time mesh (a vector of length N)
        w                           % quadrature weights for time integration
        time_step_solver_options    % options set for fsolve in time step
        verbose                     % output verbosity
    end

    methods (Abstract, Access = public)

        % Input:
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % Output:
        % f: f(y(t), z) in R^m
        % f_y: f_y(y(t), z, t) in R^{m \times m}
        % f_z: f_z(y(t), z, t) in R^{m \times n}
        [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)

        % Input
        % z: the control z in R^n
        % Output:
        % h: h(z)in R^m
        % h_z: h_z(z)in R^{m \times n}
        [h, h_z] = Initial_Condition(this, z)

        % Input:
        % v: a direction v in R^m
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % lambda: the adjoint state \lambda(t) in R^m
        % Output:
        % Mv: \lambda(t)^T f_{y, y}(y(t), z, t)v in R^m
        [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)

        % Input:
        % v: a direction v in R^n
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % lambda: the adjoint state \lambda(t) in R^m
        % Output:
        % Mv: \lambda(t)^T f_{y, z}(y(t), z, t)v in R^m
        [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)

        % Input:
        % v: a direction v in R^m
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % lambda: the adjoint state \lambda(t) in R^m
        % Output:
        % Mv: \lambda(t)^T f_{z, y}(y(t), z, t)v in R^n
        [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)

        % Input:
        % v: a direction v in R^n
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % lambda: the adjoint state \lambda(t) in R^m
        % Output:
        % Mv: \lambda(t)^T f_{z, z}(y(t), z)v in R^n
        [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)

        % Input:
        % v: a direction v in R^n
        % z: the control z in R^n
        % lambda: the adjoint state \lambda(t) in R^m
        % Output:
        % Mv: \lambda(t)^T h_{z, z}(z)v in R^n
        [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)

    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u = zeros(this.m * this.N, 1);
            u(1:this.m) = this.Initial_Condition(z);
            for k = 2:this.N
                Im = ((k - 2) * this.m + 1):((k - 1) * this.m);      % y_{k-1} = u(Im)
                I = ((k - 1) * this.m + 1):(k * this.m);           % y_{k} = u(I)
                u(I) = this.State_Eq_Time_Step(u(Im), z, this.t_mesh(k), this.t_mesh(k - 1));
            end
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            num_vecs = size(v, 2);
            Mv = zeros(this.m * this.N, num_vecs);
            I = ((this.N - 1) * this.m + 1):(this.N * this.m);      % y_{k} = u(I)
            dt = this.t_mesh(end) - this.t_mesh(end - 1);
            Mv(I, :) = this.Linearized_Adjoint_Time_Step_Solve(v(I, :), u(I), z, this.t_mesh(end), dt);

            for k = (this.N - 1):-1:2
                Im = (k * this.m + 1):((k + 1) * this.m);         % y_{k} = u(Im)
                I = ((k - 1) * this.m + 1):(k * this.m);          % y_{k-1} = u(I)
                dt = this.t_mesh(k + 1) - this.t_mesh(k);
                Mv(I, :) = this.Linearized_Adjoint_Time_Step_Solve(v(I, :) + Mv(Im, :), u(I), z, this.t_mesh(k), dt);
            end
            I = 1:(this.m);
            Im = (this.m + 1):(2 * this.m);
            Mv(I, :) = Mv(Im, :) + v(I, :);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            [~, h_z] = this.Initial_Condition(z);
            Mv = -h_z' * v(1:this.m, :);
            for k = 2:this.N
                I = ((k - 1) * this.m + 1):(k * this.m);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                [~, ~, f_z] = this.Time_Instance_RHS(u(I), z, this.t_mesh(k));
                Mv = Mv - dt * f_z' * v(I, :);
            end
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            num_vecs = size(v, 2);
            Mv = zeros(this.m * this.N, num_vecs);
            I = 1:this.m;
            Mv(I, :) = v(I, :);
            for k = 2:this.N
                Im = ((k - 2) * this.m + 1):((k - 1) * this.m);     % y_{k-1} = u(Im)
                I = ((k - 1) * this.m + 1):(k * this.m);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                Mv(I, :) = this.Linearized_Time_Step_Solve(v(I, :) + Mv(Im, :), u(I), z, this.t_mesh(k), dt);
            end
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            num_vecs = size(v, 2);
            Mv = zeros(this.m * this.N, num_vecs);
            I = 1:this.m;
            [~, h_z] = this.Initial_Condition(z);
            Mv(I, :) = -h_z * v;
            for k = 2:this.N
                I = ((k - 1) * this.m + 1):(k * this.m);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                [~, ~, f_z] = this.Time_Instance_RHS(u(I), z, this.t_mesh(k));
                Mv(I, :) = Mv(I, :) - dt * f_z * v;
            end
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.m * this.N, num_vecs);
            for k = 2:this.N
                I = ((k - 1) * this.m + 1):(k * this.m);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                f_yy = this.Time_Instance_RHS_yy_Apply(v(I, :), u(I), z, this.t_mesh(k), lambda(I));
                Mv(I, :) = -dt * f_yy;
            end
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.m * this.N, num_vecs);
            for k = 2:this.N
                I = ((k - 1) * this.m + 1):(k * this.m);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                f_yz = this.Time_Instance_RHS_yz_Apply(v, u(I), z, this.t_mesh(k), lambda(I));
                Mv(I, :) = -dt * f_yz;
            end
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n, num_vecs);
            for k = 2:this.N
                I = ((k - 1) * this.m + 1):(k * this.m);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                f_zy = this.Time_Instance_RHS_zy_Apply(v(I, :), u(I), z, this.t_mesh(k), lambda(I));
                Mv = Mv - dt * f_zy;
            end
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = -this.Initial_Condition_zz_Apply(v, z, lambda(1:this.m));
            for k = 2:this.N
                I = ((k - 1) * this.m + 1):(k * this.m);          % y_{k} = u(I)
                dt = this.t_mesh(k) - this.t_mesh(k - 1);
                f_zz = this.Time_Instance_RHS_zz_Apply(v, u(I), z, this.t_mesh(k), lambda(I));
                Mv = Mv - dt * f_zz;
            end
        end

    end

    % Input:
    % m: the dimension of the ODE state y(t)
    % n: the dimension of the control z
    % T: the final time
    % N: the number of nodes in the time mesh
    methods (Access = public)

        function this = Dynamic_Constraint(m, n, T, N)
            this.m = m;                         % ODE state dimension
            this.n = n;                         % control dimension
            this.T = T;                         % Final time (intial time is 0)
            this.N = N;                         % Number of time nodes
            this.t_mesh = linspace(0, T, N)';   % Discrete time domain
            w = ones(N, 1);
            w(2:end - 1) = 2;
            this.w = T * w / sum(w);            % Spatial weights
            this.time_step_solver_options = optimoptions('fsolve', ...
                                                         'Display', 'none', ...
                                                         'SpecifyObjectiveGradient', true);
            this.verbose = true;
        end

        %% Finite difference tests
        % Input:
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % Output:
        % diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Jacobian_y_Check(this, y, z, t)
            [f, f_y] = this.Time_Instance_RHS(y, z, t);
            v = randn(this.m, 1);
            v = v / norm(v);
            fv = f_y * v;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_fv = zeros(this.m, p);
            diffs = zeros(p, 1);
            for k = 1:p
                [fk] = this.Time_Instance_RHS(y + h(k) * v, z, t);
                fd_fv(:, k) = (fk - f) / h(k);
                diffs(k) = norm(fd_fv(:, k) - fv) / norm(fv);
            end
            if this.verbose
                disp('State Jacobian finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

        % Input:
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % Output:
        % diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Jacobian_z_Check(this, y, z, t)
            [f, ~, f_z] = this.Time_Instance_RHS(y, z, t);
            v = randn(this.n, 1);
            v = v / norm(v);
            fv = f_z * v;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_fv = zeros(this.m, p);
            diffs = zeros(p, 1);
            for k = 1:p
                [fk] = this.Time_Instance_RHS(y, z + h(k) * v, t);
                fd_fv(:, k) = (fk - f) / h(k);
                diffs(k) = norm(fd_fv(:, k) - fv) / norm(fv);
            end
            if this.verbose
                disp('Control Jacobian finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

        % Input:
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % Output:
        % diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Hessian_yy_Check(this, y, z, t)
            v = randn(this.m, 1);
            v = v / norm(v);
            lambda = randn(this.m, 1);
            Mv = this.Time_Instance_RHS_yy_Apply(v, y, z, t, lambda);
            [~, f_y] = this.Time_Instance_RHS(y, z, t);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Mv = zeros(this.m, p);
            diffs = zeros(p, 1);
            for k = 1:p
                [~, f_yk] = this.Time_Instance_RHS(y + h(k) * v, z, t);
                fd_Mv(:, k) = (f_yk' * lambda - f_y' * lambda) / h(k);
                diffs(k) = norm(fd_Mv(:, k) - Mv) / norm(Mv);
            end
            if this.verbose
                disp('Hessian yy finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

        % Input:
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % Output:
        % diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Hessian_yz_Check(this, y, z, t)
            v = randn(this.n, 1);
            v = v / norm(v);
            lambda = randn(this.m, 1);
            Mv = this.Time_Instance_RHS_yz_Apply(v, y, z, t, lambda);
            [~, f_y] = this.Time_Instance_RHS(y, z, t);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Mv = zeros(this.m, p);
            diffs = zeros(p, 1);
            for k = 1:p
                [~, f_yk] = this.Time_Instance_RHS(y, z + h(k) * v, t);
                fd_Mv(:, k) = (f_yk' * lambda - f_y' * lambda) / h(k);
                diffs(k) = norm(fd_Mv(:, k) - Mv) / norm(Mv);
            end
            if this.verbose
                disp('Hessian yz finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

        % Input:
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % Output:
        % diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Hessian_zy_Check(this, y, z, t)
            v = randn(this.m, 1);
            v = v / norm(v);
            lambda = randn(this.m, 1);
            Mv = this.Time_Instance_RHS_zy_Apply(v, y, z, t, lambda);
            [~, ~, f_z] = this.Time_Instance_RHS(y, z, t);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Mv = zeros(this.n, p);
            diffs = zeros(p, 1);
            for k = 1:p
                [~, ~, f_zk] = this.Time_Instance_RHS(y + h(k) * v, z, t);
                fd_Mv(:, k) = (f_zk' * lambda - f_z' * lambda) / h(k);
                diffs(k) = norm(fd_Mv(:, k) - Mv) / norm(Mv);
            end
            if this.verbose
                disp('Hessian zy finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

        % Input:
        % y: the ODE state y(t) in R^m
        % z: the control z in R^n
        % t: the time in the interval [0, T]
        % Output:
        % diffs: vector of finite difference errors
        function [diffs] = Time_Instance_RHS_Hessian_zz_Check(this, y, z, t)
            v = randn(this.n, 1);
            v = v / norm(v);
            lambda = randn(this.m, 1);
            Mv = this.Time_Instance_RHS_zz_Apply(v, y, z, t, lambda);
            [~, ~, f_z] = this.Time_Instance_RHS(y, z, t);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_Mv = zeros(this.n, p);
            diffs = zeros(p, 1);
            for k = 1:p
                [~, ~, f_zk] = this.Time_Instance_RHS(y, z + h(k) * v, t);
                fd_Mv(:, k) = (f_zk' * lambda - f_z' * lambda) / h(k);
                diffs(k) = norm(fd_Mv(:, k) - Mv) / norm(Mv);
            end
            if this.verbose
                disp('Hessian zz finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

    end

    %% Time stepping functions
    methods (Access = protected)

        % Input:
        % ykm: the ODE state y(t_{k-1}) in R^m
        % z: the control z in R^n
        % tk: next time step t_k
        % tdm: previous time step t_{k-1}
        % Output:
        % yk: the state y(t_k) in R^m
        function [yk] = State_Eq_Time_Step(this, ykm, z, tk, tkm)
            dt = tk - tkm;
            yk = fsolve(@(y)this.Nonlinear_Step(y, ykm, z, tk, dt), ykm, this.time_step_solver_options);
        end

        % Input:
        % v: the direction v in R^m
        % y: the ODE state y in R^m
        % z: the control z in R^n
        % tk: the time in the interval [0, T]
        % dt: time step size
        % Output:
        % Mv: (I_m - dt*f_y(y, z, t_k))^{-1}v
        function [Mv] = Linearized_Time_Step_Solve(this, v, y, z, tk, dt)
            [~, f_y] = this.Time_Instance_RHS(y, z, tk);
            A = eye(this.m) - dt * f_y;
            Mv = linsolve(A, v);
        end

        % Input:
        % v: the direction v in R^m
        % y: the ODE state y in R^m
        % z: the control z in R^n
        % tk: the time in the interval [0, T]
        % dt: time step size
        % Output:
        % Mv: (I_m - dt*f_y(y, z, t_k)^T)^{-1}v
        function [Mv] = Linearized_Adjoint_Time_Step_Solve(this, v, y, z, tk, dt)
            [~, f_y] = this.Time_Instance_RHS(y, z, tk);
            A = eye(this.m) - dt * f_y';
            Mv = linsolve(A, v);
        end

        % Input:
        % yk: the ODE state y(t_k) in R^m
        % ykm: the ODE state y(t_{k-1}) in R^m
        % z: the control z in R^n
        % tk: the time in the interval [0, T]
        % dt: time difference t_k-t_{k-1}
        % Output:
        % f: value of the residual y_k-y_{k-1} - dt f(y_k, z, t_k) in R^m
        % Jac: state Jacobian of the residual in R^{m x m}
        function [f, Jac] = Nonlinear_Step(this, yk, ykm, z, tk, dt)
            [val, val_y] = this.Time_Instance_RHS(yk, z, tk);
            f = yk - ykm - dt * val;
            Jac = eye(this.m, this.m);
            Jac = Jac - dt * val_y;
        end

    end

end
