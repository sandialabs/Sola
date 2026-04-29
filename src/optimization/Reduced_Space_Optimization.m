%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Reduced_Space_Optimization < handle
    % Solve a constrained optimization problem
    %  \min_{u,z} J(u,z)
    %  subject to
    %  c(u,z) = 0
    % by solving the equivalent, unconstrained reduced-space problem
    % min_z \hat{J}(z) := J(S(z),z)

    properties
        obj                 % Instance of a subclass of :class:`Objective`.
        con                 % Instance of a subclass of :class:`Constraint`.
        opt_tol             % Optimality tolerance.
        fun_tol             % Function tolerance.
        iteration_limit     % Maximum number of iterations.
        step_tol            % Step tolerance.
        max_cg_iter         % Maximum number of conjugate gradient iterations.
        cg_tol              % Conjugate gradient tolerance.
        verbose             % Verbosity.
        Gauss_Newton_Hess   % Use Gauss-Newton approximation of the Hessian.
        use_trust_region    % Use trust region for the optimization.
        z_lb                % Lower bounds for control.
        z_ub                % Upper bounds for control.
    end

    methods (Access = public)

        %% Constructor

        function this = Reduced_Space_Optimization(obj, con)
            arguments
                obj Objective
                con Constraint
            end

            this.obj = obj;
            this.con = con;
            this.opt_tol = 10^-8;
            this.fun_tol = 10^-6;
            this.iteration_limit = 10^3;
            this.step_tol = 10^-6;
            this.max_cg_iter = 50;
            this.cg_tol = 10^-4;
            this.verbose = true;
            this.Gauss_Newton_Hess = false;
            this.use_trust_region = true;
            this.z_lb = [];
            this.z_ub = [];
        end

        %% Optimization functions

        function [u, z] = Optimize(this, z0)

            verb = 'iter-detailed';
            if this.verbose == false
                verb = 'none';
            end
            if this.use_trust_region
                HessMultFcn = @(hessian_data, v) this.Jhat_hessVec(hessian_data, v);
                if length(this.z_lb) + length(this.z_ub) == 0
                    options = optimoptions(@fminunc, ...
                                           'Display', verb, ...
                                           'Algorithm', 'trust-region', ...
                                           'SpecifyObjectiveGradient', true, ...
                                           'OptimalityTolerance', this.opt_tol, ...
                                           'FunctionTolerance', this.fun_tol, ...
                                           'MaxIterations', this.iteration_limit, ...
                                           'StepTolerance', this.step_tol, ...
                                           'SubproblemAlgorithm', 'cg', ...
                                           'MaxPCGIter', this.max_cg_iter, ...
                                           'TolPCG', this.cg_tol, ...
                                           'HessianMultiplyFcn', HessMultFcn);
                    z = fminunc(@(z)this.Jhat(z), z0, options);
                else
                    options = optimoptions(@fmincon, ...
                                           'Display', verb, ...
                                           'Algorithm', 'trust-region-reflective', ...
                                           'SpecifyObjectiveGradient', true, ...
                                           'OptimalityTolerance', this.opt_tol, ...
                                           'FunctionTolerance', this.fun_tol, ...
                                           'MaxIterations', this.iteration_limit, ...
                                           'StepTolerance', this.step_tol, ...
                                           'SubproblemAlgorithm', 'cg', ...
                                           'MaxPCGIter', this.max_cg_iter, ...
                                           'TolPCG', this.cg_tol, ...
                                           'HessianMultiplyFcn', HessMultFcn);
                    z = fmincon(@(z)this.Jhat(z), z0, [], [], [], [], this.z_lb, this.z_ub, [], options);
                end
                u = this.con.State_Solve(z);
            else
                if length(this.z_lb) + length(this.z_ub) == 0
                    options = optimoptions(@fminunc, ...
                                           'Display', verb, ...
                                           'SpecifyObjectiveGradient', true, ...
                                           'OptimalityTolerance', this.opt_tol, ...
                                           'FunctionTolerance', this.fun_tol, ...
                                           'MaxIterations', this.iteration_limit, ...
                                           'StepTolerance', this.step_tol);
                    z = fminunc(@(z)this.Jhat(z), z0, options);
                    u = this.con.State_Solve(z);
                else
                    options = optimoptions(@fmincon, ...
                                           'Display', verb, ...
                                           'SpecifyObjectiveGradient', true, ...
                                           'OptimalityTolerance', this.opt_tol, ...
                                           'FunctionTolerance', this.fun_tol, ...
                                           'MaxIterations', this.iteration_limit, ...
                                           'StepTolerance', this.step_tol);
                    z = fmincon(@(z)this.Jhat(z), z0, [], [], [], [], this.z_lb, this.z_ub, [], options);
                    u = this.con.State_Solve(z);
                end
            end
        end

        function [val, grad, hessian_data] = Jhat(this, z)
            u = this.con.State_Solve(z);
            [val, grad_u, grad_z] = this.obj.J(u, z);
            lambda = this.con.c_u_Transpose_Inverse_Apply(-grad_u, u, z);
            grad = this.con.c_z_Transpose_Apply(lambda, u, z);
            grad = grad + grad_z;
            hessian_data = [u; z; lambda];
        end

        function [z_out] = Jhat_hessVec(this, hessian_data, z_in)
            % hessian_data is the concatenation of the state u,
            %   control z, and adjoint lambda

            % Extract state, control, and adjoint from hessian_data.
            p = length(z_in);
            m = (length(hessian_data) - p) / 2;
            u = hessian_data(1:m);
            z = hessian_data((m + 1):(m + p));
            lambda = hessian_data((m + p + 1):end);

            w = this.con.c_z_Apply(z_in, u, z);
            mu = this.con.c_u_Inverse_Apply(-w, u, z);
            yJ = this.obj.J_uu_Apply(mu, u, z) + this.obj.J_uz_Apply(z_in, u, z);
            xJ = this.obj.J_zu_Apply(mu, u, z) + this.obj.J_zz_Apply(z_in, u, z);
            if this.Gauss_Newton_Hess
                gamma = this.con.c_u_Transpose_Inverse_Apply(-yJ, u, z);
                xc = this.con.c_z_Transpose_Apply(gamma, u, z);
            else
                yc = this.con.c_uu_Apply(mu, u, z, lambda) + this.con.c_uz_Apply(z_in, u, z, lambda);
                gamma = this.con.c_u_Transpose_Inverse_Apply(-(yJ + yc), u, z);
                xc = this.con.c_z_Transpose_Apply(gamma, u, z) + this.con.c_zu_Apply(mu, u, z, lambda) + this.con.c_zz_Apply(z_in, u, z, lambda);
            end
            z_out = xJ + xc;
        end

        %% Finite difference tests

        function [diffs] = Finite_Difference_Gradient_Check(this, z)

            [val, grad] = this.Jhat(z);
            n = length(grad);
            dz = randn(n, 1);
            dz = dz / norm(dz);
            grad_dz = dz' * grad;
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_grad = zeros(p, 1);
            for k = 1:p
                valk = this.Jhat(z + h(k) * dz);
                fd_grad(k) = (valk - val) / h(k);
            end

            diffs = abs(grad_dz - fd_grad) / abs(grad_dz);
            if this.verbose
                disp('Gradient finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

        function [diffs] = Finite_Difference_Hessian_Check(this, z)

            [~, grad, hessian_data] = this.Jhat(z);
            n = length(grad);
            v = randn(n, 1);
            v = v / norm(v);
            Hv = this.Jhat_hessVec(hessian_data, v);
            h = 10.^(-2:-1:-6);
            p = length(h);
            fd_hv = zeros(n, p);
            diffs = zeros(p, 1);
            for k = 1:p
                [~, gradk] = this.Jhat(z + h(k) * v);
                fd_hv(:, k) = (gradk - grad) / h(k);
                diffs(k) = norm(fd_hv(:, k) - Hv) / norm(Hv);
            end
            if this.verbose
                disp('Hessian finite difference check');
                for k = 1:p
                    disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
                end
                disp(' ');
            end
        end

    end
end
