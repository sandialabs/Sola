% Solve the optimization problem
% \min_z J(S(z), z)
% where
% S(z) solves the constraint equation c(u, z)=0
% i.e. c(S(z), z)=0 for all z
% where
% u in R^{n_u}
% z in R^{n_z}
% c(u, z) in R^{n_u}

classdef Reduced_Space_Optimization < handle

    % Member properties are the default optimizer settings for iterations and tolerances
    properties
        obj
        con
        opt_tol
        fun_tol
        iteration_limit
        step_tol
        max_cg_iter
        cg_tol
        verbose
        Gauss_Newton_Hess
        use_trust_region
        z_lb
        z_ub
    end

    methods (Access = public)

        % Input:
        % obj: class of type Objective
        % con: class of type Constraint
        function this = Reduced_Space_Optimization(obj, con)
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
        % Input:
        % z0: initial iterate control z in R^{n_z}
        % Output:
        % u: the optimal state solution u in R^{n_u}
        % z: the optimal control z in R^{n_z}
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

        % Input:
        % z: the control z in R^{n_z}
        % Output:
        % val: \hat{J}(z)
        % grad: \nabla_z \hat{J}(z)
        % hessian_data: concatenation of state, control, and adjoint to pass to hessian-vector multiply function
        function [val, grad, hessian_data] = Jhat(this, z)
            u = this.con.State_Solve(z);
            [val, grad_u, grad_z] = this.obj.J(u, z);
            lambda = this.con.c_u_Transpose_Inverse_Apply(-grad_u, u, z);
            grad = this.con.c_z_Transpose_Apply(lambda, u, z);
            grad = grad + grad_z;
            hessian_data = [u; z; lambda];
        end

        % Input:
        % hessian_data: output from Jhat function containing the state u, control z, and adjoint lambda
        % v: a direction v in R^{n_z}
        % Output:
        % Hv: \nabla_{z, z} \hat{J}(z)v
        function [Hv] = Jhat_hessVec(this, hessian_data, v)
            p = length(v);
            m = (length(hessian_data) - p) / 2;
            u = hessian_data(1:m);
            z = hessian_data((m + 1):(m + p));
            lambda = hessian_data((m + p + 1):end);

            w = this.con.c_z_Apply(v, u, z);
            mu = this.con.c_u_Inverse_Apply(-w, u, z);
            yJ = this.obj.J_uu_Apply(mu, u, z) + this.obj.J_uz_Apply(v, u, z);
            xJ = this.obj.J_zu_Apply(mu, u, z) + this.obj.J_zz_Apply(v, u, z);
            if this.Gauss_Newton_Hess
                gamma = this.con.c_u_Transpose_Inverse_Apply(-yJ, u, z);
                xc = this.con.c_z_Transpose_Apply(gamma, u, z);
            else
                yc = this.con.c_uu_Apply(mu, u, z, lambda) + this.con.c_uz_Apply(v, u, z, lambda);
                gamma = this.con.c_u_Transpose_Inverse_Apply(-(yJ + yc), u, z);
                xc = this.con.c_z_Transpose_Apply(gamma, u, z) + this.con.c_zu_Apply(mu, u, z, lambda) + this.con.c_zz_Apply(v, u, z, lambda);
            end
            Hv = xJ + xc;
        end

        %% Finite difference test functions
        % Input:
        % z: the control z in R^{n_z}
        % Output:
        % diffs: vector of finite difference errors
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

        % Input:
        % z: the control z in R^{n_z}
        % Output:
        % diffs: vector of finite difference errors
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
