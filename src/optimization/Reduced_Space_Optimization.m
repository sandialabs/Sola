classdef Reduced_Space_Optimization < handle
    % Solve a constrained optimization problem
    %
    % .. math::
    %  \min_{\u,\z} J(\u,\z)
    %  \quad \text{subject to} \quad
    %  \c(\u,\z) = \0
    %
    % by solving the equivalent, unconstrained reduced-space problem
    %
    % .. math:: \min_\z \hat{J}(\z) := J(\S(\z), \z)
    %
    % where
    % :math:`\S(\z)` solves the constraint equation :math:`\c(\u, \z) = \0`,
    % i.e., :math:`\c(\S(\z),\z) = \0` for all admissible :math:`\z`.
    % Here,
    % :math:`\u \in \R^{n_u}` is the state,
    % :math:`\z \in \R^{n_z}` is the control, and
    % :math:`\c(\u,\z) \in \R^{n_u}` are the constraints.

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
            % Parameters
            % ----------
            % obj
            %   Objective function, an instance of a subclass of :class:`Objective`.
            % con
            %   Constraint equations, an instance of a subclass of :class:`Constraint`.

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
            % Solve the optimization problem via ``fminunc()``
            % (or ``fmincon()`` if ``z_lb`` and ``z_ub`` are set).
            % Optimization options are given by the object properties.
            %
            % Parameters
            % ----------
            % z0
            %   Initial guess :math:`\z_0\in\R^{n_z}` for the control.
            %
            % Returns
            % -------
            % u : vector
            %   Optimal state :math:`\u^* \in \R^{n_u}`.
            % z : vetor
            %   Optimal control :math:`\z^* \in \R^{n_z}`.

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
            % Evaluate the reduced-space objective function
            % :math:`\hat{J}(\z) = J(\S(\z),\z)` and its gradients.
            %
            % Parameters
            % ----------
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            %
            % Returns
            % -------
            % val : double
            %   Objective value :math:`\hat{J}(\z)\in\R`.
            % grad : vetor
            %   Objective gradient :math:`\grad{z}\hat{J}(\z)\in\R^{n_z}`.
            % hessian_data : vector
            %   Concatenation of the state :math:`\u\in\R^{n_u}`,
            %   control :math:`\z\in\R^{n_z}`,
            %   and adjoint :math:`\bflambda\in\R^{n_u}`
            %   to pass to :meth:`Jhat_hessVec()`.

            u = this.con.State_Solve(z);
            [val, grad_u, grad_z] = this.obj.J(u, z);
            lambda = this.con.c_u_Transpose_Inverse_Apply(-grad_u, u, z);
            grad = this.con.c_z_Transpose_Apply(lambda, u, z);
            grad = grad + grad_z;
            hessian_data = [u; z; lambda];
        end

        function [Hv] = Jhat_hessVec(this, hessian_data, v)
            % Compute the vector-Hessian-vector product
            % :math:`\bflambda\trp\grad{z,z}\hat{J}(\z)\v`
            % via an adjoint-based approach.
            %
            % Parameters
            % ----------
            % hessian_data : vector
            %   Concatenation of the state :math:`\u\in\R^{n_u}`,
            %   control :math:`\z\in\R^{n_z}`,
            %   and adjoint :math:`\bflambda\in\R^{n_u}`.
            % v : vector
            %   Control direction :math:`\v\in\R^{n_z}`.
            %
            % Returns
            % -------
            % Hv : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp\grad{z,z}\hat{J}(\z)\v\in\R^{n_z}`.

            % Extract state, control, and adjoint from hessian_data.
            p = length(v);
            m = (length(hessian_data) - p) / 2;
            u = hessian_data(1:m);
            z = hessian_data((m + 1):(m + p));
            lambda = hessian_data((m + p + 1):end);

            % Execute Algorithm 2 or 3 for computing the Hessian-vector product.
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

        %% Finite difference tests

        function [diffs] = Finite_Difference_Gradient_Check(this, z)
            % Check the implementation of :meth:`Jhat()` via finite differences.
            %
            % Parameters
            % ----------
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            %
            % Returns
            % -------
            % diffs : vector
            %   Finite difference errors.

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
            % Check the implementation of :meth:`Jhat_hessVec()` via finite differences.
            %
            % Parameters
            % ----------
            % z
            %   Control :math:`\z\in\R^{n_z}`.
            %
            % Returns
            % -------
            % diffs : vector
            %   Finite difference errors.

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
