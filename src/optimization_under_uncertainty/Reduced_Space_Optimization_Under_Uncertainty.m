%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Reduced_Space_Optimization_Under_Uncertainty < handle
    % Solve a constrained optimization problem
    %
    % .. math::
    %  \min_{\u_1,\u_2,\dots,\u_N,\z} \frac{1}{N} \sum\limits_{i=1}^N J(\u_i,\z)
    %  \quad \text{subject to} \quad
    %  \c(\u_i,\z,\theta_i) = \0
    %
    % by solving the equivalent, unconstrained reduced-space problem
    %
    % .. math:: \min_\z \hat{J}(\z) := \frac{1}{N} \sum\limits_{i=1}^N J(\S_i(\z), \z)
    %
    % where
    % :math:`\S_i(\z)` solves the constraint equation :math:`\c(\u_i, \z, \theta_i) = \0`,
    % i.e., :math:`\c(\S_i(\z),\z, \theta_i) = \0` for all admissible :math:`\z`.
    % Here,
    % :math:`\u \in \R^{n_u}` is the state,
    % :math:`\z \in \R^{n_z}` is the control,
    % :math:`\theta \in \R^{n_theta}` are the parameters, and
    % :math:`\c(\u_i,\z,\theta_i) \in \R^{n_u}` are the constraints.

    properties
        obj                 % Instance of a subclass of :class:`Objective`.
        cons                % Instance of a cell of subclasses of :class:`Parametric_Constraint`.
        N                   % Number of samples (length of cons)
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

        function this = Reduced_Space_Optimization_Under_Uncertainty(obj, cons)
            % Parameters
            % ----------
            % obj
            %   Objective function, an instance of a subclass of :class:`Objective`.
            % cons
            %   Constraint equations, a cell of subclasses of :class:`Constraint`.

            this.obj = obj;
            this.cons = cons;
            this.N = length(cons);
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
            %   Optimal state :math:`\u^* \in \R^{n_u \times N}`.
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
                u = cell(1, this.N);
                for k = 1:this.N
                    u{k} = this.cons{k}.State_Solve(z);
                end
                u = cell2mat(u);
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
                    u = cell(1, this.N);
                    for k = 1:this.N
                        u{k} = this.cons{k}.State_Solve(z);
                    end
                    u = cell2mat(u);
                else
                    options = optimoptions(@fmincon, ...
                                           'Display', verb, ...
                                           'SpecifyObjectiveGradient', true, ...
                                           'OptimalityTolerance', this.opt_tol, ...
                                           'FunctionTolerance', this.fun_tol, ...
                                           'MaxIterations', this.iteration_limit, ...
                                           'StepTolerance', this.step_tol);
                    z = fmincon(@(z)this.Jhat(z), z0, [], [], [], [], this.z_lb, this.z_ub, [], options);
                    u = cell(1, this.N);
                    for k = 1:this.N
                        u{k} = this.cons{k}.State_Solve(z);
                    end
                    u = cell2mat(u);
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

            val = zeros(1, this.N);
            grad = zeros(length(z), this.N);
            hessian_data = cell(1, this.N);

            for k = 1:this.N
                u = this.cons{k}.State_Solve(z);
                [val(k), grad_u, grad_z] = this.obj.J(u, z);
                lambda = this.cons{k}.c_u_Transpose_Inverse_Apply(-grad_u, u, z);
                grad(:, k) = this.cons{k}.c_z_Transpose_Apply(lambda, u, z);
                grad(:, k) = grad(:, k) + grad_z;
                hessian_data{k} = [u; z; lambda];
            end
            val = mean(val, 2);
            grad = mean(grad, 2);
            hessian_data = cell2mat(hessian_data);
        end

        function [z_out] = Jhat_hessVec(this, hessian_data, z_in)
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
            % z_in : vector
            %   Control direction :math:`\v\in\R^{n_z}`.
            %
            % Returns
            % -------
            % z_out : vector
            %   Vector-Hessian-vector product
            %   :math:`\bflambda\trp\grad{z,z}\hat{J}(\z)\v\in\R^{n_z}`.

            p = size(z_in, 1);
            M = size(z_in, 2);
            m = (length(hessian_data(:, 1)) - p) / 2;
            z_out_tmp = zeros(p, M, this.N);

            for i = 1:M
                for k = 1:this.N
                    % Extract state, control, and adjoint from hessian_data.
                    u = hessian_data(1:m, k);
                    z = hessian_data((m + 1):(m + p), k);
                    lambda = hessian_data((m + p + 1):end, k);

                    % Execute Algorithm 2 or 3 for computing the Hessian-vector product.
                    w = this.cons{k}.c_z_Apply(z_in(:, i), u, z);
                    mu = this.cons{k}.c_u_Inverse_Apply(-w, u, z);
                    yJ = this.obj.J_uu_Apply(mu, u, z) + this.obj.J_uz_Apply(z_in(:, i), u, z);
                    xJ = this.obj.J_zu_Apply(mu, u, z) + this.obj.J_zz_Apply(z_in(:, i), u, z);
                    if this.Gauss_Newton_Hess
                        gamma = this.cons{k}.c_u_Transpose_Inverse_Apply(-yJ, u, z);
                        xc = this.cons{k}.c_z_Transpose_Apply(gamma, u, z);
                    else
                        yc = this.cons{k}.c_uu_Apply(mu, u, z, lambda) + this.cons{k}.c_uz_Apply(z_in(:, i), u, z, lambda);
                        gamma = this.cons{k}.c_u_Transpose_Inverse_Apply(-(yJ + yc), u, z);
                        xc = this.cons{k}.c_z_Transpose_Apply(gamma, u, z) + this.cons{k}.c_zu_Apply(mu, u, z, lambda) + this.cons{k}.c_zz_Apply(z_in(:, i), u, z, lambda);
                    end
                    z_out_tmp(:, i, k) = xJ + xc;
                end
            end
            z_out = mean(z_out_tmp, 3);
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
