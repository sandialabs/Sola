classdef Rosenbrock < handle

    properties
        d
        coeff
    end

    methods

        function this = Rosenbrock(d)
            this.d = d;
            this.coeff = 1;
        end

        function [val, grad, hess] = J(this, z, theta)
            if length(theta) ~= this.d - 1
                disp('Error in specification of theta');
            end

            val = 0;
            grad = zeros(this.d, 1);
            hess = zeros(this.d, this.d);

            for i = 1:(this.d - 1)
                val = val + (z(i) - theta(i))^2 + this.coeff * (z(i + 1) - z(i)^2)^2;
                grad(i) = 2 * (z(i) - theta(i)) - 4 * this.coeff * (z(i + 1) - z(i)^2) * z(i);
                hess(i, i) = 2 - 4 * this.coeff * (z(i + 1) - 3 * z(i)^2);
                hess(i, i + 1) = -4 * this.coeff * z(i);
                if i > 1
                    grad(i) = grad(i) + 2 * this.coeff * (z(i) - z(i - 1)^2);
                    hess(i, i) = hess(i, i) + 2 * this.coeff;
                    hess(i, i - 1) = -4 * this.coeff * z(i - 1);
                end
            end
            grad(this.d) = 2 * this.coeff * (z(this.d) - z(this.d - 1)^2);
            hess(this.d, this.d) = 2 * this.coeff;
            hess(this.d, this.d - 1) = -4 * this.coeff * z(this.d - 1);
        end

        function [B] = Compute_B(this, z, theta)
            if length(theta) ~= this.d - 1
                disp('Error in specification of theta');
            end

            B = -2 * eye(this.d);
            B = B(:, 1:(this.d - 1));
        end

        function [z_opt] = Solve_Optimization_Problem(this, z0, theta, verbose)
            if nargin < 4
                verbose = 'none';
            end
            objective_fun = @(x) this.J(x, theta);
            options = optimoptions('fminunc', 'Algorithm', 'trust-region', 'SpecifyObjectiveGradient', true, 'HessianFcn', 'objective', 'Display', verbose, 'OptimalityTolerance', 1e-10, 'FunctionTolerance', 10^-10, 'MaxIterations', 5000);
            z_opt = fminunc(objective_fun, z0, options);
        end

        % Finite difference test for the gradient
        function [] = Gradient_FD_Check(this, z, theta)
            [val, grad] = this.J(z, theta);
            h = 10.^(-1:-1:-6);
            v = randn(length(z), 1);
            error = zeros(length(h), 1);
            for k = 1:length(h)
                valk = this.J(z + h(k) * v, theta);
                error(k) = abs((valk - val) / h(k) - grad' * v) / abs(grad' * v);
            end

            disp('Gradient finite difference test');
            for k = 1:length(h)
                disp(['Step size = ', num2str(h(k)), ' and error = ', num2str(error(k))]);
            end
        end

        % Finite difference test for the hessian
        function [] = Hessian_FD_Check(this, z, theta)
            [~, grad, hess] = this.J(z, theta);
            h = 10.^(-1:-1:-6);
            v = randn(length(z), 1);
            Hv = hess * v;
            error = zeros(length(h), 1);
            for k = 1:length(h)
                [~, gradk] = this.J(z + h(k) * v, theta);
                error(k) = norm((gradk - grad) / h(k) - Hv) / norm(Hv);
            end

            disp('Hessian finite difference test');
            for k = 1:length(h)
                disp(['Step size = ', num2str(h(k)), ' and error = ', num2str(error(k))]);
            end
        end

        % Finite difference test for B
        function [] = B_FD_Check(this, z, theta)
            [~, grad] = this.J(z, theta);
            h = 10.^(-1:-1:-6);
            v = randn(length(theta), 1);
            B = this.Compute_B(z, theta);
            Bv = B * v;
            error = zeros(length(h), 1);
            for k = 1:length(h)
                [~, gradk] = this.J(z, theta + h(k) * v);
                error(k) = norm((gradk - grad) / h(k) - Bv) / norm(Bv);
            end

            disp('B finite difference test');
            for k = 1:length(h)
                disp(['Step size = ', num2str(h(k)), ' and error = ', num2str(error(k))]);
            end
        end

    end
end
