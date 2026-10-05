%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Subsurface_HiFi_Constraint < Constraint
    % High-fidelity pressure equation on Omega = (0,1):
    %
    %   -d/dx( (k_HF(p) / mu) * p_x ) + ell(x) * (p-p0)^3 = q(x),
    %   p(0) = p(1) = p0,
    %
    % where k_HF(p) = k0 * exp(alpha*(p-p0)) and ell(x) is a localized
    % leakoff / fracture-storage coefficient.  The low-fidelity model omits
    % this pressure-triggered sink, which creates a sizeable nonlinear model
    % discrepancy while keeping the Jacobian coercive and well conditioned.

    properties
        lofi
        m
        k0
        alpha
        p0
        viscosity
        leakoff_coeff
        leakoff_center
        leakoff_width
        x
        M
        S
        interior
        coll_weights
        coll_points
        nodes_to_coll_points
        nodes_to_coll_grads
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            % Seed the high-fidelity nonlinear solve with the low-fidelity
            % state whenever possible.  This is typically much closer than a
            % constant-permeability linear solve for the optimized controls.
            try
                u0 = this.lofi.State_Solve(z);
            catch
                A = (this.k0 / this.viscosity) * this.S;
                b = this.M * z;
                A(1, :) = 0; A(1, 1) = 1; b(1) = this.p0;
                A(end, :) = 0; A(end, end) = 1; b(end) = this.p0;
                u0 = linsolve(A, b);
            end

            options = optimoptions('fsolve', 'Display', 'none', ...
                'OptimalityTolerance', 1.e-14, 'FunctionTolerance', 1.e-14, ...
                'SpecifyObjectiveGradient', true);
            u = fsolve(@(u)this.Constraint_Evaluation(u, z), u0, options);
        end

        function [c, c_u, c_z] = Constraint_Evaluation(this, u, z)
            [K, K_u] = this.Assemble_Diffusion_Function(u);
            [L, L_u] = this.Assemble_Leakoff_Function(u);
            c = K + L - this.M * z;
            c_u = K_u + L_u;
            c_z = -this.M;

            c([1, this.m]) = u([1, this.m]) - this.p0;
            c_u([1, this.m], :) = 0;
            c_u(1, 1) = 1; c_u(this.m, this.m) = 1;
            c_z([1, this.m], :) = 0;
        end

        function [con] = c(this, u, z)
            con = this.Constraint_Evaluation(u, z);
        end

        function [c_uu] = Constraint_Hessian(this, u, lambda)
            c_uu = this.Assemble_Diffusion_Function_Hessian(u, lambda) ...
                 + this.Assemble_Leakoff_Function_Hessian(u, lambda);
            c_uu([1, this.m], :) = 0;
            c_uu(:, [1, this.m]) = 0;
        end

        function [K, K_u] = Assemble_Diffusion_Function(this, u)
            N = this.nodes_to_coll_points;
            G = this.nodes_to_coll_grads;
            p = N * u;
            g = G * u;
            [a, a_prime] = this.Permeability_Function(p);

            K = G' * (this.coll_weights .* a .* g);
            K_u = G' * diag(this.coll_weights .* a) * G ...
                + G' * diag(this.coll_weights .* a_prime .* g) * N;
        end

        function [K_uu] = Assemble_Diffusion_Function_Hessian(this, u, lambda)
            N = this.nodes_to_coll_points;
            G = this.nodes_to_coll_grads;
            p = N * u;
            g = G * u;
            g_lambda = G * lambda;
            [~, a_prime, a_prime_prime] = this.Permeability_Function(p);

            K_uu = G' * diag(this.coll_weights .* a_prime .* g_lambda) * N ...
                 + N' * diag(this.coll_weights .* a_prime_prime .* g .* g_lambda) * N ...
                 + N' * diag(this.coll_weights .* a_prime .* g_lambda) * G;
        end

        function [L, L_u] = Assemble_Leakoff_Function(this, u)
            N = this.nodes_to_coll_points;
            p = N * u;
            r = p - this.p0;
            ell = this.Leakoff_Function(this.coll_points);
            L = N' * (this.coll_weights .* ell .* (r.^3));
            L_u = N' * diag(3 * this.coll_weights .* ell .* (r.^2)) * N;
        end

        function [L_uu] = Assemble_Leakoff_Function_Hessian(this, u, lambda)
            N = this.nodes_to_coll_points;
            p = N * u;
            r = p - this.p0;
            lambda_nodes = N * lambda;
            ell = this.Leakoff_Function(this.coll_points);
            L_uu = N' * diag(6 * this.coll_weights .* ell .* r .* lambda_nodes) * N;
        end

        function [a, a_prime, a_prime_prime] = Permeability_Function(this, p)
            a = (this.k0 / this.viscosity) * exp(this.alpha * (p - this.p0));
            a_prime = this.alpha * a;
            a_prime_prime = (this.alpha^2) * a;
        end

        function [ell] = Leakoff_Function(this, x)
            ell = this.leakoff_coeff * exp(-((x - this.leakoff_center) / this.leakoff_width).^2);
        end

        function [diff] = Finite_Difference_Diffusion_Function_Jacobian(this, u)
            [K, K_u] = this.Assemble_Diffusion_Function(u);
            h = 10.^(-1:-1:-6);
            v = randn(length(u), 1);
            v = v / norm(v);
            diff = zeros(6, 1);
            for k = 1:6
                [K_pert, ~] = this.Assemble_Diffusion_Function(u + h(k) * v);
                diff(k) = norm(K_u * v - (K_pert - K) / h(k)) / norm(K_u * v);
            end
            disp(log10(diff'));
        end

        function [diff] = Finite_Difference_Diffusion_Function_Hessian(this, u, lambda)
            K_uu = this.Assemble_Diffusion_Function_Hessian(u, lambda);
            [~, K_u] = this.Assemble_Diffusion_Function(u);
            h = 10.^(-1:-1:-6);
            v = randn(length(u), 1);
            v = v / norm(v);
            diff = zeros(6, 1);
            for k = 1:6
                [~, K_u_pert] = this.Assemble_Diffusion_Function(u + h(k) * v);
                diff(k) = norm(K_uu * v - (K_u_pert' * lambda - K_u' * lambda) / h(k)) / norm(K_uu * v);
            end
            disp(log10(diff'));
        end

        function [diff] = Finite_Difference_Leakoff_Function_Jacobian(this, u)
            [L, L_u] = this.Assemble_Leakoff_Function(u);
            h = 10.^(-1:-1:-6);
            v = randn(length(u), 1);
            v = v / norm(v);
            diff = zeros(6, 1);
            for k = 1:6
                [L_pert, ~] = this.Assemble_Leakoff_Function(u + h(k) * v);
                diff(k) = norm(L_u * v - (L_pert - L) / h(k)) / norm(L_u * v);
            end
            disp(log10(diff'));
        end

        function [diff] = Finite_Difference_Leakoff_Function_Hessian(this, u, lambda)
            L_uu = this.Assemble_Leakoff_Function_Hessian(u, lambda);
            [~, L_u] = this.Assemble_Leakoff_Function(u);
            h = 10.^(-1:-1:-6);
            v = randn(length(u), 1);
            v = v / norm(v);
            diff = zeros(6, 1);
            for k = 1:6
                [~, L_u_pert] = this.Assemble_Leakoff_Function(u + h(k) * v);
                diff(k) = norm(L_uu * v - (L_u_pert' * lambda - L_u' * lambda) / h(k)) / norm(L_uu * v);
            end
            disp(log10(diff'));
        end

        function [diff] = Finite_Difference_Constraint_Hessian(this, u, z, lambda)
            c_uu = this.Constraint_Hessian(u, lambda);
            [~, c_u, ~] = this.Constraint_Evaluation(u, z);
            h = 10.^(-1:-1:-6);
            v = randn(length(u), 1);
            v = v / norm(v);
            diff = zeros(6, 1);
            for k = 1:6
                [~, c_u_pert, ~] = this.Constraint_Evaluation(u + h(k) * v, z);
                diff(k) = norm(c_uu * v - (c_u_pert' * lambda - c_u' * lambda) / h(k)) / norm(c_uu * v);
            end
            disp(log10(diff'));
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            [~, c_u] = this.Constraint_Evaluation(u, z);
            Mv = linsolve(c_u', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            [~, ~, c_z] = this.Constraint_Evaluation(u, z);
            Mv = c_z' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            [~, c_u] = this.Constraint_Evaluation(u, z);
            Mv = linsolve(c_u, v);
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            [~, ~, c_z] = this.Constraint_Evaluation(u, z);
            Mv = c_z * v;
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            c_uu = this.Constraint_Hessian(u, lambda);
            Mv = c_uu * v;
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

    end

    methods (Access = public)

        function this = Subsurface_HiFi_Constraint(lofi_con, leakoff_coeff, leakoff_center, leakoff_width)
            if nargin < 2
                leakoff_coeff = 30;
            end
            if nargin < 3
                leakoff_center = 0.55;
            end
            if nargin < 4
                leakoff_width = 0.25;
            end

            this = this@Constraint();
            this.lofi = lofi_con;
            this.m = lofi_con.m;
            this.k0 = lofi_con.k0;
            this.alpha = lofi_con.alpha;
            this.p0 = lofi_con.p0;
            this.viscosity = lofi_con.viscosity;
            this.leakoff_coeff = leakoff_coeff;
            this.leakoff_center = leakoff_center;
            this.leakoff_width = leakoff_width;
            this.x = lofi_con.x;
            this.interior = lofi_con.interior;
            this.coll_weights = lofi_con.coll_weights;
            this.coll_points = lofi_con.coll_points;
            this.nodes_to_coll_points = lofi_con.nodes_to_coll_points;
            this.nodes_to_coll_grads = lofi_con.nodes_to_coll_grads;
            this.M = lofi_con.M;
            this.S = lofi_con.S;
        end

    end
end
