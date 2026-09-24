%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Heat_Problem_HiFi_Constraint < Constraint
    % High-fidelity steady heat/reaction model on Omega = (0, 1):
    %
    %   -d/dx( (1 + gamma * y^2) * y_x ) + mu * y^3 = u(x),   x in (0, 1)
    %    y(0) = y(1) = 0.
    %
    % The temperature-dependent conductivity k_gamma(y) = 1 + gamma * y^2
    % introduces a nonlinear diffusion term whose discrepancy from the
    % low-fidelity (constant-conductivity) model is controlled by gamma.
    % When gamma = 0 the high- and low-fidelity models coincide.  The
    % control is z (the distributed heat source u(x) in the outline) and
    % the state is u (the temperature y).

    properties
        lofi
        m
        mu
        gamma
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
            % Solve the linear (unit-conductivity) diffusion problem to seed
            % the nonlinear solver.
            A = this.S;
            b = this.M * z;
            A(1, :) = 0; A(1, 1) = 1; b(1) = 0;
            A(end, :) = 0; A(end, end) = 1; b(end) = 0;
            u0 = linsolve(A, b);

            % Execute nonlinear solve to determine the state
            options = optimoptions('fsolve', 'Display', 'none', 'OptimalityTolerance', 1.e-14, 'SpecifyObjectiveGradient', true);
            u = fsolve(@(u)this.Constraint_Evaluation(u, z), u0, options);
        end

        function [c, c_u, c_z] = Constraint_Evaluation(this, u, z)
            [K, K_u] = this.Assemble_Diffusion_Function(u);
            [R, R_u] = this.Assemble_Reaction_Function(u);
            c = K + this.mu * R - this.M * z;
            c_u = K_u + this.mu * R_u;
            c_z = -this.M;
            % Impose homogeneous Dirichlet boundary conditions.
            c([1, this.m]) = u([1, this.m]);
            c_u([1, this.m], :) = 0;
            c_u(1, 1) = 1; c_u(this.m, this.m) = 1;
            c_z([1, this.m], :) = 0;
        end

        function [c_uu] = Constraint_Hessian(this, u, lambda)
            c_uu = this.Assemble_Diffusion_Function_Hessian(u, lambda) ...
                 + this.mu * this.Assemble_Reaction_Function_Hessian(u, lambda);
            c_uu([1, this.m], :) = 0;
            c_uu(:, [1, this.m]) = 0;
        end

        function [K, K_u] = Assemble_Diffusion_Function(this, u)
            % Weak form of -d/dx( (1 + gamma * y^2) * y_x ):
            %   K_i = sum_c w_c * k(y_c) * y_x,c * G_{c,i},
            % with k(y) = 1 + gamma * y^2.
            N = this.nodes_to_coll_points;
            G = this.nodes_to_coll_grads;
            y = N * u;
            g = G * u;
            k = 1 + this.gamma * (y.^2);
            k_prime = 2 * this.gamma * y;
            K = G' * (this.coll_weights .* k .* g);
            K_u = G' * diag(this.coll_weights .* k) * G ...
                + G' * diag(this.coll_weights .* k_prime .* g) * N;
        end

        function [K_uu] = Assemble_Diffusion_Function_Hessian(this, u, lambda)
            N = this.nodes_to_coll_points;
            G = this.nodes_to_coll_grads;
            y = N * u;
            g = G * u;
            g_lambda = G * lambda;
            wp = this.coll_weights .* (2 * this.gamma) .* g_lambda;
            K_uu = G' * diag(wp .* y) * N ...
                 + N' * diag(wp .* g) * N ...
                 + N' * diag(wp .* y) * G;
        end

        function [R, R_u] = Assemble_Reaction_Function(this, u)
            u_nodes = this.nodes_to_coll_points * u;
            [R_nodes, R_prime_nodes] = this.Reaction_Function(u_nodes);
            R = this.nodes_to_coll_points' * (this.coll_weights .* R_nodes);
            R_u = this.nodes_to_coll_points' * (diag(this.coll_weights) * R_prime_nodes) * this.nodes_to_coll_points;
        end

        function [R_uu] = Assemble_Reaction_Function_Hessian(this, u, lambda)
            u_nodes = this.nodes_to_coll_points * u;
            lambda_nodes = this.nodes_to_coll_points * lambda;
            R_uu = this.nodes_to_coll_points' * diag(this.coll_weights) * this.Reaction_Function_Hessian(u_nodes, lambda_nodes) * this.nodes_to_coll_points;
        end

        function [R, R_prime] = Reaction_Function(this, u)
            if ~isvector(u)
                error("The PDE cannot currently handle multiple inputs of state. Please input one state at a time!");
            end
            R = u.^3;
            R_prime = 3 * diag(u.^2);
        end

        function [R_prime_prime] = Reaction_Function_Hessian(this, u, lambda)
            R_prime_prime = 6 * diag(u .* lambda);
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

        function this = Heat_Problem_HiFi_Constraint(lofi_con, gamma)
            this = this@Constraint();
            this.lofi = lofi_con;
            this.m = lofi_con.m;
            this.mu = lofi_con.mu;
            this.gamma = gamma;
            this.x = lofi_con.x;
            this.interior = lofi_con.interior;
            this.coll_weights = lofi_con.coll_weights;
            this.coll_points = lofi_con.coll_points;
            this.nodes_to_coll_points = lofi_con.nodes_to_coll_points;
            this.M = lofi_con.M;
            this.S = lofi_con.S;

            % Map nodal values to the state derivative at each collocation
            % point.  With piecewise-linear elements the derivative is
            % constant on each element, equal to (u_{k+1} - u_k) / h.
            m = this.m;
            h = this.x(2) - this.x(1);
            nodes_to_coll_grads = zeros(2 * (m - 1), m);
            for k = 1:(m - 1)
                map_to_coll = (1:2)' + 2 * (k - 1);
                nodes_to_coll_grads(map_to_coll, k) = -1 / h;
                nodes_to_coll_grads(map_to_coll, k + 1) = 1 / h;
            end
            this.nodes_to_coll_grads = nodes_to_coll_grads;
        end

    end
end