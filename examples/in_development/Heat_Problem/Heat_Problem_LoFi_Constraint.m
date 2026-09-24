%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Heat_Problem_LoFi_Constraint < Constraint
    % Low-fidelity steady heat/reaction model on Omega = (0, 1):
    %
    %   -y_xx + mu * y^3 = u(x),      x in (0, 1)
    %    y(0) = y(1) = 0.
    %
    % The temperature-dependent conductivity of the high-fidelity model is
    % replaced by a constant (unit) conductivity, while the cubic reaction
    % term mu * y^3 is retained.  The control is z (the distributed heat
    % source u(x) in the outline) and the state is u (the temperature y).

    properties
        m
        mu
        x
        M
        S
        interior
        coll_weights
        coll_points
        nodes_to_coll_points
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            % Solve the linear diffusion problem to seed the nonlinear solver
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
            [R, R_u] = this.Assemble_Reaction_Function(u);
            c = this.S * u + this.mu * R - this.M * z;
            c_u = this.S + this.mu * R_u;
            c_z = -this.M;
            % Impose homogeneous Dirichlet boundary conditions.
            c([1, this.m]) = u([1, this.m]);
            c_u([1, this.m], :) = 0;
            c_u(1, 1) = 1; c_u(this.m, this.m) = 1;
            c_z([1, this.m], :) = 0;
        end

        function [c_uu] = Constraint_Hessian(this, u, lambda)
            c_uu = this.mu * this.Assemble_Reaction_Function_Hessian(u, lambda);
            c_uu([1, this.m], :) = 0;
            c_uu(:, [1, this.m]) = 0;
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

        function [diff] = Finite_Difference_Reaction_Function_Jacobian(this, u)
            [R, R_u] = this.Assemble_Reaction_Function(u);
            h = 10.^(-1:-1:-6);
            v = randn(length(u), 1);
            v = v / norm(v);
            diff = zeros(6, 1);
            for k = 1:6
                [R_pert, ~] = this.Assemble_Reaction_Function(u + h(k) * v);
                diff(k) = norm(R_u * v - (R_pert - R) / h(k)) / norm(R_u * v);
            end
            disp(log10(diff'));
        end

        function [diff] = Finite_Difference_Reaction_Function_Hessian(this, u, lambda)
            R_uu = this.Assemble_Reaction_Function_Hessian(u, lambda);
            [~, R_u] = this.Assemble_Reaction_Function(u);
            h = 10.^(-1:-1:-6);
            v = randn(length(u), 1);
            v = v / norm(v);
            diff = zeros(6, 1);
            for k = 1:6
                [~, R_u_pert] = this.Assemble_Reaction_Function(u + h(k) * v);
                diff(k) = norm(R_uu * v - (R_u_pert' * lambda - R_u' * lambda) / h(k)) / norm(R_uu * v);
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

        function this = Heat_Problem_LoFi_Constraint(m, mu)
            this = this@Constraint();
            this.m = m;
            this.mu = mu;
            this.x = linspace(0, 1, m)';

            h = this.x(2) - this.x(1);
            this.interior = 2:(m - 1);

            coll_points = zeros(2 * (m - 1), 1);
            for k = 1:(m - 1)
                map_to_coll = (1:2)' + 2 * (k - 1);
                coll_points(map_to_coll) = this.x(k) + h * ((1 / sqrt(3)) * [-1; 1] + 1) / 2;
            end
            coll_weights = (h / 2) * ones(2 * (m - 1), 1);

            nodes_to_coll_points = zeros(2 * (m - 1), m);
            for k = 1:(m - 1)
                map_to_coll = (1:2)' + 2 * (k - 1);
                nodes_to_coll_points(map_to_coll(1), k) = (coll_points(map_to_coll(1)) - this.x(k + 1)) / (this.x(k) - this.x(k + 1));
                nodes_to_coll_points(map_to_coll(1), k + 1) = (coll_points(map_to_coll(1)) - this.x(k)) / (this.x(k + 1) - this.x(k));
                nodes_to_coll_points(map_to_coll(2), k) = (coll_points(map_to_coll(2)) - this.x(k + 1)) / (this.x(k) - this.x(k + 1));
                nodes_to_coll_points(map_to_coll(2), k + 1) = (coll_points(map_to_coll(2)) - this.x(k)) / (this.x(k + 1) - this.x(k));
            end

            this.coll_weights = coll_weights;
            this.coll_points = coll_points;
            this.nodes_to_coll_points = nodes_to_coll_points;

            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            S = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;
            this.S = S;
        end

    end
end
