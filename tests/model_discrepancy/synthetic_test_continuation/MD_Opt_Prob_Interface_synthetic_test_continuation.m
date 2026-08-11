%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Opt_Prob_Interface_synthetic_test_continuation < MD_Opt_Prob_Interface

    properties
        m
        x
        M
    end

    methods (Access = public)

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            z_out = 3 * diag(z.^2) * u_in;
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian(this, z_in, z)
            z_out = 3 * diag(z.^2) * z_in;
        end

        function [z_out] = Apply_Solution_Operator_z_Hessian_Adjoint(this, z_in, u_adj, z)
            z_out = 6 * z .* z_in .* u_adj;
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            u = this.State_Solve(z);
            grad_u = this.Misfit_Gradient(u, z);
            Sz_zin = 3 * z.^2 .* z_in;
            gauss_newton_term = 3 * z.^2 .* (this.M * Sz_zin);
            second_order_term = 6 * z .* z_in .* grad_u;
            z_out = gauss_newton_term + second_order_term;
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            grad_u = this.M * (u - (1 + this.x).^3);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = this.M * u_in;
        end

        function [u] = State_Solve(this, z)
            u = z.^3;
        end

        function [val, grad_u, grad_z] = Objective_Function(this, u, z)
            grad_u = this.M * (u - (1 + this.x).^3);
            val = 0.5 * (u - (1 + this.x).^3)' * grad_u;
            grad_z = zeros(size(z));
        end

    end

    methods

        function this = MD_Opt_Prob_Interface_synthetic_test_continuation(m)
            this.m = m;
            this.x = linspace(0, 1, m)';

            h = this.x(2) - this.x(1);
            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;
        end

    end

end
