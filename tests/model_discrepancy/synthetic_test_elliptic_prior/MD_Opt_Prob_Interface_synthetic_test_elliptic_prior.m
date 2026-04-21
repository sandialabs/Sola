%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Opt_Prob_Interface_synthetic_test_elliptic_prior < MD_Opt_Prob_Interface

    properties
        m
        x
        M
    end

    methods (Access = public)

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            z_out = 3 * diag(z.^2) * u_in;
        end

        % This implementation assumes that it is evaluated at the optimal z so that
        % the adjoint=0, a more general impl a term multiplied by the adjoint variable
        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            tmp1 = 3 * diag(z.^2) * z_in;
            tmp2 = this.M * tmp1;
            z_out = 3 * diag(z.^2) * tmp2;
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            grad_u = this.M * (u - (1 + this.x).^3);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = this.M * u_in;
        end

    end

    methods

        function this = MD_Opt_Prob_Interface_synthetic_test_elliptic_prior(m)
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
