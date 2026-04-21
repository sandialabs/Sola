%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Elliptic_u_Prior_Interface_synthetic_test_elliptic_prior < MD_Elliptic_u_Prior_Interface

    properties
        m  % Mesh resolution
        x  % Mesh nodes on [0,1]
        S  % Stiffness matrix
        M  % Mass matrix
        E_u  % State precision matrix elliptic operator
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

    end

    methods

        function this = MD_Elliptic_u_Prior_Interface_synthetic_test_elliptic_prior(alpha_u, m)
            this@MD_Elliptic_u_Prior_Interface(alpha_u);

            this.m = m;
            this.x = linspace(0, 1, m)';

            h = this.x(2) - this.x(1);
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

            this.E_u = (5.e-2) * S + M;

            num_sing_vals = 50;
            oversampling = 1;
            num_subspace_iters = 2;
            u_vec = zeros(m, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

        end

    end

end
