classdef MD_u_Prior_Interface_synthetic_test_with_hessian_gevp < MD_u_Prior_Interface

    properties
        m  % Mesh resolution
        x  % Mesh nodes on [0,1]
        S  % Stiffness matrix
        M  % Mass matrix
        W_d  % Discrepancy precision matrix
        W_u  % State weighting matrix
    end

    methods (Access = public)

        function [u_out] = Apply_W_d(this, u_in)
            u_out = this.W_d * u_in;
        end

        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse(this, u_in, scalar)
            u_out = linsolve(this.W_u + scalar * this.W_d, u_in);
        end

        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse_Factor(this, u_in, scalar)
            R = chol(this.W_u + scalar * this.W_d);
            u_out = linsolve(R, u_in);
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = linsolve(this.W_u, u_in);
        end

        function [u_out] = Apply_W_u_Inverse_Factor(this, u_in)
            R = chol(this.W_u);
            u_out = linsolve(R, u_in);
        end

    end

    methods

        function this = MD_u_Prior_Interface_synthetic_test_with_hessian_gevp(m)
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

            E_d = (1.e-6) * S + M;
            E_u = (2.0) * ((5.e-2) * S + M);

            this.W_d = E_d' * linsolve(M, E_d);
            this.W_u = E_u' * linsolve(M, E_u);
        end

    end

end
