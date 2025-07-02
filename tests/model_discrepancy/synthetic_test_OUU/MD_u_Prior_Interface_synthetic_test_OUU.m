classdef MD_u_Prior_Interface_synthetic_test_OUU < MD_u_Prior_Interface

    properties
        m  % Mesh resolution
        x  % Mesh nodes on [0,1]
        S  % Stiffness matrix
        M  % Mass matrix
        W_u  % State weighting matrix
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = linsolve(this.W_u + scalar * this.M, u_in);
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = linsolve(this.W_u, u_in);
        end

        % Compute samples from a mean zero Gaussian with covariance W_u^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            R = chol(this.W_u);
            u_out = linsolve(R, randn(size(R, 1), num_samples));
        end

        % Compute samples from a mean zero Gaussian with covariance (W_u+scalar*M_u)^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            R = chol(this.W_u + scalar * this.M);
            u_out = linsolve(R, randn(size(R, 1), num_samples));
        end

    end

    methods

        function this = MD_u_Prior_Interface_synthetic_test_OUU(m)
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

            E_u = (1.0) * ((5.e-2) * S + M);

            this.W_u = E_u' * linsolve(M, E_u);
        end

    end

end
