classdef MD_z_Prior_Interface_synthetic_test_with_hessian_gevp < MD_z_Prior_Interface

    properties
        m  % Mesh resolution
        x  % Mesh nodes on [0,1]
        S  % Stiffness matrix
        M  % Mass matrix
        W_z  % Control weighting matrix
    end

    methods (Access = public)

        function [z_out] = Apply_W_z_Inverse(this, z_in)
            z_out = linsolve(this.W_z, z_in);
        end

        function [z_out] = Apply_W_z_Inverse_Factor(this, z_in)
            R = chol(this.W_z);
            z_out = linsolve(R, z_in);
        end

        function [z_out] = Apply_W_z(this, z_in)
            z_out = this.W_z * z_in;
        end

    end

    methods

        function this = MD_z_Prior_Interface_synthetic_test_with_hessian_gevp(m)
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

            E_z = (1.e2) * ((1.e-2) * S + M);
            this.W_z = E_z' * linsolve(M, E_z);
        end

    end

end
