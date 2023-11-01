classdef HDSA_MD_Interface_synthetic_test_with_hessian_gevp < HDSA_MD_Interface

    properties
        m  % Mesh resolution
        x  % Mesh nodes on [0,1]
        S  % Stiffness matrix
        M  % Mass matrix
        W_d  % Discrepancy precision matrix
        W_u  % State weighting matrix
        W_z  % Control weighting matrix
    end

    methods (Access = public)

        function [u_out] = Apply_W_d(this, u_in)
            u_out = this.W_d * u_in;
        end

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

        function [u_opt] = Load_Optimal_u(this)
            u_opt = (1 + this.x).^3;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = 1 + this.x;
        end

        function [Z] = Load_Z_Data(this)
            Z = zeros(this.m, 2);
            Z(:, 1) = 1 + this.x;
            Z(:, 2) = this.x + this.x.^2;
        end

        function [D] = Load_d_Data(this)
            Z = this.Load_Z_Data();
            D = .2 * (Z.^3);
        end

    end

    methods

        function this = HDSA_MD_Interface_synthetic_test_with_hessian_gevp(m)
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
            E_z = (1.e2) * ((1.e-2) * S + M);

            this.W_d = E_d' * linsolve(M, E_d);
            this.W_u = E_u' * linsolve(M, E_u);
            this.W_z = E_z' * linsolve(M, E_z);
        end

    end

end
