classdef MD_Elliptic_z_Prior_Interface_Mass_Spring < MD_Elliptic_z_Prior_Interface

    properties
        sabl_opt
        E_z
        P_z
        M
    end

    methods

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(this.P_z' * this.E_z * this.P_z, z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(this.P_z' * this.E_z' * this.P_z, z_in);
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.P_z' * this.M * this.P_z * z_in;
        end

        % Compute samples from a mean zero Gaussian with covariance W_z^{-1}
        function [z_out] = Sample_with_Covariance_W_z_Inverse(this, num_samples)
            R = chol(this.P_z' * this.M * this.P_z);
            z_out = sqrt(this.alpha_z) * linsolve(this.P_z' * this.E_z * this.P_z, R' * randn(size(R, 1), num_samples));
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z(this, z_in)
            z_out = this.P_z' * this.E_z * this.P_z * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = this.P_z' * this.E_z' * this.P_z * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = linsolve(this.P_z' * this.M * this.P_z, z_in);
        end

        function this = MD_Elliptic_z_Prior_Interface_Mass_Spring(alpha_z, sabl_opt)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);

            this.sabl_opt = sabl_opt;

            n_t = sabl_opt.con.n_t;
            h = sabl_opt.con.t_mesh(2) - sabl_opt.con.t_mesh(1);
            M = diag(4 * ones(1, n_t)) + diag(ones(1, n_t - 1), 1) + diag(ones(1, n_t - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            S = diag(2 * ones(1, n_t)) + (-1) * diag(ones(1, n_t - 1), 1) + (-1) * diag(ones(1, n_t - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;

            this.P_z = sabl_opt.con.mass_spring_coupled.P_z;
            this.E_z = (1.e-1) * S + this.M;
            this.E_z(:,1) = 0;
            this.E_z(1,:) = 0;
            this.E_z(1,1) = 50;
        end

    end

end
