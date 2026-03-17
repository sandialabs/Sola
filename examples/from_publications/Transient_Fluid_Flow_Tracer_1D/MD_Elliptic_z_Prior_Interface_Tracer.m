classdef MD_Elliptic_z_Prior_Interface_Tracer < MD_Elliptic_z_Prior_Interface

    properties
        E_z
        M
    end

    methods

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(this.E_z, z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(this.E_z', z_in);
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z(this, z_in)
            z_out = this.E_z * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = this.E_z' * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = linsolve(this.M, z_in);
        end

        % Compute samples from a mean zero Gaussian with covariance W_z^{-1}
        function [z_out] = Sample_with_Covariance_W_z_Inverse(this, num_samples)
            R = chol(this.M);
            z_out = sqrt(this.alpha_z) * linsolve(this.E_z, R' * randn(size(R, 1), num_samples));
        end

        function this = MD_Elliptic_z_Prior_Interface_Tracer(alpha_u, sabl_opt)
            this@MD_Elliptic_z_Prior_Interface(alpha_u);

            S = sabl_opt.con.S_z;
            this.M = sabl_opt.con.M_z;
            this.E_z = (3.e-2) * S + this.M;
        end

    end

end
