%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Elliptic_z_Prior_Interface_Diff_React < MD_Elliptic_z_Prior_Interface

    properties
        E_z
        M
    end

    methods

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = this.E_z \ z_in;
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = this.E_z' \ z_in;
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
            z_out = this.M \ z_in;
        end

        % Compute samples from a mean zero Gaussian with covariance W_z^{-1}
        function [z_out] = Sample_with_Covariance_W_z_Inverse(this, num_samples)
            mass_mat_sqrt = M_z_Sqrt(this);
            z_in = randn(size(this.E_z, 1), num_samples);
            tmp = mass_mat_sqrt.Matrix_Sqrt_Apply(z_in);
            z_out = sqrt(this.alpha_z) * this.Apply_E_z_Inverse(tmp);
        end

        function this = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, sola_opt)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);

            S = sola_opt.con.diff_react_lofi.pde_meshing.S;
            this.M = sola_opt.con.diff_react_lofi.pde_meshing.M;
            this.E_z = (1.e-1) * S + this.M;
        end

    end

end
