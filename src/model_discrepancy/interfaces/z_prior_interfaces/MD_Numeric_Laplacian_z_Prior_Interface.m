classdef MD_Numeric_Laplacian_z_Prior_Interface < MD_Elliptic_z_Prior_Interface

    properties
        beta_z
        S
        M
        hyperparams
        E_z
        M_sqrt
    end

    methods ( Access = public)

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = this.E_z \ z_in;
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = this.E_z' \ z_in;
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M * z_in;
        end

        function [z_out] = Apply_E_z(this, z_in)
            z_out = this.E_z * z_in;
        end

        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = this.E_z' * z_in;
        end

        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = this.M_z' \ z_in;
        end

        function [z_out] = Sample_with_Covariance_W_z_Inverse(this, num_samples)
            omega = randn(size(this.S,1),num_samples);
            vec = this.M_sqrt.Matrix_Apply(omega);
            z_out = sqrt(this.alpha_z) * this.Apply_E_z_Inverse(vec);
        end

        function [] = Set_beta_z(this,beta_z_new)
            this.beta_z = beta_z_new;
            this.E_z = this.beta_z * this.S + this.M;
        end

        function this = MD_Numeric_Laplacian_z_Prior_Interface(S,M,hyperparams)
            this@MD_Elliptic_z_Prior_Interface(hyperparams.alpha_z)
            this.S = S;
            this.M = M;
            this.hyperparams = hyperparams;
            this.M_sqrt = M_z_Sqrt(this);

            if hyperparams.beta_z == 0.0
                this.hyperparams.Determine_beta_z();
            end
            this.Set_beta_z(this.hyperparams.beta_z);

            if this.hyperparams.alpha_z == 0.0
                this.hyperparams.Determine_alpha_z(this);
            end
            this.Set_alpha_z(this.hyperparams.alpha_z);
        end

    end

end
