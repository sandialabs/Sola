%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Numeric_Laplacian_z_Prior_Interface < MD_Elliptic_z_Prior_Interface

    properties
        beta_z
        S
        M
        z_hyperparam_interface
        determine_z_hyperparams
        E_z
        M_sqrt
    end

    %% Implementation of base class virtual functions
    methods (Access = public)

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
            z_out = this.M' \ z_in;
        end

        function [z_out] = Sample_with_Covariance_W_z_Acute_Inverse(this, num_samples)
            omega = randn(size(this.S, 1), num_samples);
            vec = this.M_sqrt.Matrix_Sqrt_Apply(omega);
            z_out = this.Apply_E_z_Inverse(vec);
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface)
            arguments
                S (:, :) {mustBeNumeric}
                M (:, :) {mustBeNumeric}
                data_interface MD_Data_Interface
                z_hyperparam_interface MD_z_Hyperparameter_Interface
                u_prior_interface MD_u_Prior_Interface
            end
            this@MD_Elliptic_z_Prior_Interface(z_hyperparam_interface.alpha_z);
            this.S = S;
            this.M = M;
            this.z_hyperparam_interface = z_hyperparam_interface;
            this.determine_z_hyperparams = MD_Determine_z_Hyperparameters(data_interface, z_hyperparam_interface, u_prior_interface);
            this.M_sqrt = M_z_Sqrt(this);

            if this.z_hyperparam_interface.beta_z == 0.0
                this.determine_z_hyperparams.Determine_beta_z();
            end
            this.Set_beta_z(this.z_hyperparam_interface.beta_z);

            if this.z_hyperparam_interface.alpha_z == 0.0
                this.determine_z_hyperparams.Determine_alpha_z(this);
            end
            this.Set_alpha_z(this.z_hyperparam_interface.alpha_z);
        end

        function [] = Set_beta_z(this, beta_z_new)
            this.beta_z = beta_z_new;
            this.E_z = this.beta_z * this.S + this.M;
        end

    end

end
