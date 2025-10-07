classdef MD_Lumped_Mass_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    properties
        M
        M_lumped_diag
        S
        E_u
        W_u_acute


        u_hyperparam_interface
        determine_u_hyperparams
        beta_u
        n_u
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            A = this.W_u_acute + scalar * this.M;
            u_out = A \ u_in;
        end


        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)
            u_out = this.W_u_acute \ u_in;
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Inverse(this, num_samples)
            omega = randn(this.n_u,num_samples);
            u_out = this.E_u \ (diag(sqrt(this.M_lumped_diag)) * omega);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            omega = randn(this.n_u,num_samples);
            A = this.W_u_acute + scalar * this.M;
            R = chol(A);
            u_out = R \ omega;
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Lumped_Mass_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface)
            arguments
                S (:, :) {mustBeNumeric}
                M (:, :) {mustBeNumeric}
                data_interface MD_Data_Interface
                u_hyperparam_interface MD_u_Hyperparameter_Interface
            end
            this@MD_Scaled_u_Prior_Interface(u_hyperparam_interface.alpha_u);

            this.M = M;
            this.S = S;
            this.u_hyperparam_interface = u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface, u_hyperparam_interface);

            if this.u_hyperparam_interface.beta_u == 0.0
                this.determine_u_hyperparams.Determine_beta_u();
            end
            this.Set_beta_u(this.u_hyperparam_interface.beta_u);

            this.Assemble_Operators();

            if ~this.u_hyperparam_interface.is_transient
                if this.u_hyperparam_interface.alpha_u == 0.0
                    this.determine_u_hyperparams.Determine_alpha_u(this);
                end
                this.Set_alpha_u(this.u_hyperparam_interface.alpha_u);
            end
            this.n_u = size(this.M,1);
        end

        function [] = Set_beta_u(this, beta_u_new)
            this.beta_u = beta_u_new;
        end

        function [] = Assemble_Operators(this)
            this.M_lumped_diag = sum(this.M,2);
            this.E_u = this.beta_u * this.S + this.M;
            this.W_u_acute = this.E_u' * diag(1./this.M_lumped_diag) * this.E_u;
        end

    end

end
