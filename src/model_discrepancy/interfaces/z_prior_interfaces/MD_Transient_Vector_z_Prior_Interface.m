%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Transient_Vector_z_Prior_Interface < MD_Elliptic_z_Prior_Interface

    properties
        beta_t
        S
        M
        n_t
        num_controls
        z_hyperparam_interface
        determine_z_hyperparams
        E_t
        V
        lambda
    end

    %% Implementation of base class virtual functions
    methods (Access = public)

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = kron(this.V * diag(1 ./ sqrt(this.lambda)) * this.V', eye(this.num_controls)) * z_in;
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = kron(this.V * diag(1 ./ sqrt(this.lambda)) * this.V', eye(this.num_controls)) * z_in;
        end

        function [z_out] = Apply_E_z(this, z_in)
            z_out = kron(this.V * diag(sqrt(this.lambda)) * this.V', eye(this.num_controls)) * z_in;
        end

        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = kron(this.V * diag(sqrt(this.lambda)) * this.V', eye(this.num_controls)) * z_in;
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = kron(this.M, eye(this.num_controls)) * z_in;
        end

        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = linsolve(kron(this.M, eye(this.num_controls)), z_in);
        end

        function [z_out] = Sample_with_Covariance_W_z_Acute_Inverse(this, num_samples)
            z_out = zeros(this.n_t * this.num_controls, num_samples);

            % This piece of code is more efficient if num_samples > num_controls
            for k = 1:this.num_controls
                I = k:this.num_controls:(this.n_t * this.num_controls);
                z_out(I, :) = this.V * diag(sqrt(this.lambda)) * randn(this.n_t, num_samples);
            end

            % This piece of code replicates the HdsaLib implementation
            % for k = 1:num_samples
            %      omega = randn(this.num_controls,this.n_t)';
            %      tmp = this.V * diag(sqrt(this.lambda)) * omega;
            %      tmp = tmp';
            %      z_out(:,k) = tmp(:);
            %  end
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Transient_Vector_z_Prior_Interface(S, M, num_controls, data_interface, z_hyperparam_interface, u_prior_interface)
            arguments
                S (:, :) {mustBeNumeric}
                M (:, :) {mustBeNumeric}
                num_controls (1, 1) {mustBeNumeric}
                data_interface MD_Data_Interface
                z_hyperparam_interface MD_z_Hyperparameter_Interface
                u_prior_interface MD_u_Prior_Interface
            end
            this@MD_Elliptic_z_Prior_Interface(z_hyperparam_interface.alpha_z);
            this.S = S;
            this.M = M;
            this.n_t = size(M, 1);
            this.num_controls = num_controls;
            this.z_hyperparam_interface = z_hyperparam_interface;
            this.determine_z_hyperparams = MD_Determine_z_Hyperparameters(data_interface, z_hyperparam_interface, u_prior_interface);

            if z_hyperparam_interface.beta_t == 0.0
                this.determine_z_hyperparams.Determine_beta_t();
            end
            this.Set_beta_t(this.z_hyperparam_interface.beta_t);

            if this.z_hyperparam_interface.alpha_z == 0.0
                this.determine_z_hyperparams.Determine_alpha_z(this);
            end
            this.Set_alpha_z(this.z_hyperparam_interface.alpha_z);
        end

        function [] = Set_beta_t(this, beta_t_new)
            this.beta_t = beta_t_new;
            this.E_t = this.beta_t * this.S + this.M;
            [this.V, this.lambda] = eig(this.E_t, this.M, 'vector');
            this.V = this.V(:, end:-1:1);
            this.lambda = this.lambda(end:-1:1);
        end

    end

end
