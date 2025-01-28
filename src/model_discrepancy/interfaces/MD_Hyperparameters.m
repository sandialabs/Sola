classdef MD_Hyperparameters < handle

    properties
        data_interface

        data_noise_percent
        discrepancy_uncertainty_percent
        alpha_u_init
        num_z_samples
        discrepancy_percent_z_variation
        W_u_inv_spectral_gap

        alpha_d
        alpha_u
        beta_u
        alpha_z
        beta_z
        gsvd_num_sing_vals
        gsvd_oversampling
        gsvd_num_subspace_iter
    end

    methods (Access = public)

        % Defaults to return empty array
        % Overload function to load node data
        function [nodes] = Load_Node_Data(this)
            nodes = [];
            disp('Load_Node_Data must be implemented to use Determine_alpha_u and/or Determine_beta_u');
        end

    end


    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_alpha_d(this,alpha_d_new)
            this.alpha_d = alpha_d_new;
        end

        function [] = Determine_alpha_d(this)
            alpha_d_new = (this.data_noise_percent * mean(abs(this.data_interface.D(:))))^2;
            this.Set_alpha_d(alpha_d_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_alpha_u(this,alpha_u_new)
            this.alpha_u = alpha_u_new;
        end

        function [] = Determine_alpha_u(this, u_prior_interface)

            delta = this.data_interface.Load_d_Data();
            delta_norm = sqrt(delta(:,1)' * u_prior_interface.Apply_M_u(delta(:,1)));
            delta_norm = this.discrepancy_uncertainty_percent * mean(diag(delta_norm));

            alpha_u_new = delta_norm/sum(u_prior_interface.sing_vals.^2);
            this.Set_alpha_u(alpha_u_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_beta_u(this,beta_u_new)
            this.beta_u = beta_u_new;
        end

        function [] = Determine_beta_u(this)
            nodes = this.Load_Node_Data();
            if size(nodes,2) == 1
                correlation_length = computeCorrelationLength_1D(nodes(:,1),this.data_interface.D(:,1));
                beta_u_new = correlation_length^2/18;
            elseif size(nodes,2) == 2
                correlation_length = computeCorrelationLength_2D(nodes(:,1),nodes(:,2),this.data_interface.D(:,1));
                beta_u_new = correlation_length^2/8;
            end
            this.Set_beta_u(beta_u_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_alpha_z(this,alpha_z_new)
            this.alpha_z = alpha_z_new;
        end

        function [] = Determine_alpha_z(this, z_prior_interface)
            z_prior_interface.Set_alpha_z(1.0);
            samples = z_prior_interface.Sample_with_Covariance_W_z_Inverse(this.num_z_samples);
            tmp = z_prior_interface.Apply_M_z(this.data_interface.z_opt);
            tmp = sqrt(tmp'*this.data_interface.z_opt);
            for k = 1:this.num_z_samples
                samples(:,k) = tmp*samples(:,k)/sqrt(samples(:,k)'*z_prior_interface.Apply_M_z(samples(:,k)));
            end
            sample_norms = zeros(this.num_z_samples,1);
            for k = 1:this.num_z_samples
                E_z = z_prior_interface.Apply_E_z_Inverse(samples(:,k));
                sample_norms(k) = E_z'*z_prior_interface.Apply_M_z(E_z);
            end
            alpha_z_new = (this.discrepancy_percent_z_variation.^2) / mean(sample_norms);
            this.Set_alpha_z(alpha_z_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_beta_z(this,beta_z_new)
            this.beta_z = beta_z_new;
        end

        function [] = Determine_beta_z(this)
            nodes = this.Load_Node_Data();
            if size(nodes,2) == 1
                correlation_length = computeCorrelationLength_1D(nodes(:,1),this.data_interface.Z(:,1));
                beta_z_new = correlation_length^2/18;
            elseif size(nodes,2) == 2
                correlation_length = computeCorrelationLength_2D(nodes(:,1),nodes(:,2),this.data_interface.Z(:,1));
                beta_z_new = correlation_length^2/8;
            end
            this.Set_beta_z(beta_z_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_GSVD_Hyperparameters(this,gsvd_num_sing_vals_new,gsvd_oversampling_new,gsvd_num_subspace_iter_new)
            this.gsvd_num_sing_vals = gsvd_num_sing_vals_new;
            this.gsvd_oversampling = gsvd_oversampling_new;
            this.gsvd_num_subspace_iter = gsvd_num_subspace_iter_new;
        end

        function [] = Determine_GSVD_Hyperparameters(this,m)
            num_sing_vals = round(sqrt((1/(4*pi^2*this.beta_u)) * (-1 + (1/sqrt(this.W_u_inv_spectral_gap))*(4*pi^2+this.beta_u+1))));
            this.gsvd_num_sing_vals = min(num_sing_vals,m);
            this.gsvd_oversampling = min(10,m-this.gsvd_num_sing_vals);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function this = MD_Hyperparameters(data_interface)
            this.data_interface = data_interface;

            this.data_noise_percent = 0.001;
            this.num_z_samples = 50;
            this.discrepancy_percent_z_variation = 0.2;
            this.W_u_inv_spectral_gap = 1.e-6;

            this.alpha_d = 0.0;
            this.alpha_u = 0.0;
            this.beta_u = 0.0;
            this.alpha_z = 0.0;
            this.beta_z = 0.0;
            this.gsvd_num_sing_vals = 0;
            this.gsvd_oversampling = 0;
            this.gsvd_num_subspace_iter = 1;
        end
    end
end