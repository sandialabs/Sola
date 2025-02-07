classdef MD_Hyperparameters < handle

    properties
        data_interface

        data_noise_percent
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
            delta_norm = delta(:,1)' * u_prior_interface.Apply_M_u(delta(:,1));

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
                correlation_lengths = zeros(size(this.data_interface.D,2),1);
                for k = 1:length(correlation_lengths)
                    correlation_lengths(k) = computeCorrelationLength_1D(nodes(:,1),this.data_interface.D(:,k));
                end
                beta_u_new = mean(correlation_lengths)^2/12;
            elseif size(nodes,2) == 2
                correlation_lengths = zeros(size(this.data_interface.D,2),1);
                for k = 1:length(correlation_lengths)
                    correlation_lengths(k) = computeCorrelationLength_2D(nodes(:,1),nodes(:,2),this.data_interface.D(:,k));
                end
                beta_u_new = mean(correlation_lengths)^2/8;
            else
                disp('Determine_beta_u error: Dimensions greater than 2 are not supported.')
            end
            this.Set_beta_u(beta_u_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_alpha_z(this,alpha_z_new)
            this.alpha_z = alpha_z_new;
        end

        function [] = Determine_alpha_z(this, z_prior_interface)
            tmp = z_prior_interface.Apply_M_z(this.data_interface.z_opt);
            zopt_norm = tmp'*this.data_interface.z_opt;

            nodes = this.Load_Node_Data();
            if size(nodes,2) == 1
                Lx = max(nodes(:,1))-min(nodes(:,1));
                n = length(nodes(:,1))-1;
                e = 1 + z_prior_interface.beta_z * (pi/Lx)^2 * (0:n).^2;
                e = e';
            elseif size(nodes,2) == 2
                Lx = max(nodes(:,1))-min(nodes(:,1));
                Ly = max(nodes(:,2))-min(nodes(:,2));
                n = sqrt(length(nodes(:,1)))-1;
                e = 1 + z_prior_interface.beta_z * pi^2 * ( kron(((0:n).^2)',ones(n+1,1))/Lx^2 + kron(ones(n+1,1),((0:n).^2)')/Ly^2 );
            else
                disp('Determine_alpha_z error: Dimensions greater than 2 are not supported.')
            end

            evals = sort(1./e,'descend');
            I = find(evals < 1.e-2);
            rank = n;
            if ~isempty(I)
                rank = I(1);
            end
            evals = evals(1:rank);
            samples = 1000;
            nu = randn(samples,rank).^2;
            tmp = mean(( nu * evals.^4 ) ./ ( nu * evals.^2 ));
            
            alpha_z_new = (this.discrepancy_percent_z_variation.^2) / (tmp * zopt_norm);
            this.Set_alpha_z(alpha_z_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_beta_z(this,beta_z_new)
            this.beta_z = beta_z_new;
        end

        function [] = Determine_beta_z(this)
            nodes = this.Load_Node_Data();
            if size(nodes,2) == 1
                correlation_lengths = zeros(size(this.data_interface.Z,2),1);
                for k = 1:length(correlation_lengths)
                    correlation_lengths(k) = computeCorrelationLength_1D(nodes(:,1),this.data_interface.Z(:,k));
                end
                beta_z_new = mean(correlation_lengths)^2/12;
            elseif size(nodes,2) == 2
                correlation_lengths = zeros(size(this.data_interface.Z,2),1);
                for k = 1:length(correlation_lengths)
                    correlation_lengths(k) = computeCorrelationLength_2D(nodes(:,1),nodes(:,2),this.data_interface.Z(:,k));
                end
                beta_z_new = mean(correlation_lengths)^2/8;
            else
                disp('Determine_beta_z error: Dimensions greater than 2 are not supported.')
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
            this.discrepancy_percent_z_variation = 1.0;
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