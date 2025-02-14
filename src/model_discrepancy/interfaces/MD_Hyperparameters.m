classdef MD_Hyperparameters < handle

    properties
        data_interface
        is_transient

        data_noise_percent
        discrepancy_percent_z_variation
        W_u_inv_spectral_gap

        alpha_d
        alpha_u
        beta_u
        alpha_z
        beta_z
        alpha_t
        beta_t

        gsvd_num_sing_vals
        gsvd_oversampling
        gsvd_num_subspace_iter

        d1_norm_sq
        d_pert_norm_sq
        z_pert_norm_sq
    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%% Functions to be overloaded to enable some functionality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = [];
            disp('Load_Spatial_Node_Data is required for automate hyperparameters')
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = [];
            disp('Load_Time_Node_Data is required for automate hyperparameters for transient problems')
        end

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
            this.d1_norm_sq = delta(:,1)' * u_prior_interface.Apply_M_u(delta(:,1));

            N = size(this.data_interface.D,2);
            if N > 1
                this.d_pert_norm_sq = zeros(N-1,1);
                for k = 2:N
                    v = this.data_interface.D(:,k) - this.data_interface.D(:,1);
                    this.d_pert_norm_sq(k-1) = v'*u_prior_interface.Apply_M_u(v);
                end
            end

            if this.is_transient
                u_op_trace = sum(u_prior_interface.spatial_prior_cov.sing_vals.^2) * sum(u_prior_interface.transient_prior_cov.evals);
            else
                u_op_trace = sum(u_prior_interface.sing_vals.^2);
            end

            alpha_u_new = this.d1_norm_sq/u_op_trace;
            this.Set_alpha_u(alpha_u_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_beta_u(this,beta_u_new)
            this.beta_u = beta_u_new;
        end

        function [] = Determine_beta_u(this)
            nodes = this.Load_Spatial_Node_Data();
            n_y = size(nodes,1);
            n_t = size(this.data_interface.D(:,1),1)/n_y;
            N = size(this.data_interface.D,2);
            if size(nodes,2) == 1
                correlation_lengths = zeros(N,n_t);
                for i = 1:N
                    di = reshape(this.data_interface.D(:,i),n_y,n_t);
                    for j = 1:n_t
                        correlation_lengths(i,j) = computeCorrelationLength_1D(nodes(:,1),di(:,j));
                    end
                end
                beta_u_new = mean(correlation_lengths(:),'omitnan')^2/12;
            elseif size(nodes,2) == 2
                correlation_lengths = zeros(N,n_t);
                for i = 1:N
                    di = reshape(this.data_interface.D(:,i),n_y,n_t);
                    for j = 1:n_t
                        correlation_lengths(i,j) = computeCorrelationLength_2D(nodes(:,1),nodes(:,2),di(:,j));
                    end
                end
                beta_u_new = mean(correlation_lengths(:),'omitnan')^2/8;
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

            nodes = this.Load_Spatial_Node_Data();
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

            N = size(this.data_interface.Z,2);
            if N > 1
                this.z_pert_norm_sq = zeros(N-1,1);
                for k = 2:N
                    v = this.data_interface.Z(:,k) - this.data_interface.Z(:,1);
                    this.z_pert_norm_sq(k-1) = v'*z_prior_interface.Apply_M_z(v);
                end
            end
            
            this.discrepancy_percent_z_variation = sqrt(max(this.d_pert_norm_sq./this.z_pert_norm_sq))*1.25;
            alpha_z_new = (this.discrepancy_percent_z_variation.^2) / (tmp * zopt_norm);
            this.Set_alpha_z(alpha_z_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_beta_z(this,beta_z_new)
            this.beta_z = beta_z_new;
        end

        function [] = Determine_beta_z(this)
            nodes = this.Load_Spatial_Node_Data();
            if size(nodes,2) == 1
                correlation_lengths = zeros(size(this.data_interface.Z,2),1);
                for k = 1:length(correlation_lengths)
                    correlation_lengths(k) = computeCorrelationLength_1D(nodes(:,1),this.data_interface.Z(:,k));
                end
                beta_z_new = mean(correlation_lengths,'omitnan')^2/12;
            elseif size(nodes,2) == 2
                correlation_lengths = zeros(size(this.data_interface.Z,2),1);
                for k = 1:length(correlation_lengths)
                    correlation_lengths(k) = computeCorrelationLength_2D(nodes(:,1),nodes(:,2),this.data_interface.Z(:,k));
                end
                beta_z_new = mean(correlation_lengths,'omitnan')^2/8;
            else
                disp('Determine_beta_z error: Dimensions greater than 2 are not supported.')
            end
            this.Set_beta_z(beta_z_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_beta_t(this,beta_t_new)
            this.beta_t = beta_t_new;
        end

        function [] = Determine_beta_t(this)
            time_nodes = this.Load_Time_Node_Data();
            n_t = length(time_nodes);
            n_y = size(this.data_interface.D(:,1),1)/n_t;
            N = size(this.data_interface.D,2);

            correlation_lengths = zeros(N,n_y);
            for i = 1:N
                di = reshape(this.data_interface.D(:,i),n_y,n_t)';
                for j = 1:n_y
                    correlation_lengths(i,j) = computeCorrelationLength_1D(time_nodes,di(:,j));
                end
            end
            beta_t_new = mean(correlation_lengths(:),'omitnan')^2/4;

            this.Set_beta_t(beta_t_new);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [] = Set_GSVD_Hyperparameters(this,gsvd_num_sing_vals_new,gsvd_oversampling_new,gsvd_num_subspace_iter_new)
            this.gsvd_num_sing_vals = gsvd_num_sing_vals_new;
            this.gsvd_oversampling = gsvd_oversampling_new;
            this.gsvd_num_subspace_iter = gsvd_num_subspace_iter_new;
        end

        function [] = Determine_GSVD_Hyperparameters(this)
            nodes = this.Load_Spatial_Node_Data();
            m = size(nodes,1);
            d = size(nodes,2);
            modes_per_dim = zeros(d,1);
            for k = 1:d
                L = max(nodes(:,k)) - min(nodes(:,k));
                tmp = (L/pi)^2 * (1/this.beta_u) * (1/this.W_u_inv_spectral_gap - 1);
                modes_per_dim(k) = round(sqrt(tmp));
            end
            num_sing_vals = prod(modes_per_dim);
            this.gsvd_num_sing_vals = min(num_sing_vals,m);
            this.gsvd_oversampling = min(10,m-this.gsvd_num_sing_vals);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function this = MD_Hyperparameters(data_interface,is_transient)
            this.data_interface = data_interface;
            this.is_transient = is_transient;

            this.data_noise_percent = 0.001;
            this.discrepancy_percent_z_variation = 1.0;
            this.W_u_inv_spectral_gap = 1.e-4;

            this.alpha_d = 0.0;
            this.alpha_u = 0.0;
            this.beta_u = 0.0;
            this.alpha_z = 0.0;
            this.beta_z = 0.0;
            this.alpha_t = 0.0;
            this.beta_t = 0.0;
            
            this.gsvd_num_sing_vals = 0;
            this.gsvd_oversampling = 0;
            this.gsvd_num_subspace_iter = 1;
        end
    end
end