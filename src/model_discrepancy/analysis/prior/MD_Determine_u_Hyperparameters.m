classdef MD_Determine_u_Hyperparameters < handle

    properties
        data_interface
        u_hyperparam_interface

        component_id
        is_transient
    end

    methods

        function [] = Determine_alpha_u(this, u_prior_interface)
            I = this.data_interface.Separate_State_Components(this.component_id);
            delta1 = this.data_interface.D(I, 1);
            this.u_hyperparam_interface.d1_norm_sq = delta1' * u_prior_interface.Apply_M_u(delta1);

            N = size(this.data_interface.D, 2);
            if N > 1
                this.u_hyperparam_interface.d_pert_norm_sq = zeros(N - 1, 1);
                for k = 2:N
                    v = this.data_interface.D(I, k) - delta1;
                    this.u_hyperparam_interface.d_pert_norm_sq(k - 1) = v' * u_prior_interface.Apply_M_u(v);
                end
            end

            if this.is_transient
                u_op_trace = sum(u_prior_interface.spatial_prior_cov.sing_vals.^2) * sum(u_prior_interface.transient_prior_cov.evals);
            else
                u_op_trace = sum(u_prior_interface.sing_vals.^2);
            end

            alpha_u_new = this.u_hyperparam_interface.d1_norm_sq / u_op_trace;
            this.u_hyperparam_interface.Set_alpha_u(alpha_u_new);
        end

        function [] = Determine_beta_u(this)
            I = this.data_interface.Separate_State_Components(this.component_id);
            nodes = this.u_hyperparam_interface.Load_Spatial_Node_Data();
            nodes = nodes{this.component_id};
            n_y = size(nodes, 1);
            n_t = size(this.data_interface.D(I, 1), 1) / n_y;
            N = size(this.data_interface.D, 2);
            initial_guess = 0;

            if size(nodes, 2) == 1
                corr_len_fun = @(nodes, d, initial_guess) computeCorrelationLength_1D(nodes(:, 1), d, initial_guess);
                normalization = 12;
            elseif size(nodes, 2) == 2
                corr_len_fun = @(nodes, d, initial_guess) computeCorrelationLength_2D(nodes(:, 1), nodes(:, 2), d, initial_guess);
                normalization = 8;
            else
                disp('Determine_beta_u error: Dimensions greater than 2 are not supported.');
            end

            correlation_lengths = zeros(N, 1);
            for i = 1:N
                di = mean(reshape(this.data_interface.D(I, i), n_y, n_t),2);
                correlation_lengths(i) = corr_len_fun(nodes, di, initial_guess);
                initial_guess = correlation_lengths(i);
            end
            beta_u_new = mean(correlation_lengths, 'omitnan')^2 / normalization;

            this.u_hyperparam_interface.Set_beta_u(beta_u_new);
        end

        function [] = Determine_alpha_t(this, u_prior_interface)
            time_nodes = this.u_hyperparam_interface.Load_Time_Node_Data();
            I = this.data_interface.Separate_State_Components(this.component_id);
            n_t = length(time_nodes);
            n_y = size(this.data_interface.D(I, 1), 1) / n_t;
            N = size(this.data_interface.D, 2);
            alpha_t_new = zeros(n_t, N);
            for i = 1:N
                di = reshape(this.data_interface.D(I, i), n_y, n_t);
                tmp = diag(di' * u_prior_interface.spatial_prior_cov.Apply_M_u(di));
                tmp = tmp / max(tmp);
                tmp = tmp + (1.e-2);
                tmp = tmp / (1 + 1.e-2);
                alpha_t_new(:, i) = tmp;
            end

            alpha_t_new = mean(alpha_t_new, 2);

            this.u_hyperparam_interface.Set_alpha_t(alpha_t_new);
        end

        function [] = Determine_beta_t(this)
            time_nodes = this.u_hyperparam_interface.Load_Time_Node_Data();
            I = this.data_interface.Separate_State_Components(this.component_id);
            n_t = length(time_nodes);
            n_y = size(this.data_interface.D(I, 1), 1) / n_t;
            N = size(this.data_interface.D, 2);

            correlation_lengths = zeros(N, n_y);
            initial_guess = 0;
            for i = 1:N
                di = reshape(this.data_interface.D(I, i), n_y, n_t)';
                for j = 1:n_y
                    correlation_lengths(i, j) = computeCorrelationLength_1D(time_nodes, di(:, j), initial_guess);
                    initial_guess = correlation_lengths(i, j);
                end
                initial_guess = correlation_lengths(i, 1);
            end
            beta_t_new = mean(correlation_lengths(:), 'omitnan')^2 / 4;

            this.u_hyperparam_interface.Set_beta_t(beta_t_new);
        end

        function [] = Determine_alpha_d(this)
            I = this.data_interface.Separate_State_Components(this.component_id);
            tmp = this.data_interface.D(I, :);
            alpha_d_new = (this.u_hyperparam_interface.data_noise_percent * mean(abs(tmp(:))))^2;
            this.u_hyperparam_interface.Set_alpha_d(alpha_d_new);
        end

        function [] = Determine_GSVD_Hyperparameters(this)
            nodes = this.u_hyperparam_interface.Load_Spatial_Node_Data();
            nodes = nodes{this.component_id};
            m = size(nodes, 1);
            d = size(nodes, 2);
            modes_per_dim = zeros(d, 1);
            for k = 1:d
                L = max(nodes(:, k)) - min(nodes(:, k));
                tmp = (L / pi)^2 * (1 / this.u_hyperparam_interface.beta_u) * (1 / this.u_hyperparam_interface.W_u_inv_spectral_gap - 1);
                modes_per_dim(k) = round(sqrt(tmp));
            end
            num_sing_vals = prod(modes_per_dim);
            gsvd_num_sing_vals = min(num_sing_vals, m);
            gsvd_oversampling = min(10, m - gsvd_num_sing_vals);
            this.u_hyperparam_interface.Set_GSVD_Hyperparameters(gsvd_num_sing_vals, gsvd_oversampling, 1);
        end

        function [] = Determine_Data_Centering(this)
            I = this.data_interface.Separate_State_Components(this.component_id);
            data_shift = mean(this.data_interface.D(I, 1)) * ones(size(this.data_interface.D, 1), 1);
            this.data_interface.data_shift = data_shift;
            this.data_interface.D = this.data_interface.D - data_shift;
        end

        function this = MD_Determine_u_Hyperparameters(data_interface, u_hyperparam_interface)
            arguments
                data_interface MD_Data_Interface
                u_hyperparam_interface MD_u_Hyperparameter_Interface
            end

            this.data_interface = data_interface;
            this.u_hyperparam_interface = u_hyperparam_interface;

            this.component_id = u_hyperparam_interface.component_id;
            this.is_transient = u_hyperparam_interface.is_transient;

            this.data_interface.Load_Data();
            if this.u_hyperparam_interface.center_data
                this.Determine_Data_Centering();
            end
            this.Determine_alpha_d();
        end

    end
end
