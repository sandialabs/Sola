classdef MD_u_Hyperparameter_Interface < handle

    properties
        is_transient
        center_data
        adapt_time_variance
        component_id
        trace_estimator_sample_size

        alpha_u
        beta_u
        alpha_t
        beta_t
        alpha_d

        gsvd_num_sing_vals
        gsvd_oversampling
        gsvd_num_subspace_iter

        time_variance_inflation
        data_noise_percent
        W_u_inv_spectral_gap

    end

    %% Virtual functions for user implementation
    methods

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = [];
            disp('Load_Spatial_Node_Data is required for hyperparameter algorithm-based initialization');
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = [];
            disp('Load_Time_Node_Data is required for hyperparameter algorithm-based initialization');
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_u_Hyperparameter_Interface(is_transient, center_data, adapt_time_variance, component_id)
            arguments
                is_transient {logical}
                center_data {logical} = false
                adapt_time_variance {logical} = false
                component_id (1, 1) {mustBeNumeric} = 1
            end

            this.is_transient = is_transient;
            this.center_data = center_data;
            this.adapt_time_variance = adapt_time_variance;
            this.component_id = component_id;
            this.trace_estimator_sample_size = 0;

            this.alpha_u = 0.0;
            this.beta_u = 0.0;
            this.alpha_t = 1.0;
            this.beta_t = 0.0;
            this.alpha_d = 0.0;

            this.gsvd_num_sing_vals = 0;
            this.gsvd_oversampling = 0;
            this.gsvd_num_subspace_iter = 1;

            this.time_variance_inflation = .01;
            this.data_noise_percent = 0.001;
            this.W_u_inv_spectral_gap = 1.e-4;
        end

        function [] = Set_alpha_u(this, alpha_u_new)
            this.alpha_u = alpha_u_new;
        end

        function [] = Set_beta_u(this, beta_u_new)
            this.beta_u = beta_u_new;
        end

        function [] = Set_alpha_t(this, alpha_t_new)
            this.alpha_t = alpha_t_new;
        end

        function [] = Set_beta_t(this, beta_t_new)
            this.beta_t = beta_t_new;
        end

        function [] = Set_alpha_d(this, alpha_d_new)
            this.alpha_d = alpha_d_new;
        end

        function [] = Set_GSVD_Hyperparameters(this, gsvd_num_sing_vals_new, gsvd_oversampling_new, gsvd_num_subspace_iter_new)
            this.gsvd_num_sing_vals = gsvd_num_sing_vals_new;
            this.gsvd_oversampling = gsvd_oversampling_new;
            this.gsvd_num_subspace_iter = gsvd_num_subspace_iter_new;
        end

    end

end
