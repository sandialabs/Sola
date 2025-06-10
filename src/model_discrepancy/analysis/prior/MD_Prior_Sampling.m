classdef MD_Prior_Sampling < handle

    properties
        data_interface
        u_prior_interface
        z_prior_interface
        u_opt
        z_opt
        z_pert_subsample_factor

        delta_samples_z_opt
        delta_samples_z_pert
        z_pert
        z_pert_evals
        z_pert_corr
        z_pert_corr_binned
        z_pert_corr_bin_means
        delta_mag
        delta_corr
        d1_mag
        d1_corr
        temporal_mag
        temporal_corr_len
        delta_pert_mag
        delta_pert_mag_binned

        delta_z_opt_time_evol
        data_time_evol
    end

    %% Constructor and user access functions
    methods

        function this = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface)
            arguments
                data_interface MD_Data_Interface
                u_prior_interface MD_u_Prior_Interface
                z_prior_interface MD_z_Prior_Interface
            end
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.u_opt = this.data_interface.u_opt;
            this.z_opt = this.data_interface.z_opt;
            this.z_pert_subsample_factor = 1;
        end

        function [] = Generate_Prior_Discrepancy_Sample_Data(this, num_samps, econ)
            if nargin < 3
                econ = false;
            end
            this.Generate_Prior_Discrepancy_z_opt_Sample_Data(num_samps, econ);
            this.Generate_Prior_Discrepancy_z_pert_Sample_Data(econ);
        end

        function [] = Generate_Prior_Discrepancy_z_opt_Sample_Data(this, num_samps, econ)
            if nargin < 3
                econ = false;
            end

            this.delta_samples_z_opt = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps)  + this.data_interface.data_shift;
            if ~econ
                this.Compute_Delta_z_opt_Metrics();
            end
            if this.u_prior_interface.u_hyperparam_interface.is_transient
                this.Compute_Temporal_Data(econ);
            end
        end

        function [] = Generate_Prior_Discrepancy_z_pert_Sample_Data(this, econ)
            if nargin < 2
                econ = false;
            end

            this.Compute_z_pert_Data(econ);
            if ~econ
                this.Compute_Delta_z_pert_Metrics();
            end
        end

        function [delta_samples] = Prior_Discrepancy_Samples_at_z_opt(this, num_samps)
            delta_samples = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps) + this.data_interface.data_shift;
        end

        function [delta_samples, delta_zopt_samples] = Prior_Discrepancy_Samples(this, z, num_samps)
            Z = z - this.z_opt;
            Mz_Z = this.z_prior_interface.Apply_M_z(Z);
            Sigma = Mz_Z' * this.z_prior_interface.Apply_W_z_Inverse(Mz_Z);
            p = size(Z, 2);
            R = chol(Sigma);

            delta_zopt_samples = zeros(size(this.data_interface.D, 1), num_samps);
            delta_samples = cell(num_samps, 1);
            for k = 1:num_samps
                u_vec = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(p + 1);
                delta_samples{k} = u_vec(:, 1:p) * R + u_vec(:, p + 1) + this.data_interface.data_shift;
                delta_zopt_samples(:, k) = u_vec(:, p + 1) + this.data_interface.data_shift;
            end
        end

    end

    %% Functions to compute data for prior visualization
    methods

        function [] = Compute_Temporal_Data(this, econ)
            if nargin < 2
                econ = false;
            end

            num_components = length(this.delta_mag);
            if isempty(this.delta_z_opt_time_evol)
                this.delta_z_opt_time_evol = cell(num_components, 1);
            end
            num_samps = size(this.delta_samples_z_opt, 2);

            t = this.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
            n_t = length(t);
            n_y = length(this.delta_samples_z_opt(:, 1)) / n_t;
            for component_id = 1:num_components
                this.delta_z_opt_time_evol{component_id} = zeros(num_samps, n_t);
            end

            for i = 1:num_samps
                di = reshape(this.delta_samples_z_opt(:, i), n_y, n_t);
                tmp = this.u_prior_interface.Apply_M_u(di);
                if num_components > 1
                    for component_id = 1:num_components
                        J = this.u_prior_interface.I{component_id};
                        J = J(1:this.u_prior_interface.u_prior_interface_cell{component_id}.transient_prior_cov.n_y);
                        this.delta_z_opt_time_evol{component_id}(i, :) = sqrt(diag(di(J, :)' * tmp(J, :)));
                    end
                else
                    this.delta_z_opt_time_evol{1}(i, :) = sqrt(diag(di' * tmp));
                end
            end

            this.data_time_evol = cell(num_components, 1);
            d = reshape(this.data_interface.D(:, 1), n_y, n_t);
            tmp = this.u_prior_interface.Apply_M_u(d);
            if num_components > 1
                for component_id = 1:num_components
                    J = this.u_prior_interface.I{component_id};
                    J = J(1:this.u_prior_interface.u_prior_interface_cell{component_id}.transient_prior_cov.n_y);
                    this.data_time_evol{component_id} = sqrt(diag(d(J, :)' * tmp(J, :)));
                end
            else
                this.data_time_evol{1} = sqrt(diag(d' * tmp));
            end

            if ~econ
                this.temporal_mag = cell(num_components, 1);
                this.temporal_corr_len = cell(num_components, 1);
                for component_id = 1:num_components
                    this.temporal_corr_len{component_id} = zeros(num_samps, 1);
                    initial_guess = 0;
                    for i = 1:num_samps
                        this.temporal_corr_len{component_id}(i) = Compute_Correlation_Length_1D(t, this.delta_z_opt_time_evol{component_id}(i, :), initial_guess);
                        initial_guess = this.temporal_corr_len{component_id}(i);
                    end
                    this.temporal_mag{component_id} = max(abs(this.delta_z_opt_time_evol{component_id}), [], 2);
                end
            end
        end

        function [] = Compute_z_pert_Data(this, econ)

            if nargin < 2
                econ = false;
            end

            if ~econ

                if strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'vector')
                    this.z_pert = eye(this.z_prior_interface.n_z);
                    this.z_pert_evals = ones(this.z_prior_interface.n_z, 1);
                else
                    e = this.z_prior_interface.determine_z_hyperparams.Compute_Eigenvalues(this.z_prior_interface);
                    num_perts_init = length(find(1 ./ e > .1));
                    if strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'transient vector')
                        num_perts_init = this.z_prior_interface.num_controls * num_perts_init;
                    end
                    oversampling = min(10, length(this.data_interface.z_opt) - num_perts_init);
                    num_subspace_iters = 10;
                    E_z_inv_gsvd = E_z_Inv_GSVD(this.z_prior_interface, this.z_opt);
                    [this.z_pert, ~, this.z_pert_evals] = E_z_inv_gsvd.Compute_GSVD(num_perts_init, oversampling, num_subspace_iters);
                    while this.z_pert_evals(end) > .1
                        num_perts_init = 2 * num_perts_init;
                        oversampling = min(10, length(this.data_interface.z_opt) - num_perts_init);
                        [this.z_pert, ~, this.z_pert_evals] = E_z_inv_gsvd.Compute_GSVD(num_perts_init, oversampling, num_subspace_iters);
                    end
                end

                I = find(this.z_pert_evals > .1);
                I = round(linspace(I(1), I(end), round(length(I) / this.z_pert_subsample_factor)));
                this.z_pert = this.z_pert(:, I);
                this.z_pert_evals = this.z_pert_evals(I);
                scaling = .3 * sqrt(this.z_opt' * this.z_prior_interface.Apply_M_z(this.z_opt));
                this.z_pert = scaling * this.z_pert;

                num_perts = length(this.z_pert_evals);
                num_samps = size(this.delta_samples_z_opt, 2);
                this.delta_samples_z_pert = cell(num_perts, 1);
                for k = 1:num_perts
                    this.delta_samples_z_pert{k} = scaling * sqrt(this.z_prior_interface.alpha_z) * this.z_pert_evals(k) * this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
                end

            else
                if strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'spatial field')
                    this.z_pert = zeros(length(this.z_opt), 2);
                    this.z_pert_evals = zeros(2, 1);
                    scaling = .3 * sqrt(this.z_opt' * this.z_prior_interface.Apply_M_z(this.z_opt));

                    v = 0 * this.z_opt + 1;
                    tmp = sqrt(v' * this.z_prior_interface.Apply_M_z(v));
                    this.z_pert(:, 1) = (scaling / tmp) * v;
                    this.z_pert_evals(1) = 1;

                    x = this.z_prior_interface.z_hyperparam_interface.x;
                    n = size(x, 2);
                    L = zeros(n, 1);
                    for k = 1:n
                        x0 = min(x(:, k));
                        x1 = max(x(:, k));
                        L(k) = x1 - x0;
                        x(:, k) = 2 * pi * (x(:, k) - x0) / (x1 - x0);
                    end

                    if n == 1
                        omega = round((3 / (2 * pi)) * (L(1) / sqrt(this.z_prior_interface.beta_z)));
                        v = cos(omega * x);
                    elseif n == 2
                        omega = round((3 / (2 * pi)) * (1 / sqrt(this.z_prior_interface.beta_z)) / sqrt(1 / L(1)^2 + 1 / L(2)^2));
                        v = cos(omega * x(:, 1)) .* cos(omega * x(:, 2));
                    elseif n == 3
                        omega = round((3 / (2 * pi)) * (1 / sqrt(this.z_prior_interface.beta_z)) / sqrt(1 / L(1)^2 + 1 / L(2)^2 + 1 / L(3)^2));
                        v = cos(omega * x(:, 1)) .* cos(omega * x(:, 2)) .* cos(omega * x(:, 3));
                    end
                    tmp = sqrt(v' * this.z_prior_interface.Apply_M_z(v));
                    v = v / tmp;
                    this.z_pert(:, 2) = scaling * v;
                    this.z_pert_evals(2) = 1 / (v' * this.z_prior_interface.Apply_E_z(v));

                    num_samps = size(this.delta_samples_z_opt, 2);
                    this.delta_samples_z_pert = cell(2, 1);
                    for k = 1:2
                        this.delta_samples_z_pert{k} = scaling * sqrt(this.z_prior_interface.alpha_z) * this.z_pert_evals(k) * this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
                    end
                else
                    disp('Need to implement econ mode for z_type = "vector"');
                end

            end

        end

        function [] = Compute_Delta_z_opt_Metrics(this)

            nodes = this.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
            num_components = length(nodes);

            this.delta_mag = cell(num_components, 1);
            this.delta_corr = cell(num_components, 1);
            for k = 1:num_components
                I = this.data_interface.Separate_State_Components(k);
                n_y = size(nodes{k}, 1);
                n_t = size(this.delta_samples_z_opt(I, :), 1) / n_y;

                num_samps = size(this.delta_samples_z_opt(I, :), 2);
                this.delta_mag{k} = max(abs(this.delta_samples_z_opt(I, :)))';
                this.d1_mag{k} = max(abs(this.data_interface.D(I, 1) + this.data_interface.data_shift(I)));

                initial_guess = 0;
                correlation_lengths = zeros(num_samps, 1);
                if size(nodes{k}, 2) == 1
                    corr_len_fun = @(nodes, d, initial_guess) Compute_Correlation_Length_1D(nodes(:, 1), d, initial_guess);
                elseif size(nodes{k}, 2) == 2
                    corr_len_fun = @(nodes, d, initial_guess) Compute_Correlation_Length_2D(nodes(:, 1), nodes(:, 2), d, initial_guess);
                end

                for i = 1:num_samps
                    di = mean(reshape(this.delta_samples_z_opt(I, i), n_y, n_t), 2);
                    correlation_lengths(i) = corr_len_fun(nodes{k}, di, initial_guess);
                    initial_guess = correlation_lengths(i);
                end
                this.delta_corr{k} = correlation_lengths;
                d = mean(reshape(this.data_interface.D(I, 1) + this.data_interface.data_shift(I), n_y, n_t), 2);
                this.d1_corr{k} = corr_len_fun(nodes{k}, d, mean(correlation_lengths));
            end
        end

        function Compute_Delta_z_pert_Metrics(this)

            num_components = length(this.delta_mag);
            num_perts = size(this.z_pert, 2);

            this.delta_pert_mag = cell(num_components, num_perts);
            for k = 1:num_components
                I = this.data_interface.Separate_State_Components(k);
                for j = 1:num_perts
                    this.delta_pert_mag{k, j} = max(abs(this.delta_samples_z_pert{j}(I, :)))';
                end
            end

            if ~strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'vector')

                if strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'spatial field')
                    nodes = this.z_prior_interface.z_hyperparam_interface.Load_Spatial_Node_Data();
                    if size(nodes, 2) == 1
                        corr_len_fun = @(nodes, z, initial_guess) Compute_Correlation_Length_1D(nodes(:, 1), z, initial_guess);
                    elseif size(nodes, 2) == 2
                        corr_len_fun = @(nodes, z, initial_guess) Compute_Correlation_Length_2D(nodes(:, 1), nodes(:, 2), z, initial_guess);
                    end
                elseif strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'transient vector')
                    nodes = this.z_prior_interface.z_hyperparam_interface.Load_Time_Node_Data();
                    corr_len_fun = @(t, z, initial_guess) Compute_Correlation_Length_Transient_Vector(t, z, initial_guess);
                end
                this.z_pert_corr = zeros(num_perts, 1);
                initial_guess = 0;
                for j = 1:num_perts
                    this.z_pert_corr(j) = corr_len_fun(nodes, this.z_pert(:, j), initial_guess);
                    initial_guess = this.z_pert_corr(j);
                end

                bin_width = 0.05;
                bin_centers = (min(this.z_pert_corr) - bin_width / 2):bin_width:(max(this.z_pert_corr) + bin_width / 2);
                num_bins = length(bin_centers);
                this.z_pert_corr_binned = cell(num_bins, 1);
                this.delta_pert_mag_binned = cell(num_components, num_bins);
                for i = 1:num_bins
                    this.z_pert_corr_binned{i} = [];
                    for k = 1:num_components
                        this.delta_pert_mag_binned{k, i} = [];
                    end
                end
                for j = 1:num_perts
                    [~, i] = min(abs(bin_centers - this.z_pert_corr(j)));
                    this.z_pert_corr_binned{i} = [this.z_pert_corr_binned{i}; this.z_pert_corr(j)];
                    for k = 1:num_components
                        this.delta_pert_mag_binned{k, i} = [this.delta_pert_mag_binned{k, i}; this.delta_pert_mag{k, j}];
                    end
                end
                J = ~cellfun('isempty', this.z_pert_corr_binned);
                this.z_pert_corr_binned = this.z_pert_corr_binned(J);
                this.z_pert_corr_bin_means = cellfun(@(x)mean(x), this.z_pert_corr_binned);
                this.delta_pert_mag_binned = this.delta_pert_mag_binned(:, J);
            end
        end

    end

end
