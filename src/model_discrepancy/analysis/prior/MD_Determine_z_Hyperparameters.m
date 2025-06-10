classdef MD_Determine_z_Hyperparameters < handle

    properties
        data_interface
        z_hypeparam_interface
        u_prior_interface
        z_type
    end

    %% Constructor
    methods

        function this = MD_Determine_z_Hyperparameters(data_interface, z_hypeparam_interface, u_prior_interface)
            arguments
                data_interface MD_Data_Interface
                z_hypeparam_interface MD_z_Hyperparameter_Interface
                u_prior_interface MD_u_Prior_Interface
            end

            this.data_interface = data_interface;
            this.z_hypeparam_interface = z_hypeparam_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_type = this.z_hypeparam_interface.z_type;
        end

    end

    %% Functions to determine hyperparameters
    methods

        function [] = Determine_alpha_z(this, z_prior_interface)
            tmp = z_prior_interface.Apply_M_z(this.data_interface.z_opt);
            zopt_norm = sqrt(tmp' * this.data_interface.z_opt);

            e = this.Compute_Eigenvalues(z_prior_interface);
            evals = sort(1 ./ e, 'descend');
            I = find(evals < 1.e-3);
            rank = length(evals);
            if ~isempty(I)
                rank = I(1) - 1;
            end
            evals = evals(1:rank);
            samples = 1000;
            nu = randn(rank, samples)';
            nu = nu.^2;
            tmp = mean((nu * evals.^4) ./ (nu * evals.^2));

            if this.z_hypeparam_interface.discrepancy_percent_z_variation == 1
                if this.z_hypeparam_interface.num_state_solves > 0
                    u_nom = this.z_hypeparam_interface.State_Solve(this.data_interface.z_opt);
                    z_samples = z_prior_interface.Sample_with_Covariance_W_z_Acute_Inverse(this.z_hypeparam_interface.num_state_solves);
                    mags = sqrt(diag(z_samples' * z_prior_interface.Apply_M_z(z_samples)));
                    z_samples = zopt_norm * z_samples * diag(1 ./ mags);
                    u_samples = this.z_hypeparam_interface.State_Solve(z_samples + this.data_interface.z_opt);
                    e = u_samples - u_nom;
                    e_norm_sqr = diag(e' * this.u_prior_interface.Apply_M_u(e));
                    d1 = this.data_interface.D(:, 1);
                    d1_norm_sq = d1' * this.u_prior_interface.Apply_M_u(d1);
                    this.z_hypeparam_interface.discrepancy_percent_z_variation = sqrt(mean(e_norm_sqr) / d1_norm_sq);
                else
                    N = size(this.data_interface.Z, 2);
                    if N > 1
                        z1_norm_sq = this.data_interface.Z(:, 1)' * z_prior_interface.Apply_M_z(this.data_interface.Z(:, 1));
                        z_pert_norm_sq = zeros(N - 1, 1);
                        for k = 2:N
                            v = this.data_interface.Z(:, k) - this.data_interface.Z(:, 1);
                            z_pert_norm_sq(k - 1) = v' * z_prior_interface.Apply_M_z(v);
                        end
                        d1 = this.data_interface.D(:, 1);
                        d1_norm_sq = d1' * this.u_prior_interface.Apply_M_u(d1);
                        d_pert_norm_sq = zeros(N - 1, 1);
                        for k = 2:N
                            d = this.data_interface.D(:, k) - d1;
                            d_pert_norm_sq(k - 1) = d' * this.u_prior_interface.Apply_M_u(d);
                        end
                        this.z_hypeparam_interface.discrepancy_percent_z_variation = mean(sqrt((d_pert_norm_sq / d1_norm_sq) ./ (z_pert_norm_sq / z1_norm_sq)));
                    end
                end
            end

            alpha_z_new = (this.z_hypeparam_interface.discrepancy_percent_z_variation.^2) / (tmp * zopt_norm^2);
            this.z_hypeparam_interface.Set_alpha_z(alpha_z_new);
        end

        function [e] = Compute_Eigenvalues(this, z_prior_interface)

            if strcmp(this.z_type, 'spatial field')
                nodes = this.z_hypeparam_interface.Load_Spatial_Node_Data();
                if size(nodes, 2) == 1
                    Lx = max(nodes(:, 1)) - min(nodes(:, 1));
                    n = length(nodes(:, 1)) - 1;
                    e = 1 + z_prior_interface.beta_z * (pi / Lx)^2 * (0:n).^2;
                    e = e';
                elseif size(nodes, 2) == 2
                    Lx = max(nodes(:, 1)) - min(nodes(:, 1));
                    Ly = max(nodes(:, 2)) - min(nodes(:, 2));
                    n = round(sqrt(length(nodes(:, 1)))) - 1;
                    e = 1 + z_prior_interface.beta_z * pi^2 * (kron(((0:n).^2)', ones(n + 1, 1)) / Lx^2 + kron(ones(n + 1, 1), ((0:n).^2)') / Ly^2);
                else
                    disp('Determine_alpha_z error: Dimensions greater than 2 are not supported.');
                end
            end

            if strcmp(this.z_type, 'transient vector')
                t = this.z_hypeparam_interface.Load_Time_Node_Data();
                L = max(t) - min(t);
                n = length(t) - 1;
                e = 1 + z_prior_interface.beta_t * (pi / L)^2 * (0:n).^2;
                e = sqrt(e');
            end

            if strcmp(this.z_type, 'vector')
                e = ones(length(this.data_interface.z_opt), 1);
            end
        end

        function [] = Determine_beta_z(this)
            nodes = this.z_hypeparam_interface.Load_Spatial_Node_Data();
            initial_guess = 0;
            if size(nodes, 2) == 1
                correlation_lengths = zeros(size(this.data_interface.Z, 2), 1);
                for k = 1:length(correlation_lengths)
                    correlation_lengths(k) = Compute_Correlation_Length_1D(nodes(:, 1), this.data_interface.Z(:, k), initial_guess);
                    initial_guess = correlation_lengths(k);
                end
                beta_z_new = mean(correlation_lengths, 'omitnan')^2 / 12;
            elseif size(nodes, 2) == 2
                correlation_lengths = zeros(size(this.data_interface.Z, 2), 1);
                for k = 1:length(correlation_lengths)
                    correlation_lengths(k) = Compute_Correlation_Length_2D(nodes(:, 1), nodes(:, 2), this.data_interface.Z(:, k), initial_guess);
                    initial_guess = correlation_lengths(k);
                end
                beta_z_new = mean(correlation_lengths, 'omitnan')^2 / 8;
            else
                disp('Determine_beta_z error: Dimensions greater than 2 are not supported.');
            end
            this.z_hypeparam_interface.Set_beta_z(beta_z_new);
        end

        function [] = Determine_beta_t(this)
            t = this.z_hypeparam_interface.Load_Time_Node_Data();
            initial_guess = 0;
            N = size(this.data_interface.Z, 2);
            n_t = length(t);
            num_comp = size(this.data_interface.Z, 1) / n_t;
            correlation_lengths = zeros(N, num_comp);
            for k = 1:N
                z_data = reshape(this.data_interface.Z(:, k), num_comp, n_t)';
                for j = 1:num_comp
                    correlation_lengths(k, j) = Compute_Correlation_Length_1D(t, z_data(:, j), initial_guess);
                    initial_guess = correlation_lengths(k, j);
                end
            end
            beta_t_new = mean(correlation_lengths(:), 'omitnan')^2 / 4;

            this.z_hypeparam_interface.Set_beta_t(beta_t_new);
        end

    end
end
