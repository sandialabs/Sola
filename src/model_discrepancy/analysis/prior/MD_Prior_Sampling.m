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

    methods

        function this = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface)
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.u_opt = this.data_interface.u_opt;
            this.z_opt = this.data_interface.z_opt;
            this.z_pert_subsample_factor = 1;
        end

        %%
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

        %%
        function [] = Generate_Prior_Discrepancy_Sample_Data(this, num_samps)
            this.Generate_Prior_Discrepancy_z_opt_Sample_Data(num_samps);
            this.Generate_Prior_Discrepancy_z_pert_Sample_Data(num_samps);
        end

        function [] = Generate_Prior_Discrepancy_z_opt_Sample_Data(this, num_samps)
            this.delta_samples_z_opt = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps)  + this.data_interface.data_shift;
            this.Compute_Delta_z_opt_Metrics();
            if this.u_prior_interface.u_hyperparam_interface.is_transient
                this.Compute_temporal_data();
            end
        end

        function [] = Generate_Prior_Discrepancy_z_pert_Sample_Data(this, num_samps)
            this.Compute_z_pert();
            this.Compute_Delta_z_pert_Metrics();
        end

        %%
        function [] = Visualization_for_Prior_Discrepancy_at_z_pert(this, component_id)

            if ~strcmp(this.z_prior_interface.z_hyperparam_interface.z_type,'vector')

                I = this.data_interface.Separate_State_Components(component_id);
                num_perts = length(this.delta_samples_z_pert);
                nodes = this.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
                nodes = nodes{component_id};
                spatial_dim = size(nodes, 2);
                is_transient = this.u_prior_interface.u_hyperparam_interface.is_transient;
                if is_transient
                    t = this.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
                    t_z = this.z_prior_interface.z_hyperparam_interface.Load_Time_Node_Data();
                else
                    t = [];
                    t_z = [];
                end

                z_range = [min(this.z_pert(:)), max(this.z_pert(:))];

                if strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'transient vector')
                    generate_z_plot = @(nodes, t, u, z_range, name) this.Plot_Transient_Vector(t, u, z_range, name);
                elseif strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'spatial field')
                    if spatial_dim == 1
                        generate_z_plot = @(nodes, t, u, z_range, name) this.Plot_Stationary_1D(nodes, t, u, z_range, name);
                    elseif spatial_dim == 2
                        generate_z_plot = @(nodes, t, u, z_range, name) this.Plot_Stationary_2D(nodes, t, u, z_range, name);
                    end
                end

                if is_transient && (spatial_dim == 1)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Transient_1D(nodes, t, u, range, name);
                elseif is_transient && (spatial_dim == 2)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Transient_2D(nodes, t, u, range, name);
                elseif ~is_transient && (spatial_dim == 1)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Stationary_1D(nodes, t, u, range, name);
                elseif ~is_transient && (spatial_dim == 2)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Stationary_2D(nodes, t, u, range, name);
                end

                z_fig = figure;
                generate_z_plot(nodes, t_z, this.data_interface.Z(:, 1), z_range, '$z_1$');

                delta_fig = figure;
                c_range = [min(this.data_interface.D(I, 1) + this.data_interface.data_shift(I)),max(this.data_interface.D(I, 1) + this.data_interface.data_shift(I))];
                generate_delta_plot(nodes, t, this.data_interface.D(I, 1) + this.data_interface.data_shift(I), c_range, '$d_1$');

                margin = max(this.z_pert_corr)*.05;
                xmin = min(this.z_pert_corr)-margin;
                xmax = max(this.z_pert_corr)+margin;
                yrange = [0, 0];
                for k = 1:num_perts
                    yrange(1) = min(yrange(1), min(this.delta_pert_mag{component_id, k}));
                    yrange(2) = max(yrange(2), max(this.delta_pert_mag{component_id, k}));
                end
                ymin = yrange(1) * .9;
                ymax = yrange(2) * 1.1;
                colors = lines(length(this.z_pert_corr_binned));
                scatter_fig = figure;
                hold on;
                for k = 1:length(this.z_pert_corr_binned)
                    violinplot(this.z_pert_corr_bin_means(k),this.delta_pert_mag_binned{component_id, k},'DensityDirection','positive','DensityWidth',.05,'DensityScale','width','FaceColor','black')
                    plot(this.z_pert_corr_bin_means(k), this.delta_pert_mag_binned{component_id, k}, 'o','Color',colors(k,:),'MarkerSize',8);
                end
                xlabel('$\Delta z$ Correlation Length','Interpreter','latex');
                ylabel('$\delta(\Delta z)$ Magnitude','Interpreter','latex');
                xlim([xmin, xmax]);
                ylim([ymin, ymax]);
                set(gca, 'fontsize', 24);

                user_continue = true;
                bnd_check = @(x, y) (x < xmin) || (x > xmax) || (y < ymin) || (y > ymax);

                disp('Please select point number for rendering. Select a point outside of the plotting domain to terminate');
                while user_continue

                    figure(scatter_fig);
                    [x, y] = ginput(1);

                    if bnd_check(x, y)
                        user_continue = false;
                    else

                        bin = find(abs(this.z_pert_corr_bin_means - x) < .01);
                        J = find(this.z_pert_corr==this.z_pert_corr_binned{bin}(1));
                        for k = 2:length(this.z_pert_corr_binned{bin})
                            J = [J;find(this.z_pert_corr==this.z_pert_corr_binned{bin}(k))];
                        end
                        i = J(1);
                        [val, j] = min(abs(this.delta_pert_mag{component_id,i} - y));
                        for k = 2:length(J)
                            [valk, jk] = min(abs(this.delta_pert_mag{component_id,J(k)} - y));
                            if valk < val
                                val = valk;
                                j = jk;
                                i = J(k);
                            end
                        end

                        figure(scatter_fig);
                        if length(scatter_fig.Children(1).Children) > num_perts
                            delete(scatter_fig.Children(1).Children(1));
                        end
                        plot(this.z_pert_corr_bin_means(bin), this.delta_pert_mag{component_id,i}(j), 'x', 'MarkerSize', 20, 'color', 'black');

                        figure(z_fig);
                        generate_z_plot(nodes, t_z, this.z_pert(:, i), z_range,'$\Delta z_k$');

                        figure(delta_fig);
                        c_range = [min(this.delta_samples_z_pert{i}(I, j)),max(this.delta_samples_z_pert{i}(I, j))];
                        generate_delta_plot(nodes, t, this.delta_samples_z_pert{i}(I, j), c_range, '$\delta(\Delta z_k,\theta_i)$');

                        figure(delta_fig);
                        figure(z_fig);
                        figure(scatter_fig);
                    end
                end
                disp('Concluded interactive visualization');

            else

                I = this.data_interface.Separate_State_Components(component_id);
                num_perts = length(this.delta_samples_z_pert);
                nodes = this.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
                nodes = nodes{component_id};
                spatial_dim = size(nodes, 2);
                is_transient = this.u_prior_interface.u_hyperparam_interface.is_transient;
                if is_transient
                    t = this.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
                else
                    t = [];
                    t_z = [];
                end

                z_range = [min(this.z_pert(:)), max(this.z_pert(:))];
                generate_z_plot = @(nodes, t, u, z_range, name) this.Plot_Vector(u, z_range, name);

                if is_transient && (spatial_dim == 1)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Transient_1D(nodes, t, u, range, name);
                elseif is_transient && (spatial_dim == 2)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Transient_2D(nodes, t, u, range, name);
                elseif ~is_transient && (spatial_dim == 1)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Stationary_1D(nodes, t, u, range, name);
                elseif ~is_transient && (spatial_dim == 2)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Stationary_2D(nodes, t, u, range, name);
                end

                z_fig = figure;
                generate_z_plot(nodes, t_z, this.data_interface.Z(:, 1), z_range, '$z_1$');

                if is_transient && (spatial_dim == 1)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Transient_1D(nodes, t, u, range, name);
                elseif is_transient && (spatial_dim == 2)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Transient_2D(nodes, t, u, range, name);
                elseif ~is_transient && (spatial_dim == 1)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Stationary_1D(nodes, t, u, range, name);
                elseif ~is_transient && (spatial_dim == 2)
                    generate_delta_plot = @(nodes, t, u, range, name) this.Plot_Stationary_2D(nodes, t, u, range, name);
                end

                delta_fig = figure;
                c_range = [min(this.data_interface.D(I, 1) + this.data_interface.data_shift(I)),max(this.data_interface.D(I, 1) + this.data_interface.data_shift(I))];
                generate_delta_plot(nodes, t, this.data_interface.D(I, 1) + this.data_interface.data_shift(I), c_range, '$d_1$');

                xmin = 0;
                xmax = this.z_prior_interface.n_z + 1;
                yrange = [0, 0];
                for k = 1:num_perts
                    yrange(1) = min(yrange(1), min(this.delta_pert_mag{component_id, k}));
                    yrange(2) = max(yrange(2), max(this.delta_pert_mag{component_id, k}));
                end
                ymin = yrange(1) * .9;
                ymax = yrange(2) * 1.1;
                colors = lines(this.z_prior_interface.n_z);
                scatter_fig = figure;
                hold on;
                for k = 1:this.z_prior_interface.n_z
                    violinplot(k,this.delta_pert_mag{component_id, k},'DensityDirection','positive','DensityWidth',.05,'DensityScale','width','FaceColor','black')
                    plot(k, this.delta_pert_mag{component_id, k}, 'o','Color',colors(k,:),'MarkerSize',8);
                end
                xlabel('$z$ Component','Interpreter','latex');
                ylabel('$\delta(\Delta z)$ Magnitude','Interpreter','latex');
                xlim([xmin, xmax]);
                ylim([ymin, ymax]);
                set(gca, 'fontsize', 24);

                user_continue = true;
                bnd_check = @(x, y) (x < xmin) || (x > xmax) || (y < ymin) || (y > ymax);

                disp('Please select point number for rendering. Select a point outside of the plotting domain to terminate');
                while user_continue

                    figure(scatter_fig);
                    [x, y] = ginput(1);

                    if bnd_check(x, y)
                        user_continue = false;
                    else

                        i = find(abs((1:this.z_prior_interface.n_z)' - x) < .1);
                        [~, j] = min(abs(this.delta_pert_mag{component_id,i} - y));

                        figure(scatter_fig);
                        if length(scatter_fig.Children(1).Children) > num_perts
                            delete(scatter_fig.Children(1).Children(1));
                        end
                        plot(i, this.delta_pert_mag{component_id,i}(j), 'x', 'MarkerSize', 20, 'color', 'black');

                        figure(z_fig);
                        generate_z_plot(nodes, t_z, this.z_pert(:, i), z_range,'$\Delta z_k$');

                        figure(delta_fig);
                        c_range = [min(this.delta_samples_z_pert{i}(I, j)),max(this.delta_samples_z_pert{i}(I, j))];
                        generate_delta_plot(nodes, t, this.delta_samples_z_pert{i}(I, j), c_range, '$\delta(\Delta z_k,\theta_i)$');

                        figure(delta_fig);
                        figure(z_fig);
                        figure(scatter_fig);
                    end
                end
                disp('Concluded interactive visualization');

            end

        end

        function [] = Visualization_for_Prior_Discrepancy_at_z_opt(this, component_id)

            I = this.data_interface.Separate_State_Components(component_id);
            nodes = this.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
            nodes = nodes{component_id};
            spatial_dim = size(nodes, 2);
            is_transient = this.u_prior_interface.u_hyperparam_interface.is_transient;
            if is_transient
                t = this.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
            else
                t = [];
            end
            range = [min(min(this.delta_samples_z_opt(I,:))), max(max(this.delta_samples_z_opt(I,:)))];

            if is_transient && (spatial_dim == 1)
                generate_plot = @(nodes, t, u, name) this.Plot_Transient_1D(nodes, t, u, range, name);
            elseif is_transient && (spatial_dim == 2)
                generate_plot = @(nodes, t, u, name) this.Plot_Transient_2D(nodes, t, u, range, name);
            elseif ~is_transient && (spatial_dim == 1)
                generate_plot = @(nodes, t, u, name) this.Plot_Stationary_1D(nodes, t, u, range, name);
            elseif ~is_transient && (spatial_dim == 2)
                generate_plot = @(nodes, t, u, name) this.Plot_Stationary_2D(nodes, t, u, range, name);
            end

            data_fig = figure;
            generate_plot(nodes, t, this.data_interface.D(I, 1) + this.data_interface.data_shift(I), '$d_1$');

            sample_fig = figure;
            generate_plot(nodes, t, this.data_interface.D(I, 1) + this.data_interface.data_shift(I), '$d_1$');

            xmin = min([this.d1_corr{component_id};this.delta_corr{component_id}]) * .9;
            xmax = max([this.d1_corr{component_id};this.delta_corr{component_id}]) * 1.1;
            ymin = min([this.d1_mag{component_id};this.delta_mag{component_id}]) * .9;
            ymax = max([this.d1_mag{component_id};this.delta_mag{component_id}]) * 1.1;
            scatter_fig = figure;
            hold on;
            plot(this.delta_corr{component_id}, this.delta_mag{component_id}, 'o');
            plot(this.d1_corr{component_id},this.d1_mag{component_id},'*');
            xlabel('Correlation Length');
            ylabel('Magnitude');
            legend({'Prior Samples','Discrepancy Data'})
            xlim([xmin, xmax]);
            ylim([ymin, ymax]);
            set(gca, 'fontsize', 24);

            user_continue = true;
            bnd_check = @(x, y) (x < xmin) || (x > xmax) || (y < ymin) || (y > ymax);

            disp('Please select point number for rendering. Select a point outside of the plotting domain to terminate');
            while user_continue

                figure(scatter_fig);
                [x, y] = ginput(1);

                if bnd_check(x, y)
                    user_continue = false;
                else
                    [~, i] = min(vecnorm([this.delta_corr{component_id}, this.delta_mag{component_id}] - [x, y], 2, 2));
                    figure(sample_fig);
                    generate_plot(nodes, t, this.delta_samples_z_opt(I, i), '$\delta(\tilde{z},\theta_i)$');

                    figure(scatter_fig);
                    if length(scatter_fig.Children(1).Children) > 1
                        delete(scatter_fig.Children(1).Children(1));
                    end
                    plot(this.delta_corr{component_id}(i),this.delta_mag{component_id}(i), 'x', 'MarkerSize', 20, 'color', 'black');
                    legend({'Prior Samples','Discrepancy Data'})

                    figure(sample_fig);
                    figure(data_fig);
                    figure(scatter_fig);
                end
            end
            disp('Concluded interactive visualization');

        end

        function [] = Visualization_for_Prior_Time_Evolution(this, component_id, interactive_viz)

            delta_min = .95 * min([this.delta_z_opt_time_evol{component_id}(:);this.data_time_evol{component_id}]);
            delta_max = 1.05 * max([this.delta_z_opt_time_evol{component_id}(:);this.data_time_evol{component_id}]);
            t = this.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();

            sample_fig = figure;
            hold on;
            plot(t, this.delta_z_opt_time_evol{component_id}, 'LineWidth', 3, 'Color', [.9, .9, .9]);
            plot(t, this.data_time_evol{component_id}, 'LineWidth', 3, 'Color', 'red');
            xlabel('Time');
            ylabel('Discrepancy Spatial Norm');
            ylim([delta_min, delta_max]);
            set(gca, 'fontsize', 24);

            data_fig = figure;
            hold on;
            plot(t, this.data_time_evol{component_id}, 'LineWidth', 3, 'Color', 'red');
            xlabel('Time');
            ylabel('Discrepancy Spatial Norm');
            ylim([delta_min, delta_max]);
            set(gca, 'fontsize', 24);

            xmid = median(this.temporal_corr_len{component_id});
            ymid = median(this.temporal_mag{component_id});
            xmin = min(this.temporal_corr_len{component_id}) - 0.3 * xmid;
            xmax = max(this.temporal_corr_len{component_id}) + 0.3 * xmid;
            ymin = min(this.temporal_mag{component_id}) - 0.3 * ymid;
            ymax = max(this.temporal_mag{component_id}) + 0.3 * ymid;
            scatter_fig = figure;
            hold on;
            plot(this.temporal_corr_len{component_id}, this.temporal_mag{component_id}, 'o');
            xlabel('Temporal Correlation Length');
            ylabel('Magnitude');
            xlim([xmin, xmax]);
            ylim([ymin, ymax]);
            set(gca, 'fontsize', 24);

            user_continue = interactive_viz;
            bnd_check = @(x, y) (x < xmin) || (x > xmax) || (y < ymin) || (y > ymax);
            colors = lines(20);
            current_selection = 1;

            disp('Please select point number for rendering. Select a point outside of the plotting domain to terminate');
            while user_continue && (current_selection <= 20)

                [x, y] = ginput(1);

                if bnd_check(x, y)
                    user_continue = false;
                else

                    [~, i] = min(vecnorm([this.temporal_corr_len{component_id}, this.temporal_mag{component_id}] - [x, y], 2, 2));
                    plot(data_fig.CurrentAxes, t, this.delta_z_opt_time_evol{component_id}(i, :), 'LineWidth', 3, 'Color', colors(current_selection, :));
                    plot(scatter_fig.CurrentAxes, this.temporal_corr_len{component_id}(i), this.temporal_mag{component_id}(i), 'x', 'MarkerSize', 15, 'Color', colors(current_selection, :));

                    current_selection = current_selection + 1;
                    figure(sample_fig);
                    figure(data_fig);
                    figure(scatter_fig);
                end
            end
            disp('Concluded interactive visualization');
        end

        %%
        function [] = Compute_temporal_data(this)
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
                        this.delta_z_opt_time_evol{component_id}(i,:) = sqrt(diag(di(J,:)'*tmp(J,:)));
                    end
                else
                    this.delta_z_opt_time_evol{1}(i,:) = sqrt(diag(di' * tmp));
                end
            end

            this.data_time_evol = cell(num_components, 1);
            d = reshape(this.data_interface.D(:, 1), n_y, n_t);
            tmp = this.u_prior_interface.Apply_M_u(d);
            if num_components > 1
                for component_id = 1:num_components
                    J = this.u_prior_interface.I{component_id};
                    J = J(1:this.u_prior_interface.u_prior_interface_cell{component_id}.transient_prior_cov.n_y);
                    this.data_time_evol{component_id} = sqrt(diag(d(J,:)' * tmp(J,:)));
                end
            else
                this.data_time_evol{1} = sqrt(diag(d' * tmp));
            end

            this.temporal_mag = cell(num_components,1);
            this.temporal_corr_len = cell(num_components,1);
            for component_id = 1:num_components
                this.temporal_corr_len{component_id} = zeros(num_samps, 1);
                initial_guess = 0;
                for i = 1:num_samps
                    this.temporal_corr_len{component_id}(i) = computeCorrelationLength_1D(t, this.delta_z_opt_time_evol{component_id}(i, :), initial_guess);
                    initial_guess = this.temporal_corr_len{component_id}(i);
                end
                this.temporal_mag{component_id} = max(abs(this.delta_z_opt_time_evol{component_id}), [], 2);
            end
        end

        function [] = Compute_z_pert(this)


            if strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'vector')
                this.z_pert = eye(this.z_prior_interface.n_z);
                this.z_pert_evals = ones(this.z_prior_interface.n_z,1);
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
            I = round(linspace(I(1),I(end),round(length(I)/this.z_pert_subsample_factor)));
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
        end

        %%
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
                this.d1_mag{k} = max(abs(this.data_interface.D(I,1)+this.data_interface.data_shift(I)));

                initial_guess = 0;
                correlation_lengths = zeros(num_samps, 1);
                if size(nodes{k}, 2) == 1
                    corr_len_fun = @(nodes, d, initial_guess) computeCorrelationLength_1D(nodes(:, 1), d, initial_guess);
                elseif size(nodes{k}, 2) == 2
                    corr_len_fun = @(nodes, d, initial_guess) computeCorrelationLength_2D(nodes(:, 1), nodes(:, 2), d, initial_guess);
                end

                for i = 1:num_samps
                    di = mean(reshape(this.delta_samples_z_opt(I, i), n_y, n_t),2);
                    correlation_lengths(i) = corr_len_fun(nodes{k}, di, initial_guess);
                    initial_guess = correlation_lengths(i);
                end
                this.delta_corr{k} = correlation_lengths;
                d = mean(reshape(this.data_interface.D(I, 1)+this.data_interface.data_shift(I), n_y, n_t),2);
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
                        corr_len_fun = @(nodes, z, initial_guess) computeCorrelationLength_1D(nodes(:, 1), z, initial_guess);
                    elseif size(nodes, 2) == 2
                        corr_len_fun = @(nodes, z, initial_guess) computeCorrelationLength_2D(nodes(:, 1), nodes(:, 2), z, initial_guess);
                    end
                elseif strcmp(this.z_prior_interface.z_hyperparam_interface.z_type, 'transient vector')
                    nodes = this.z_prior_interface.z_hyperparam_interface.Load_Time_Node_Data();
                    corr_len_fun = @(t,z,initial_guess) computeCorrelationLength_transientVector(t,z,initial_guess);
                end
                this.z_pert_corr = zeros(num_perts,1);
                initial_guess = 0;
                for j = 1:num_perts
                    this.z_pert_corr(j) = corr_len_fun(nodes,this.z_pert(:,j),initial_guess);
                    initial_guess = this.z_pert_corr(j);
                end

                bin_width = 0.05;
                bin_centers = (min(this.z_pert_corr)-bin_width/2):bin_width:(max(this.z_pert_corr)+bin_width/2);
                num_bins = length(bin_centers);
                this.z_pert_corr_binned = cell(num_bins,1);
                this.delta_pert_mag_binned = cell(num_components, num_bins);
                for i = 1:num_bins
                    this.z_pert_corr_binned{i} = [];
                    for k = 1:num_components
                        this.delta_pert_mag_binned{k,i} = [];
                    end
                end
                for j = 1:num_perts
                    [~,i] = min(abs(bin_centers-this.z_pert_corr(j)));
                    this.z_pert_corr_binned{i} = [this.z_pert_corr_binned{i};this.z_pert_corr(j)];
                    for k = 1:num_components
                        this.delta_pert_mag_binned{k,i} = [this.delta_pert_mag_binned{k,i};this.delta_pert_mag{k,j}];
                    end
                end
                J = ~cellfun('isempty', this.z_pert_corr_binned);
                this.z_pert_corr_binned = this.z_pert_corr_binned(J);
                this.z_pert_corr_bin_means = cellfun(@(x)mean(x), this.z_pert_corr_binned);
                this.delta_pert_mag_binned = this.delta_pert_mag_binned(:,J);
            end
        end

        %%
        function [] = Plot_Transient_1D(this, nodes, t, u, range, name)
            n_t = length(t);
            n_y = length(nodes);
            u = reshape(u, n_y, n_t);
            m = round(n_t / 2);

            subplot(3, 1, 1);
            plot(nodes, u(:, 1), 'LineWidth', 3);
            xlabel('$x$','Interpreter','latex');
            ylabel(name,'Interpreter','latex');
            ylim(range);
            title(['Time t = ', num2str(t(1))]);
            set(gca, 'fontsize', 24);

            subplot(3, 1, 2);
            plot(nodes, u(:, m), 'LineWidth', 3);
            xlabel('$x$','Interpreter','latex');
            ylabel(name,'Interpreter','latex');
            ylim(range);
            title(['Time t = ', num2str(t(m))]);
            set(gca, 'fontsize', 24);

            subplot(3, 1, 3);
            plot(nodes, u(:, end), 'LineWidth', 3);
            xlabel('$x$','Interpreter','latex');
            ylabel(name,'Interpreter','latex');
            ylim(range);
            title(['Time t = ', num2str(t(end))]);
            set(gca, 'fontsize', 24);
        end

        function [] = Plot_Transient_2D(this, nodes, t, u, range, name)
            n_y = size(nodes,1);
            n_t = length(t);
            m = round(sqrt(n_y)*1.25);
            xl = linspace(min(nodes(:,1)), max(nodes(:,1)), m)';
            yl = linspace(min(nodes(:,2)), max(nodes(:,2)), m)';
            [X, Y] = meshgrid(xl, yl);
            u = reshape(u,n_y,n_t);
            count = 1;
            for k = round(linspace(1,n_t,3))
                subplot(3, 1, count);
                F = scatteredInterpolant(nodes(:,1), nodes(:,2), u(:,k));
                f = F(X, Y);
                surf(X, Y, f);
                xlabel('$x$','Interpreter','latex')
                ylabel('$y$','Interpreter','latex')
                view(2);
                shading interp;
                clim(range);
                colorbar();
                title([name,' at time ',num2str(t(k))],'Interpreter','latex');
                set(gca, 'fontsize', 24);
                count = count + 1;
            end

        end

        function [] = Plot_Stationary_1D(this, nodes, t, u, range, name)
            plot(nodes,u)
            xlabel('$x$','Interpreter','latex')
            ylabel(name,'Interpreter','latex')
            view(2);
            ylim(range);
            set(gca, 'fontsize', 24);
        end

        function [] = Plot_Stationary_2D(this, nodes, t, u, range, name)
            m = round(sqrt(size(nodes,1))*1.25);
            xl = linspace(min(nodes(:,1)), max(nodes(:,1)), m)';
            yl = linspace(min(nodes(:,2)), max(nodes(:,2)), m)';
            [X, Y] = meshgrid(xl, yl);
            F = scatteredInterpolant(nodes(:,1), nodes(:,2), u);
            f = F(X, Y);
            clf;
            surf(X, Y, f);
            xlabel('$x$','Interpreter','latex')
            ylabel('$y$','Interpreter','latex')
            view(2);
            shading interp;
            clim(range);
            colorbar();
            title(name,'Interpreter','latex');
            set(gca, 'fontsize', 24);

        end

        function [] = Plot_Vector(this, u, range, name)
            clf;
            plot(1:length(u), u, 'o', 'MarkerSize',10);
            xlabel('Component');
            ylabel(name,'Interpreter','latex');
            xlim([0,length(u)+1])
            range(1) = min(range(1),min(u) - .05*abs(min(u)));
            range(2) = max(range(2),max(u) + .05*abs(max(u)));
            ylim(range);
            set(gca, 'fontsize', 24);
        end

        function [] = Plot_Transient_Vector(this, t, u, range, name)
            clf;
            hold on;
            z = reshape(u, [], length(t));
            for k = 1:size(z, 1)
                plot(t, z(k, :), 'LineWidth', 3);
            end
            xlabel('Time');
            ylabel(name,'Interpreter','latex');
            ylim(range);
            set(gca, 'fontsize', 24);
        end

    end

end
