%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Prior_Visualization < handle

    properties
        prior_sampling
    end

    %% Constructor
    methods

        function this = MD_Prior_Visualization(prior_sampling)
            arguments
                prior_sampling MD_Prior_Sampling
            end
            this.prior_sampling = prior_sampling;
        end

    end

    %% Functions to generate visualizations
    methods

        function [] = Visualization_for_Prior_Discrepancy_at_z_pert(this, component_id)

            if ~strcmp(this.prior_sampling.z_prior_interface.z_hyperparam_interface.z_type, 'vector')

                I = this.prior_sampling.data_interface.Separate_State_Components(component_id);
                num_perts = length(this.prior_sampling.delta_samples_z_pert);
                nodes = this.prior_sampling.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
                nodes = nodes{component_id};
                spatial_dim = size(nodes, 2);
                is_transient = this.prior_sampling.u_prior_interface.u_hyperparam_interface.is_transient;
                if is_transient
                    t = this.prior_sampling.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
                    t_z = this.prior_sampling.z_prior_interface.z_hyperparam_interface.Load_Time_Node_Data();
                else
                    t = [];
                    t_z = [];
                end

                z_range = [min(this.prior_sampling.z_pert(:)), max(this.prior_sampling.z_pert(:))];

                if strcmp(this.prior_sampling.z_prior_interface.z_hyperparam_interface.z_type, 'transient vector')
                    generate_z_plot = @(nodes, t, u, z_range, name) this.Plot_Transient_Vector(t, u, z_range, name);
                elseif strcmp(this.prior_sampling.z_prior_interface.z_hyperparam_interface.z_type, 'spatial field')
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
                generate_z_plot(nodes, t_z, this.prior_sampling.data_interface.Z(:, 1), z_range, '$z_1$');

                delta_fig = figure;
                c_range = [min(this.prior_sampling.data_interface.D(I, 1) + this.prior_sampling.data_interface.data_shift(I)), max(this.prior_sampling.data_interface.D(I, 1) + this.prior_sampling.data_interface.data_shift(I))];
                generate_delta_plot(nodes, t, this.prior_sampling.data_interface.D(I, 1) + this.prior_sampling.data_interface.data_shift(I), c_range, '$d_1$');

                margin = max(this.prior_sampling.z_pert_corr) * .05;
                xmin = min(this.prior_sampling.z_pert_corr) - margin;
                xmax = max(this.prior_sampling.z_pert_corr) + margin;
                yrange = [0, 0];
                for k = 1:num_perts
                    yrange(1) = min(yrange(1), min(this.prior_sampling.delta_pert_mag{component_id, k}));
                    yrange(2) = max(yrange(2), max(this.prior_sampling.delta_pert_mag{component_id, k}));
                end
                ymin = yrange(1) * .9;
                ymax = yrange(2) * 1.1;
                colors = lines(length(this.prior_sampling.z_pert_corr_binned));
                scatter_fig = figure;
                hold on;
                for k = 1:length(this.prior_sampling.z_pert_corr_binned)
                    violinplot(this.prior_sampling.z_pert_corr_bin_means(k), this.prior_sampling.delta_pert_mag_binned{component_id, k}, 'DensityDirection', 'positive', 'DensityWidth', .05, 'DensityScale', 'width', 'FaceColor', 'black');
                    plot(this.prior_sampling.z_pert_corr_bin_means(k), this.prior_sampling.delta_pert_mag_binned{component_id, k}, 'o', 'Color', colors(k, :), 'MarkerSize', 8);
                end
                xlabel('$\Delta z$ Correlation Length', 'Interpreter', 'latex');
                ylabel('$\delta(\Delta z)$ Magnitude', 'Interpreter', 'latex');
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

                        bin = find(abs(this.prior_sampling.z_pert_corr_bin_means - x) < .01);
                        J = find(this.prior_sampling.z_pert_corr == this.prior_sampling.z_pert_corr_binned{bin}(1));
                        for k = 2:length(this.prior_sampling.z_pert_corr_binned{bin})
                            J = [J; find(this.prior_sampling.z_pert_corr == this.prior_sampling.z_pert_corr_binned{bin}(k))];
                        end
                        i = J(1);
                        [val, j] = min(abs(this.prior_sampling.delta_pert_mag{component_id, i} - y));
                        for k = 2:length(J)
                            [valk, jk] = min(abs(this.prior_sampling.delta_pert_mag{component_id, J(k)} - y));
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
                        plot(this.prior_sampling.z_pert_corr_bin_means(bin), this.prior_sampling.delta_pert_mag{component_id, i}(j), 'x', 'MarkerSize', 20, 'color', 'black');

                        figure(z_fig);
                        generate_z_plot(nodes, t_z, this.prior_sampling.z_pert(:, i), z_range, '$\Delta z_k$');

                        figure(delta_fig);
                        c_range = [min(this.prior_sampling.delta_samples_z_pert{i}(I, j)), max(this.prior_sampling.delta_samples_z_pert{i}(I, j))];
                        generate_delta_plot(nodes, t, this.prior_sampling.delta_samples_z_pert{i}(I, j), c_range, '$\delta(\Delta z_k,\theta_i)$');

                        figure(delta_fig);
                        figure(z_fig);
                        figure(scatter_fig);
                    end
                end
                disp('Concluded interactive visualization');

            else

                I = this.prior_sampling.data_interface.Separate_State_Components(component_id);
                num_perts = length(this.prior_sampling.delta_samples_z_pert);
                nodes = this.prior_sampling.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
                nodes = nodes{component_id};
                spatial_dim = size(nodes, 2);
                is_transient = this.prior_sampling.u_prior_interface.u_hyperparam_interface.is_transient;
                if is_transient
                    t = this.prior_sampling.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
                else
                    t = [];
                    t_z = [];
                end

                z_range = [min(this.prior_sampling.z_pert(:)), max(this.prior_sampling.z_pert(:))];
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
                generate_z_plot(nodes, t_z, this.prior_sampling.data_interface.Z(:, 1), z_range, '$z_1$');

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
                c_range = [min(this.prior_sampling.data_interface.D(I, 1) + this.prior_sampling.data_interface.data_shift(I)), max(this.prior_sampling.data_interface.D(I, 1) + this.prior_sampling.data_interface.data_shift(I))];
                generate_delta_plot(nodes, t, this.prior_sampling.data_interface.D(I, 1) + this.prior_sampling.data_interface.data_shift(I), c_range, '$d_1$');

                xmin = 0;
                xmax = this.prior_sampling.z_prior_interface.n_z + 1;
                yrange = [0, 0];
                for k = 1:num_perts
                    yrange(1) = min(yrange(1), min(this.prior_sampling.delta_pert_mag{component_id, k}));
                    yrange(2) = max(yrange(2), max(this.prior_sampling.delta_pert_mag{component_id, k}));
                end
                ymin = yrange(1) * .9;
                ymax = yrange(2) * 1.1;
                colors = lines(this.prior_sampling.z_prior_interface.n_z);
                scatter_fig = figure;
                hold on;
                for k = 1:this.prior_sampling.z_prior_interface.n_z
                    violinplot(k, this.prior_sampling.delta_pert_mag{component_id, k}, 'DensityDirection', 'positive', 'DensityWidth', .05, 'DensityScale', 'width', 'FaceColor', 'black');
                    plot(k, this.prior_sampling.delta_pert_mag{component_id, k}, 'o', 'Color', colors(k, :), 'MarkerSize', 8);
                end
                xlabel('$z$ Component', 'Interpreter', 'latex');
                ylabel('$\delta(\Delta z)$ Magnitude', 'Interpreter', 'latex');
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

                        i = find(abs((1:this.prior_sampling.z_prior_interface.n_z)' - x) < .1);
                        [~, j] = min(abs(this.prior_sampling.delta_pert_mag{component_id, i} - y));

                        figure(scatter_fig);
                        if length(scatter_fig.Children(1).Children) > num_perts
                            delete(scatter_fig.Children(1).Children(1));
                        end
                        plot(i, this.prior_sampling.delta_pert_mag{component_id, i}(j), 'x', 'MarkerSize', 20, 'color', 'black');

                        figure(z_fig);
                        generate_z_plot(nodes, t_z, this.prior_sampling.z_pert(:, i), z_range, '$\Delta z_k$');

                        figure(delta_fig);
                        c_range = [min(this.prior_sampling.delta_samples_z_pert{i}(I, j)), max(this.prior_sampling.delta_samples_z_pert{i}(I, j))];
                        generate_delta_plot(nodes, t, this.prior_sampling.delta_samples_z_pert{i}(I, j), c_range, '$\delta(\Delta z_k,\theta_i)$');

                        figure(delta_fig);
                        figure(z_fig);
                        figure(scatter_fig);
                    end
                end
                disp('Concluded interactive visualization');

            end

        end

        function [] = Visualization_for_Prior_Discrepancy_at_z_opt(this, component_id)

            I = this.prior_sampling.data_interface.Separate_State_Components(component_id);
            nodes = this.prior_sampling.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
            nodes = nodes{component_id};
            spatial_dim = size(nodes, 2);
            is_transient = this.prior_sampling.u_prior_interface.u_hyperparam_interface.is_transient;
            if is_transient
                t = this.prior_sampling.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
            else
                t = [];
            end
            range = [min(min(this.prior_sampling.delta_samples_z_opt(I, :))), max(max(this.prior_sampling.delta_samples_z_opt(I, :)))];

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
            generate_plot(nodes, t, this.prior_sampling.data_interface.D(I, 1) + this.prior_sampling.data_interface.data_shift(I), '$d_1$');

            sample_fig = figure;
            generate_plot(nodes, t, this.prior_sampling.data_interface.D(I, 1) + this.prior_sampling.data_interface.data_shift(I), '$d_1$');

            xmin = min([this.prior_sampling.d1_corr{component_id}; this.prior_sampling.delta_corr{component_id}]) * .9;
            xmax = max([this.prior_sampling.d1_corr{component_id}; this.prior_sampling.delta_corr{component_id}]) * 1.1;
            ymin = min([this.prior_sampling.d1_mag{component_id}; this.prior_sampling.delta_mag{component_id}]) * .9;
            ymax = max([this.prior_sampling.d1_mag{component_id}; this.prior_sampling.delta_mag{component_id}]) * 1.1;
            scatter_fig = figure;
            hold on;
            plot(this.prior_sampling.delta_corr{component_id}, this.prior_sampling.delta_mag{component_id}, 'o');
            plot(this.prior_sampling.d1_corr{component_id}, this.prior_sampling.d1_mag{component_id}, '*');
            xlabel('Correlation Length');
            ylabel('Magnitude');
            legend({'Prior Samples', 'Discrepancy Data'});
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
                    [~, i] = min(vecnorm([this.prior_sampling.delta_corr{component_id}, this.prior_sampling.delta_mag{component_id}] - [x, y], 2, 2));
                    figure(sample_fig);
                    generate_plot(nodes, t, this.prior_sampling.delta_samples_z_opt(I, i), '$\delta(\tilde{z},\theta_i)$');

                    figure(scatter_fig);
                    if length(scatter_fig.Children(1).Children) > 1
                        delete(scatter_fig.Children(1).Children(1));
                    end
                    plot(this.prior_sampling.delta_corr{component_id}(i), this.prior_sampling.delta_mag{component_id}(i), 'x', 'MarkerSize', 20, 'color', 'black');
                    legend({'Prior Samples', 'Discrepancy Data'});

                    figure(sample_fig);
                    figure(data_fig);
                    figure(scatter_fig);
                end
            end
            disp('Concluded interactive visualization');

        end

        function [] = Visualization_for_Prior_Time_Evolution(this, component_id, interactive_viz)

            delta_min = .95 * min([this.prior_sampling.delta_z_opt_time_evol{component_id}(:); this.prior_sampling.data_time_evol{component_id}]);
            delta_max = 1.05 * max([this.prior_sampling.delta_z_opt_time_evol{component_id}(:); this.prior_sampling.data_time_evol{component_id}]);
            t = this.prior_sampling.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();

            sample_fig = figure;
            hold on;
            plot(t, this.prior_sampling.delta_z_opt_time_evol{component_id}, 'LineWidth', 3, 'Color', [.9, .9, .9]);
            plot(t, this.prior_sampling.data_time_evol{component_id}, 'LineWidth', 3, 'Color', 'red');
            xlabel('Time');
            ylabel('Discrepancy Spatial Norm');
            ylim([delta_min, delta_max]);
            set(gca, 'fontsize', 24);

            data_fig = figure;
            hold on;
            plot(t, this.prior_sampling.data_time_evol{component_id}, 'LineWidth', 3, 'Color', 'red');
            xlabel('Time');
            ylabel('Discrepancy Spatial Norm');
            ylim([delta_min, delta_max]);
            set(gca, 'fontsize', 24);

            xmid = median(this.prior_sampling.temporal_corr_len{component_id});
            ymid = median(this.prior_sampling.temporal_mag{component_id});
            xmin = min(this.prior_sampling.temporal_corr_len{component_id}) - 0.3 * xmid;
            xmax = max(this.prior_sampling.temporal_corr_len{component_id}) + 0.3 * xmid;
            ymin = min(this.prior_sampling.temporal_mag{component_id}) - 0.3 * ymid;
            ymax = max(this.prior_sampling.temporal_mag{component_id}) + 0.3 * ymid;
            scatter_fig = figure;
            hold on;
            plot(this.prior_sampling.temporal_corr_len{component_id}, this.prior_sampling.temporal_mag{component_id}, 'o');
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

                    [~, i] = min(vecnorm([this.prior_sampling.temporal_corr_len{component_id}, this.prior_sampling.temporal_mag{component_id}] - [x, y], 2, 2));
                    plot(data_fig.CurrentAxes, t, this.prior_sampling.delta_z_opt_time_evol{component_id}(i, :), 'LineWidth', 3, 'Color', colors(current_selection, :));
                    plot(scatter_fig.CurrentAxes, this.prior_sampling.temporal_corr_len{component_id}(i), this.prior_sampling.temporal_mag{component_id}(i), 'x', 'MarkerSize', 15, 'Color', colors(current_selection, :));

                    current_selection = current_selection + 1;
                    figure(sample_fig);
                    figure(data_fig);
                    figure(scatter_fig);
                end
            end
            disp('Concluded interactive visualization');
        end

    end

    %% Plotting functions
    methods

        function [] = Plot_Transient_1D(this, nodes, t, u, range, name)
            n_t = length(t);
            n_y = length(nodes);
            u = reshape(u, n_y, n_t);
            m = round(n_t / 2);

            subplot(3, 1, 1);
            plot(nodes, u(:, 1), 'LineWidth', 3);
            xlabel('$x$', 'Interpreter', 'latex');
            ylabel(name, 'Interpreter', 'latex');
            ylim(range);
            title(['Time t = ', num2str(t(1))]);
            set(gca, 'fontsize', 24);

            subplot(3, 1, 2);
            plot(nodes, u(:, m), 'LineWidth', 3);
            xlabel('$x$', 'Interpreter', 'latex');
            ylabel(name, 'Interpreter', 'latex');
            ylim(range);
            title(['Time t = ', num2str(t(m))]);
            set(gca, 'fontsize', 24);

            subplot(3, 1, 3);
            plot(nodes, u(:, end), 'LineWidth', 3);
            xlabel('$x$', 'Interpreter', 'latex');
            ylabel(name, 'Interpreter', 'latex');
            ylim(range);
            title(['Time t = ', num2str(t(end))]);
            set(gca, 'fontsize', 24);
        end

        function [] = Plot_Transient_2D(this, nodes, t, u, range, name)
            n_y = size(nodes, 1);
            n_t = length(t);
            m = round(sqrt(n_y) * 1.25);
            xl = linspace(min(nodes(:, 1)), max(nodes(:, 1)), m)';
            yl = linspace(min(nodes(:, 2)), max(nodes(:, 2)), m)';
            [X, Y] = meshgrid(xl, yl);
            u = reshape(u, n_y, n_t);
            count = 1;
            for k = round(linspace(1, n_t, 3))
                subplot(3, 1, count);
                F = scatteredInterpolant(nodes(:, 1), nodes(:, 2), u(:, k));
                f = F(X, Y);
                surf(X, Y, f);
                xlabel('$x$', 'Interpreter', 'latex');
                ylabel('$y$', 'Interpreter', 'latex');
                view(2);
                shading interp;
                clim(range);
                colorbar();
                title([name, ' at time ', num2str(t(k))], 'Interpreter', 'latex');
                set(gca, 'fontsize', 24);
                count = count + 1;
            end

        end

        function [] = Plot_Stationary_1D(this, nodes, t, u, range, name)
            plot(nodes, u);
            xlabel('$x$', 'Interpreter', 'latex');
            ylabel(name, 'Interpreter', 'latex');
            view(2);
            ylim(range);
            set(gca, 'fontsize', 24);
        end

        function [] = Plot_Stationary_2D(this, nodes, t, u, range, name)
            m = round(sqrt(size(nodes, 1)) * 1.25);
            xl = linspace(min(nodes(:, 1)), max(nodes(:, 1)), m)';
            yl = linspace(min(nodes(:, 2)), max(nodes(:, 2)), m)';
            [X, Y] = meshgrid(xl, yl);
            F = scatteredInterpolant(nodes(:, 1), nodes(:, 2), u);
            f = F(X, Y);
            clf;
            surf(X, Y, f);
            xlabel('$x$', 'Interpreter', 'latex');
            ylabel('$y$', 'Interpreter', 'latex');
            view(2);
            shading interp;
            clim(range);
            colorbar();
            title(name, 'Interpreter', 'latex');
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
            ylabel(name, 'Interpreter', 'latex');
            ylim(range);
            set(gca, 'fontsize', 24);
        end

        function [] = Plot_Vector(this, u, range, name)
            clf;
            plot(1:length(u), u, 'o', 'MarkerSize', 10);
            xlabel('Component');
            ylabel(name, 'Interpreter', 'latex');
            xlim([0, length(u) + 1]);
            range(1) = min(range(1), min(u) - .05 * abs(min(u)));
            range(2) = max(range(2), max(u) + .05 * abs(max(u)));
            ylim(range);
            set(gca, 'fontsize', 24);
        end

    end

end
