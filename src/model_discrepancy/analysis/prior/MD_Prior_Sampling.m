classdef MD_Prior_Sampling < handle

    properties
        data_interface
        u_prior_interface
        z_prior_interface
        u_opt
        z_opt

        delta_samples_z_opt
        delta_samples_z_pert
        z_pert
        z_pert_evals
        delta_mag
        delta_corr
        delta_pert_mag

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

        function [] = Visualization_for_Prior_Discrepancy_at_z_pert(this,component_id)

            I = this.data_interface.Separate_State_Components(component_id);
            num_perts = length(this.delta_samples_z_pert);
            nodes = this.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
            nodes = nodes{component_id};
            spatial_dim = size(nodes,2);
            is_transient = this.u_prior_interface.u_hyperparam_interface.is_transient;
            if is_transient
                t = this.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
            end
            delta_range = [min(this.delta_samples_z_pert{1}(:)) , max(this.delta_samples_z_pert{1}(:))];
            z_range = [min(this.z_pert(:)),max(this.z_pert(:))];
            for k = 2:num_perts
                delta_range(1) = min(delta_range(1),min(this.delta_samples_z_pert{k}(:)));
                delta_range(2) = max(delta_range(2),max(this.delta_samples_z_pert{k}(:)));
            end

            if strcmp(this.z_prior_interface.z_hyperparam_interface.z_type,'transient vector')
                generate_z_plot = @(t,u) this.Plot_Transient_Vector(t,u,z_range);
            end

            if is_transient && (spatial_dim == 1)
                generate_delta_plot = @(nodes,t,u) this.Plot_Transient_1D(nodes,t,u,delta_range);
            elseif is_transient && (spatial_dim == 2)
                generate_delta_plot = @(nodes,t,u) this.Plot_Transient_2D(nodes,t,u,delta_range);
            elseif ~is_transient && (spatial_dim == 1)
                generate_delta_plot = @(nodes,t,u) this.Plot_Stationary_1D(nodes,t,u,delta_range);
            elseif ~is_transient && (spatial_dim == 2)
                generate_delta_plot = @(nodes,t,u) this.Plot_Stationary_2D(nodes,t,u,delta_range);
            end

            z_fig = figure;
            generate_z_plot(t,this.data_interface.Z(:,1))

            delta_fig = figure;
            generate_delta_plot(nodes,t,this.data_interface.D(I,1))

            xmin = 0.05;
            xmax = 1.05;
            yrange = [0,0];
            for k = 1:num_perts
                yrange(1) = min(yrange(1),min(this.delta_pert_mag{component_id,k}));
                yrange(2) = max(yrange(2),max(this.delta_pert_mag{component_id,k}));
            end
            ymin = yrange(1) * .9;
            ymax = yrange(2) * 1.1;
            num_samps = size(this.delta_pert_mag{component_id,1},1);
            scatter_fig = figure;
            hold on
            for k = 1:num_perts
                plot(this.z_pert_evals(k)*ones(num_samps,1),this.delta_pert_mag{component_id,k},'o')
            end
            xlabel('z Perturbation Correlation Length')
            ylabel('Delta Perturbation Magnitude')
            xlim([xmin,xmax])
            ylim([ymin,ymax])
            set(gca, 'fontsize', 24);


            user_continue = true;
            bnd_check = @(x,y) (x < xmin) || (x > xmax) || (y < ymin) || (y > ymax);

            disp('Please select point number for rendering. Select a point outside of the plotting domain to terminate')
            while user_continue 

                figure(scatter_fig);
                [x,y] = ginput(1);

                if bnd_check(x,y)
                    user_continue = false;
                else

                    [~,i] = min(abs(this.z_pert_evals-x));
                    [val,j] = min(abs(this.delta_pert_mag{i}-y));
                    J = find(abs(this.z_pert_evals-this.z_pert_evals(i))<1.e-8);
                    for k = 1:length(J)
                        [valk,jk] = min(abs(this.delta_pert_mag{J(k)}-y));
                        if valk < val
                            val = valk;
                            j = jk;
                            i = J(k);
                        end
                    end

                    figure(scatter_fig);
                    if length(scatter_fig.Children(1).Children) > num_perts
                        delete(scatter_fig.Children(1).Children(1))
                    end
                    plot(this.z_pert_evals(i),this.delta_pert_mag{i}(j),'x','MarkerSize',20,'color','black')

                    figure(z_fig);
                    generate_z_plot(t,this.z_pert(:,i));

                    figure(delta_fig);
                    generate_delta_plot(nodes,t,this.delta_samples_z_pert{i}(I,j));

                    figure(delta_fig);
                    figure(z_fig);
                    figure(scatter_fig);
                end
            end
            disp('Concluded interactive visualization')


        end

        function [] = Visualization_for_Prior_Discrepancy_at_z_opt(this,component_id)

            I = this.data_interface.Separate_State_Components(component_id);
            nodes = this.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
            nodes = nodes{component_id};
            spatial_dim = size(nodes,2);
            is_transient = this.u_prior_interface.u_hyperparam_interface.is_transient;
            if is_transient
                t = this.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
            end
            range = [min(this.delta_samples_z_opt(:)) , max(this.delta_samples_z_opt(:))];

            if is_transient && (spatial_dim == 1)
                generate_plot = @(nodes,t,u) this.Plot_Transient_1D(nodes,t,u,range);
            elseif is_transient && (spatial_dim == 2)
                generate_plot = @(nodes,t,u) this.Plot_Transient_2D(nodes,t,u,range);
            elseif ~is_transient && (spatial_dim == 1)
                generate_plot = @(nodes,t,u) this.Plot_Stationary_1D(nodes,t,u,range);
            elseif ~is_transient && (spatial_dim == 2)
                generate_plot = @(nodes,t,u) this.Plot_Stationary_2D(nodes,t,u,range);
            end

            data_fig = figure;
            generate_plot(nodes,t,this.data_interface.D(I,1))

            sample_fig = figure;
            generate_plot(nodes,t,this.data_interface.D(I,1))

            xmin = min(this.delta_corr{component_id}) * .9;
            xmax = max(this.delta_corr{component_id}) * 1.1;
            ymin = min(this.delta_mag{component_id}) * .9;
            ymax = max(this.delta_mag{component_id}) * 1.1;
            scatter_fig = figure;
            hold on
            plot(this.delta_corr{component_id},this.delta_mag{component_id},'o')
            xlabel('Correlation Length')
            ylabel('Magnitude')
            xlim([xmin,xmax])
            ylim([ymin,ymax])
            set(gca, 'fontsize', 24);


            user_continue = true;
            bnd_check = @(x,y) (x < xmin) || (x > xmax) || (y < ymin) || (y > ymax);

            disp('Please select point number for rendering. Select a point outside of the plotting domain to terminate')
            while user_continue 

                figure(scatter_fig);
                [x,y] = ginput(1);

                if bnd_check(x,y)
                    user_continue = false;
                else
                    [~,i] = min(vecnorm([this.delta_corr{component_id},this.delta_mag{component_id}]-[x,y],2,2));
                    figure(sample_fig);
                    generate_plot(nodes,t,this.delta_samples_z_opt(:,i));

                    figure(sample_fig);
                    figure(data_fig);
                    figure(scatter_fig);
                end
            end
            disp('Concluded interactive visualization')


        end

        function [] = Visualization_for_Prior_Time_Evolution(this,component_id)
            
            num_components = length(this.delta_mag);
            if isempty(this.delta_z_opt_time_evol)
                this.delta_z_opt_time_evol = cell(num_components,1);
            end
            
            num_samps = size(this.delta_samples_z_opt,2);
            I = this.data_interface.Separate_State_Components(component_id);
            t = this.u_prior_interface.u_hyperparam_interface.Load_Time_Node_Data();
            n_t = length(t);
            n_y = length(this.delta_samples_z_opt(I,1)) / n_t;
            this.delta_z_opt_time_evol{component_id} = zeros(num_samps,n_t);

            for i = 1:num_samps
                di = reshape(this.delta_samples_z_opt(I, i), n_y, n_t);
                this.delta_z_opt_time_evol{component_id}(i,:) = sqrt(diag(di' * this.u_prior_interface.Apply_M_u(di)));
            end

            d = reshape(this.data_interface.D(I,1), n_y, n_t);
            this.data_time_evol = sqrt(diag(d' * this.u_prior_interface.Apply_M_u(d)));

            temporal_corr_len = zeros(num_samps,1);
            initial_guess = 0;
            for i = 1:num_samps
                temporal_corr_len(i) = computeCorrelationLength_1D(t,this.delta_z_opt_time_evol{component_id}(i,:),initial_guess);
                initial_guess = temporal_corr_len(i);
            end
            temporal_mag = max(abs(this.delta_z_opt_time_evol{component_id}),[],2);

            delta_min = .95*min(this.delta_z_opt_time_evol{component_id}(:));
            delta_max = 1.05*max(temporal_mag);

            sample_fig = figure;
            hold on
            plot(t,this.delta_z_opt_time_evol{component_id},'LineWidth',3,'Color',[.9,.9,.9])
            plot(t,this.data_time_evol,'LineWidth',3,'Color','red')
            xlabel('Time')
            ylabel('Discrepancy Spatial Norm')
            ylim([delta_min,delta_max])
            set(gca, 'fontsize', 24);

            data_fig = figure;
            hold on
            plot(t,this.data_time_evol,'LineWidth',3,'Color','red')
            xlabel('Time')
            ylabel('Discrepancy Spatial Norm')
            ylim([delta_min,delta_max])
            set(gca, 'fontsize', 24);

            xmid = median(temporal_corr_len);
            ymid = median(temporal_mag);
            xmin = min(temporal_corr_len) - 0.3 * xmid;
            xmax = max(temporal_corr_len) + 0.3 * xmid;
            ymin = min(temporal_mag) - 0.3 * ymid;
            ymax = max(temporal_mag) + 0.3 * ymid;
            scatter_fig = figure;
            hold on
            plot(temporal_corr_len,temporal_mag,'o')
            xlabel('Temporal Correlation Length')
            ylabel('Magnitude')
            xlim([xmin,xmax])
            ylim([ymin,ymax])
            set(gca, 'fontsize', 24);

            user_continue = true;
            bnd_check = @(x,y) (x < xmin) || (x > xmax) || (y < ymin) || (y > ymax);
            colors = lines(20);
            current_selection = 1;
            
            disp('Please select point number for rendering. Select a point outside of the plotting domain to terminate')
            while user_continue && (current_selection <= 20)

                [x,y] = ginput(1);

                if bnd_check(x,y)
                    user_continue = false;
                else

                    [~,i] = min(vecnorm([temporal_corr_len,temporal_mag]-[x,y],2,2));
                    plot(data_fig.CurrentAxes,t,this.delta_z_opt_time_evol{component_id}(i,:),'LineWidth',3,'Color',colors(current_selection,:));
                    plot(scatter_fig.CurrentAxes,temporal_corr_len(i),temporal_mag(i),'x','MarkerSize',15,'Color',colors(current_selection,:));

                    current_selection = current_selection + 1;
                    figure(sample_fig);
                    figure(data_fig);
                    figure(scatter_fig);
                end
            end
            disp('Concluded interactive visualization')
        end

        function [] = Generate_Prior_Discrepancy_Sample_Data(this, num_samps)

            this.delta_samples_z_opt = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps)  + this.data_interface.data_shift;
            this.Compute_Delta_z_opt_Metrics();

            this.Compute_z_pert();
            this.Compute_Delta_z_pert_Metrics();
        end

        function [] = Compute_z_pert(this)
            e = this.z_prior_interface.determine_z_hyperparams.Compute_Eigenvalues(this.z_prior_interface);
            num_perts_init = length(find(1./e>.1));
            if strcmp(this.z_prior_interface.z_hyperparam_interface.z_type,'transient vector')
                num_perts_init = this.z_prior_interface.num_controls * num_perts_init;
            end
            oversampling = min(10,length(this.data_interface.z_opt)-num_perts_init);
            num_subspace_iters = 1;
            E_z_inv_gsvd = E_z_Inv_GSVD(this.z_prior_interface, this.z_opt);
            [this.z_pert, ~, this.z_pert_evals] = E_z_inv_gsvd.Compute_GSVD(num_perts_init, oversampling, num_subspace_iters);
            while this.z_pert_evals(end) > .1
                num_perts_init = 2 * num_perts_init;
                oversampling = min(10,length(this.data_interface.z_opt)-num_perts_init);
                [this.z_pert, ~, this.z_pert_evals] = E_z_inv_gsvd.Compute_GSVD(num_perts_init, oversampling, num_subspace_iters);
            end
            I = find(this.z_pert_evals>.1);
            this.z_pert = this.z_pert(:,I);
            this.z_pert_evals = this.z_pert_evals(I);
            scaling = .3 * sqrt(this.z_opt' * this.z_prior_interface.Apply_M_z(this.z_opt));
            this.z_pert = scaling * this.z_pert;

            num_perts = length(this.z_pert_evals);
            num_samps = size(this.delta_samples_z_opt,2);
            this.delta_samples_z_pert = cell(num_perts, 1);
            for k = 1:num_perts
                this.delta_samples_z_pert{k} = scaling * sqrt(this.z_prior_interface.alpha_z) * this.z_pert_evals(k) * this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
            end
        end

        function [] = Compute_Delta_z_opt_Metrics(this)

            nodes = this.u_prior_interface.u_hyperparam_interface.Load_Spatial_Node_Data();
            num_components = length(nodes);

            this.delta_mag = cell(num_components,1);
            this.delta_corr = cell(num_components,1);
            for k = 1:num_components
                I = this.data_interface.Separate_State_Components(k);
                n_y = size(nodes{k}, 1);
                n_t = size(this.delta_samples_z_opt(I,:),1) / n_y;

                num_samps = size(this.delta_samples_z_opt(I,:),2);
                this.delta_mag{k} = max(abs(this.delta_samples_z_opt(I,:)))';

                initial_guess = 0;
                correlation_lengths = zeros(num_samps, n_t);
                if size(nodes, 2) == 1
                    corr_len_fun = @(nodes,d,initial_guess) computeCorrelationLength_1D(nodes(:, 1), d, initial_guess);
                elseif size(nodes,2) == 2
                    corr_len_fun = @(nodes,d,initial_guess) computeCorrelationLength_2D(nodes(:, 1), nodes(:,2), d, initial_guess);
                end

                for i = 1:num_samps
                    di = reshape(this.delta_samples_z_opt(I, i), n_y, n_t);
                    for j = 1:n_t
                        correlation_lengths(i, j) = corr_len_fun(nodes{k}, di(:, j), initial_guess);
                        initial_guess = correlation_lengths(i, j);
                    end
                    initial_guess = correlation_lengths(i, 1);
                end
                this.delta_corr{k} = mean(correlation_lengths,2);
            end
        end

        function Compute_Delta_z_pert_Metrics(this)
            num_components = length(this.delta_mag);
            num_perts = size(this.z_pert,2);
            this.delta_pert_mag = cell(num_components,num_perts);
            for k = 1:num_components
                I = this.data_interface.Separate_State_Components(k);
                for j = 1:num_perts
                    this.delta_pert_mag{k,j} = max(abs(this.delta_samples_z_pert{j}(I,:)))';
                end
            end
        end

        function [] = Plot_Transient_1D(this,nodes,t,u,range)
            n_t = length(t);
            n_y = length(nodes);
            u = reshape(u,n_y,n_t);
            m = round(n_t/2);

            subplot(3,1,1)
            plot(nodes,u(:,1),'LineWidth',3)
            xlabel('x')
            ylabel('Discrepancy')
            ylim(range)
            title(['Time t = ',num2str(t(1))])
            set(gca, 'fontsize', 24);

            subplot(3,1,2)
            plot(nodes,u(:,m),'LineWidth',3)
            xlabel('x')
            ylabel('Discrepancy')
            ylim(range)
            title(['Time t = ',num2str(t(m))])
            set(gca, 'fontsize', 24);

            subplot(3,1,3)
            plot(nodes,u(:,end),'LineWidth',3)
            xlabel('x')
            ylabel('Discrepancy')
            ylim(range)
            title(['Time t = ',num2str(t(end))])
            set(gca, 'fontsize', 24);
        end

        function [] = Plot_Transient_2D(this,nodes,t,u,range)
           
        end

        function [] = Plot_Stationary_1D(this,nodes,t,u,range)

        end

        function [] = Plot_Stationary_2D(this,nodes,t,u,range)

        end

        function [] = Plot_Transient_Vector(this,t,u,range)
            clf;
            hold on;
            z = reshape(u,[],length(t));
            for k = 1:size(z,1)
                plot(t,z(k,:),'LineWidth',3)
            end
            xlabel('Time')
            ylabel('Controller')
            ylim(range)
            set(gca, 'fontsize', 24);
        end

    end

end

