%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Prior_Continuation_Sampling < handle

    % Continuation-based sampling of the prior discrepancy-induced solution
    % update.  Each prior sample of the discrepancy,
    %
    %   delta_k(z) = data_shift + delta0_k + X_k * V' (z - z_opt),
    %
    % is introduced through a pseudo-time parameter t in [0, 1] and the
    % optimal solution is tracked from t = 0 (z = z_opt) to t = 1 via
    % Pseudo_Time_Continuation.  The realizations (delta0_k, X_k) are drawn
    % once in Generate_Prior_Samples; samples are then tracked independently
    % and may be processed in parallel (options.use_parallel_sampling).

    properties
        data_interface
        u_prior_interface
        z_prior_interface
        hessian_analysis
        opt_prob_interface
        u_opt
        z_opt
        r
        options

        delta0
        X
    end

    methods

        function this = MD_Prior_Continuation_Sampling(data_interface, u_prior_interface, z_prior_interface, hessian_analysis)
            arguments
                data_interface MD_Data_Interface
                u_prior_interface MD_u_Prior_Interface
                z_prior_interface MD_z_Prior_Interface
                hessian_analysis MD_Hessian_Analysis
            end
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.hessian_analysis = hessian_analysis;
            this.opt_prob_interface = hessian_analysis.opt_prob_interface;
            this.u_opt = data_interface.u_opt;
            this.z_opt = data_interface.z_opt;
            if isempty(hessian_analysis.evals)
                this.r = length(this.z_opt);
            else
                this.r = length(hessian_analysis.evals);
            end
            this.options = struct();
            this.options.num_steps = 3;
            this.options.use_parallel_sampling = false;
        end

        % Draw prior realizations (delta0_k, X_k), k = 1, ..., num_samples.
        function [] = Generate_Prior_Samples(this, num_samples)
            this.delta0 = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samples);
            this.X = cell(num_samples, 1);
            for k = 1:num_samples
                sampler = MD_Breve_Beta_Sampler(this.hessian_analysis, this.z_prior_interface, [], ...
                                                this.u_prior_interface, length(this.u_opt), 1e-10, true);
                this.X{k} = sampler.Get_Materialized_Matrix();
            end
        end

        function [u_ks, z_ks, beta_ks] = Prior_Update_Samples(this, num_samples, options)
            if nargin >= 3 && ~isempty(options)
                fields = fieldnames(options);
                for i = 1:length(fields)
                    this.options.(fields{i}) = options.(fields{i});
                end
            end
            if nargin < 2 || isempty(num_samples)
                num_samples = size(this.delta0, 2);
            end
            if isempty(this.delta0) || size(this.delta0, 2) < num_samples
                this.Generate_Prior_Samples(num_samples);
            end

            u_ks = zeros(length(this.u_opt), num_samples);
            z_ks = zeros(length(this.z_opt), num_samples);
            beta_ks = zeros(this.r, num_samples);

            use_parallel_sampling = this.options.use_parallel_sampling;
            if use_parallel_sampling && ~this.Parallel_Sampling_Available()
                warning('MD_Prior_Continuation_Sampling:ParallelUnavailable', ...
                    'Parallel sampling requested, but the Parallel Computing Toolbox is not available. Falling back to serial sampling.');
                use_parallel_sampling = false;
            end

            if use_parallel_sampling
                parfor k = 1:num_samples
                    [u_k, z_k, beta_k] = this.Prior_Update_Core(k);
                    u_ks(:, k) = u_k;
                    z_ks(:, k) = z_k;
                    beta_ks(:, k) = beta_k;
                end
            else
                for k = 1:num_samples
                    [u_ks(:, k), z_ks(:, k), beta_ks(:, k)] = this.Prior_Update_Core(k);
                end
            end
        end

        function [tf] = Parallel_Sampling_Available(this) %#ok<MANU>
            tf = exist('parallel.pool.DataQueue', 'class') == 8;
            if tf && exist('license', 'builtin') == 5
                tf = license('test', 'Distrib_Computing_Toolbox');
            end
        end

        function [u_k, z_k, beta_k] = Prior_Update_Core(this, sample_idx)
            % Use a separate operator object per sample so that parallel workers do not share state.
            sen_op_k = MD_Prior_Continuation_Sensitivity_Operators(this.hessian_analysis, this.z_opt, ...
                this.data_interface.data_shift, this.delta0, this.X);
            qn_prec = MD_Quasi_Newton_Preconditioner(this.hessian_analysis);
            pt_cont = Pseudo_Time_Continuation(zeros(this.r, 1), sen_op_k, qn_prec);
            theta_traj = MD_Discrepancy_Parameter_Trajectory(this.options.num_steps, sample_idx);
            beta_tmp = pt_cont.Pseudo_Time_Continuation_Forward_Euler(theta_traj);
            beta_k = beta_tmp(:, end);
            u_k = sen_op_k.current_u;
            z_k = sen_op_k.current_z;
        end

    end

end
