%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Hybrid_Continuation_Update < handle

    % Hybrid posterior sampling:
    %
    %   z_k = F_Nc(theta_bar) + grad_theta F(theta_bar) (theta_k - theta_bar)
    %       = z_bar - V H_beta^{-1} B_k,
    %
    % where z_bar = z_opt + V beta_bar is obtained by pseudo-time continuation
    % on the posterior mean discrepancy, H_beta is the beta-space Hessian at
    % (beta_bar, theta_bar), and B_k is the mixed derivative applied to the
    % discrepancy perturbation theta_k - theta_bar (see
    % MD_Continuation_Sensitivity_Operators.Apply_B_hybrid).
    %
    % One continuation run is performed (for the mean); each sample then
    % costs a single preconditioned CG solve with the Hessian at the mean.
    % The quasi-Newton preconditioner built during the mean continuation is
    % reused for these solves.

    properties
        post_sampling
        hessian_analysis
        opt_prob_interface
        u_opt
        z_opt
        num_continuation_steps
        r
        discard_cache

        sen_op
        qn_prec
        pt_cont
        theta_traj
        u_bar
        z_bar
        beta_bar
    end

    methods

        function this = MD_Hybrid_Continuation_Update(post_sampling, hessian_analysis, num_continuation_steps, discard_cache)
            arguments
                post_sampling MD_Posterior_Sampling
                hessian_analysis MD_Hessian_Analysis
                num_continuation_steps (1, 1) {mustBeNumeric}
                discard_cache = true
            end
            this.post_sampling = post_sampling;
            this.hessian_analysis = hessian_analysis;
            this.opt_prob_interface = hessian_analysis.opt_prob_interface;
            this.u_opt = post_sampling.data_interface.u_opt;
            this.z_opt = post_sampling.data_interface.z_opt;
            if isempty(hessian_analysis.evals)
                this.r = length(this.z_opt);
            else
                this.r = length(hessian_analysis.evals);
            end
            this.num_continuation_steps = num_continuation_steps;
            this.discard_cache = discard_cache;
        end

        function [u, z, beta] = Posterior_Update_Mean(this)
            % Continuation for the posterior mean (cached after first call).
            if isempty(this.beta_bar)
                this.sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.discard_cache);
                this.qn_prec = MD_Quasi_Newton_Preconditioner(this.hessian_analysis);
                this.pt_cont = Pseudo_Time_Continuation(zeros(this.r, 1), this.sen_op, this.qn_prec);
                this.theta_traj = MD_Discrepancy_Parameter_Trajectory(this.num_continuation_steps, 0);
                beta_traj = this.pt_cont.Pseudo_Time_Continuation_Forward_Euler(this.theta_traj);

                this.beta_bar = beta_traj(:, end);
                this.z_bar = this.z_opt + this.hessian_analysis.Apply_V(this.beta_bar);
                this.u_bar = this.opt_prob_interface.State_Solve(this.z_bar);
            end
            u = this.u_bar;
            z = this.z_bar;
            beta = this.beta_bar;
        end

        function [u_ks, z_ks, beta_ks] = Posterior_Update_Samples(this, solve_state)
            % solve_state: if true (default), u_k = S(z_k) for each sample;
            % otherwise u_ks is returned empty to save state solves.
            arguments
                this
                solve_state (1, 1) logical = true
            end
            this.Posterior_Update_Mean();

            num_samples = this.post_sampling.post_data.num_samples;
            N = this.theta_traj.Get_Number_of_Timesteps();
            z_ks = zeros(length(this.z_opt), num_samples);
            beta_ks = zeros(this.r, num_samples);
            if solve_state
                u_ks = zeros(length(this.u_opt), num_samples);
            else
                u_ks = [];
            end

            for sample_idx = 1:num_samples
                B_k = this.sen_op.Apply_B_hybrid(this.beta_bar, this.theta_traj, N, sample_idx);
                dbeta = this.pt_cont.Apply_Inverse_Hessian(B_k, this.beta_bar, this.theta_traj, N);
                beta_k = this.beta_bar - dbeta;

                beta_ks(:, sample_idx) = beta_k;
                z_ks(:, sample_idx) = this.z_opt + this.hessian_analysis.Apply_V(beta_k);
                if solve_state
                    u_ks(:, sample_idx) = this.opt_prob_interface.State_Solve(z_ks(:, sample_idx));
                end
            end
        end

    end

end
