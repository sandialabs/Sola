%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Demonstrates the practical benefit of two continuation improvements.
%
%   Experiment 1: mean-to-sample vs nominal-to-sample continuation.
%       Nominal-to-sample starts every sample path at beta = 0 and must
%       traverse the whole nonlinear update.  Mean-to-sample starts at the
%       posterior-mean update beta_bar (computed once) and only continues
%       the sample-minus-mean discrepancy increment.  Claim: for the same
%       number of steps it is far more accurate, i.e. it reaches a given
%       accuracy with fewer steps / less wall time, even after paying for
%       the mean solve.
%
%   Experiment 2: a priori time-step heuristic (with backtracking safeguard)
%       vs fixed uniform steps.  Claim: for the same amount of work
%       (accepted steps + backtracked steps) the adaptive grid is more
%       accurate, and it spends its work where it is needed (per sample)
%       instead of using one hand-tuned step count for everything.
%
% Accuracy is measured against a high-step reference solve of the same
% endpoint problem (relative z distance) and by the endpoint optimality
% residual ||grad_beta J||, which does not depend on any reference.
% All solves share one MD_Continuation_Sensitivity_Operators instance with
% materialized breve samples, so every method sees identical sample
% realizations.

OED_Setup;
oed_results = load(fullfile(article_dir, 'oed-results.mat'));
num_data_points = 3;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

% ------------------------
% User-adjustable settings
% ------------------------
num_samples = 20;
export_figures = false;      % passed to Plot_Continuation_Improvements

% Experiment 1
fixed_mean_steps = 7;                  % mean solve used to start mean-to-sample
step_counts = [1, 2, 3, 5, 7, 10];     % N_c sweep for both sample methods
reference_steps = 25;                  % mean-to-sample steps for the sample reference
show_N_c = 2;                          % N_c used for per-sample scatter and bands (in step_counts)

% Experiment 2
mean_fixed_sweep = [1, 2, 3, 5, 7, 10, 15];
mean_epsilon_nl = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002];      % sweep for the mean heuristic
sample_fixed_sweep = [1, 2, 3, 5, 7, 10];
sample_epsilon_nl = [0.01, 0.005, 0.002, 0.001, 0.0005]; % sweep for the sample heuristic
show_sample_epsilon = 0.002;           % epsilon_nl used for single-run panels (nearest in sample_epsilon_nl)
epsilon_curv = 1;       
sample_max_steps = 15;                 % small epsilon_nl needs more than the default 15 steps
mean_reference_steps = 40;

fprintf('\n============================================================\n');
fprintf('Continuation improvements: mean-to-sample and adaptive steps\n');
fprintf('============================================================\n');

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

if isempty(md_hessian_analysis.evals)
    r = length(data_interface.z_opt);
else
    r = length(md_hessian_analysis.evals);
end
z_opt = data_interface.z_opt;

% discard_cache = false and materialized breve samples: all methods see the
% same sample realizations and trial evaluations are side-effect free.
sen_op = MD_Continuation_Sensitivity_Operators(md_post_sampling, md_hessian_analysis, false, true);
cont_defaults = MD_Continuation_Update(md_post_sampling, md_hessian_analysis).Default_Options();
mean_bt_options = cont_defaults.mean.time_step.options;
sample_bt_options = cont_defaults.samples.time_step.options;
mean_bt_options.epsilon_curv = epsilon_curv;
sample_bt_options.epsilon_curv = epsilon_curv;
sample_bt_options.max_steps = sample_max_steps;
sample_indices = 1:num_samples;

%% References
fprintf('\nMean reference (%d fixed steps)...\n', mean_reference_steps);
mean_ref = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, mean_reference_steps, [], mean_bt_options);

%% Experiment 1: mean-to-sample vs nominal-to-sample
fprintf('\n--- Experiment 1: mean-to-sample vs nominal-to-sample ---\n');
mean_start = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, fixed_mean_steps, [], mean_bt_options);
beta_bar = mean_start.beta;
fprintf('Mean solve (%d steps): %.2fs, ||grad|| = %.3e\n', fixed_mean_steps, mean_start.time, mean_start.grad_norm);

fprintf('Sample reference: mean-to-sample, N_c = %d (%d samples)...\n', reference_steps, num_samples);
sample_ref = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', reference_steps, beta_bar, [], [], sample_indices);
fprintf('  reference median/max endpoint ||grad|| = %.3e / %.3e\n', ...
        median(sample_ref.grad_norm, 'omitnan'), max(sample_ref.grad_norm, [], 'omitnan'));

num_N = length(step_counts);
exp1 = struct('N_c', num2cell(step_counts), 'nominal', [], 'mts', []);
for jj = 1:num_N
    N_c = step_counts(jj);
    exp1(jj).nominal = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'nominal', N_c, [], [], [], sample_indices);
    exp1(jj).mts = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', N_c, beta_bar, [], [], sample_indices);
    exp1(jj).nominal.ref_err = Rel_Col_Err(exp1(jj).nominal.z, sample_ref.z);
    exp1(jj).mts.ref_err = Rel_Col_Err(exp1(jj).mts.z, sample_ref.z);
    fprintf('N_c = %2d | median ref err: nominal %.2e, MTS %.2e | median ||grad||: nominal %.2e, MTS %.2e | time/sample: %.2fs vs %.2fs\n', ...
            N_c, median(exp1(jj).nominal.ref_err, 'omitnan'), median(exp1(jj).mts.ref_err, 'omitnan'), ...
            median(exp1(jj).nominal.grad_norm, 'omitnan'), median(exp1(jj).mts.grad_norm, 'omitnan'), ...
            exp1(jj).nominal.time / num_samples, exp1(jj).mts.time / num_samples);
end

%% Experiment 2a: adaptive time steps for the posterior mean
fprintf('\n--- Experiment 2a: adaptive time steps, posterior mean ---\n');
mean_fixed = [];
for jj = 1:length(mean_fixed_sweep)
    rr = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, mean_fixed_sweep(jj), [], mean_bt_options);
    rr.ref_err = Rel_Err(rr.z, mean_ref.z);
    mean_fixed = [mean_fixed, rr]; %#ok<AGROW>
    fprintf('fixed N = %2d | work %2d | ref err %.3e | ||grad|| %.3e | %.2fs\n', mean_fixed_sweep(jj), ...
            mean_fixed(jj).work, mean_fixed(jj).ref_err, mean_fixed(jj).grad_norm, mean_fixed(jj).time);
end
mean_adapt = [];
for jj = 1:length(mean_epsilon_nl)
    rr = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, [], mean_epsilon_nl(jj), mean_bt_options);
    rr.ref_err = Rel_Err(rr.z, mean_ref.z);
    mean_adapt = [mean_adapt, rr]; %#ok<AGROW>
    fprintf('adaptive eps_nl = %g | dt_prior %.3e | steps %2d, backtracks %d | work %2d | ref err %.3e | ||grad|| %.3e | %.2fs\n', ...
            mean_epsilon_nl(jj), rr.dt_prior, rr.steps, rr.backtracks, rr.work, rr.ref_err, rr.grad_norm, rr.time);
end

%% Experiment 2b: adaptive time steps for mean-to-sample continuation
fprintf('\n--- Experiment 2b: adaptive time steps, mean-to-sample ---\n');
% Same fixed-step mean start as Experiment 1, isolating the sample heuristic.
sample_fixed = cell(1, length(sample_fixed_sweep));
for jj = 1:length(sample_fixed_sweep)
    sample_fixed{jj} = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', sample_fixed_sweep(jj), beta_bar, [], [], sample_indices);
    sample_fixed{jj}.ref_err = Rel_Col_Err(sample_fixed{jj}.z, sample_ref.z);
    sample_fixed{jj}.work = sample_fixed_sweep(jj) * ones(num_samples, 1);
    fprintf('fixed N = %2d | median ref err %.3e | median ||grad|| %.3e | %.2fs/sample\n', sample_fixed_sweep(jj), ...
            median(sample_fixed{jj}.ref_err, 'omitnan'), median(sample_fixed{jj}.grad_norm, 'omitnan'), sample_fixed{jj}.time / num_samples);
end
sample_adapt = cell(1, length(sample_epsilon_nl));
for jj = 1:length(sample_epsilon_nl)
    sample_adapt{jj} = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'adaptive_mts', [], beta_bar, sample_epsilon_nl(jj), sample_bt_options, sample_indices);
    sample_adapt{jj}.ref_err = Rel_Col_Err(sample_adapt{jj}.z, sample_ref.z);
    fprintf('adaptive eps_nl = %g | dt_prior min/max %.3e/%.3e | work mean/max %.1f/%d | median ref err %.3e | median ||grad|| %.3e | %.2fs/sample | failures %d\n', ...
            sample_epsilon_nl(jj), min(sample_adapt{jj}.dt_prior), max(sample_adapt{jj}.dt_prior), ...
            mean(sample_adapt{jj}.work, 'omitnan'), max(sample_adapt{jj}.work), ...
            median(sample_adapt{jj}.ref_err, 'omitnan'), median(sample_adapt{jj}.grad_norm, 'omitnan'), ...
            sample_adapt{jj}.time / num_samples, sum(sample_adapt{jj}.failures));
end

[~, show_idx] = min(abs(log(sample_epsilon_nl) - log(show_sample_epsilon)));

%% Save
results = struct();
results.settings = struct('num_samples', num_samples, 'fixed_mean_steps', fixed_mean_steps, ...
                          'step_counts', step_counts, 'reference_steps', reference_steps, ...
                          'sample_epsilon_nl', sample_epsilon_nl, 'mean_epsilon_nl', mean_epsilon_nl, ...
                          'sample_fixed_sweep', sample_fixed_sweep, 'show_N_c', show_N_c, ...
                          'epsilon_curv', epsilon_curv, 'show_idx', show_idx);
results.z_opt = z_opt;  results.z_hifi = z_hifi;
results.mean_ref = mean_ref;  results.mean_start = mean_start;  results.sample_ref = sample_ref;
results.exp1 = exp1;
results.mean_fixed = mean_fixed;  results.mean_adapt = mean_adapt;
results.sample_fixed = sample_fixed;  results.sample_adapt = sample_adapt;
results_file = fullfile(article_dir, 'continuation-improvements-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

% Plotting lives in its own function that reads the saved file, so results
% survive any graphics problem and figures can be regenerated without rerunning.
% Plot_Continuation_Improvements(results_file, export_figures);
Plot_Continuation_Improvements_Basic;

% -------------------------------------------------------------------------
% Local helper functions
% -------------------------------------------------------------------------

function out = Run_Mean(sen_op, hessian_analysis, z_opt, r, num_steps, epsilon_nl, bt_options)
    % Posterior-mean continuation: fixed uniform steps if epsilon_nl is
    % empty, otherwise the a priori heuristic with backtracking.
    qn_prec = MD_Quasi_Newton_Preconditioner(hessian_analysis);
    pt = Pseudo_Time_Continuation(zeros(r, 1), sen_op, qn_prec);
    tic;
    if isempty(epsilon_nl)
        theta = MD_Discrepancy_Parameter_Trajectory(num_steps, 0);
        beta_traj = pt.Pseudo_Time_Continuation_Forward_Euler(theta);
        out.steps = num_steps;
        out.backtracks = 0;
        out.time_grid = linspace(0, 1, num_steps + 1);
        out.dt_prior = 1 / num_steps;
    else
        theta = MD_Discrepancy_Parameter_Trajectory(1, 0);
        bt_options.time_step_callback = @(beta, th, ti, beta_dot) ...
            sen_op.Mean_A_Priori_Time_Step_Heuristic(beta, th, ti, beta_dot, epsilon_nl, bt_options);
        [beta_traj, ~, out.time_grid, step_diags, prior_diag] = pt.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta, [], bt_options);
        out.dt_prior = prior_diag.dt;
        out.steps = length(step_diags);
        out.backtracks = sum(cellfun(@(d) d.num_backtracks, step_diags));
    end
    out.time = toc;
    out.work = out.steps + out.backtracks;
    out.beta = beta_traj(:, end);
    out.z = z_opt + hessian_analysis.Apply_V(out.beta);
    grad = sen_op.Gradient(out.beta, theta, theta.Get_Number_of_Timesteps());
    out.grad_norm = norm(grad);
end

function out = Run_Samples(sen_op, hessian_analysis, z_opt, r, method, N_c, beta_bar, epsilon_nl, bt_options, sample_indices)
    % method: 'nominal'      fixed steps from beta = 0
    %         'fixed_mts'    fixed steps from beta_bar
    %         'adaptive_mts' heuristic + backtracking from beta_bar
    num_samples = length(sample_indices);
    out.z = nan(length(z_opt), num_samples);
    out.grad_norm = nan(num_samples, 1);
    out.steps = nan(num_samples, 1);
    out.backtracks = zeros(num_samples, 1);
    out.dt_prior = nan(num_samples, 1);
    out.failures = false(num_samples, 1);

    tic;
    for ii = 1:num_samples
        idx = sample_indices(ii);
        try
            qn_prec = MD_Quasi_Newton_Preconditioner(hessian_analysis);
            switch method
                case 'nominal'
                    pt = Pseudo_Time_Continuation(zeros(r, 1), sen_op, qn_prec);
                    theta = MD_Discrepancy_Parameter_Trajectory(N_c, idx);
                    beta_traj = pt.Pseudo_Time_Continuation_Forward_Euler(theta);
                    out.steps(ii) = N_c;
                    out.dt_prior(ii) = 1 / N_c;
                case 'fixed_mts'
                    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
                    theta = MD_Discrepancy_Mean_To_Sample_Trajectory(N_c, idx);
                    beta_traj = pt.Pseudo_Time_Continuation_Forward_Euler(theta);
                    out.steps(ii) = N_c;
                    out.dt_prior(ii) = 1 / N_c;
                case 'adaptive_mts'
                    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
                    theta = MD_Discrepancy_Mean_To_Sample_Trajectory(1, idx);
                    opts = bt_options;
                    opts.time_step_callback = @(beta, th, ti, beta_dot) ...
                        sen_op.Sample_A_Priori_Time_Step_Heuristic(beta, th, ti, beta_dot, idx, epsilon_nl, bt_options);
                    [beta_traj, ~, ~, step_diags, prior_diag] = pt.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta, [], opts);
                    out.dt_prior(ii) = prior_diag.dt;
                    out.steps(ii) = length(step_diags);
                    out.backtracks(ii) = sum(cellfun(@(d) d.num_backtracks, step_diags));
                otherwise
                    error('Unknown method "%s".', method);
            end
            beta_k = beta_traj(:, end);
            out.z(:, ii) = z_opt + hessian_analysis.Apply_V(beta_k);

            % All methods are scored on the same endpoint (t = 1) sample objective.
            theta_end = MD_Discrepancy_Parameter_Trajectory(1, idx);
            out.grad_norm(ii) = norm(sen_op.Gradient(beta_k, theta_end, 1));
        catch ME
            out.failures(ii) = true;
            fprintf('  %s sample %d failed: %s\n', method, idx, ME.message);
        end
    end
    out.time = toc;
    out.work = out.steps + out.backtracks;
end

function err = Rel_Err(a, b)
    err = norm(a - b) / max(1, norm(b));
end

function errs = Rel_Col_Err(A, B)
    errs = nan(size(A, 2), 1);
    for ii = 1:size(A, 2)
        if all(isfinite(A(:, ii)))
            errs(ii) = Rel_Err(A(:, ii), B(:, ii));
        end
    end
end

function lh = Band(samples, ptile)
    lh = [prctile(samples, 50 - ptile / 2, 2), prctile(samples, 50 + ptile / 2, 2)];
end
