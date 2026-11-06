%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 2b: a priori time-step heuristic vs fixed steps for
% mean-to-sample continuation.  Uses the same fixed-step mean start as
% Experiment 1, isolating the sample heuristic.  Runs one fixed-step solve
% per entry of sample_fixed_sweep and one adaptive solve per
% sample_epsilon_nl, each over all samples.

ci_num_data_points = 1;                % this experiment uses 1 data point (other drivers use 3)
CI_Setup;
CI_Sample_Reference;
clear ci_num_data_points;

% Plots 7-9 use N_c = 1, 3, 5, 7 from this sweep.
sample_fixed_sweep = [1, 2, 3];
polish_M = 2;                          % extra Newton steps at t = 1 (counted as work)
sample_epsilon_nl = [0.05, 0.01, 0.005];  % sweep for the sample heuristic
show_sample_epsilon = 0.005;           % epsilon_nl used for single-run panels (nearest in sample_epsilon_nl)

fprintf('\n--- Experiment 2b: adaptive time steps, mean-to-sample ---\n');
sample_fixed = cell(1, length(sample_fixed_sweep));
for jj = 1:length(sample_fixed_sweep)
    sample_fixed{jj} = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', sample_fixed_sweep(jj), beta_bar, [], [], sample_indices, polish_M);
    sample_fixed{jj}.ref_err = Rel_Col_Err(sample_fixed{jj}.z, sample_ref.z);
    sample_fixed{jj}.work = (sample_fixed_sweep(jj) + polish_M) * ones(num_samples, 1);
    fprintf('fixed N = %2d | median ref err %.3e | median ||grad|| %.3e | %.2fs/sample\n', sample_fixed_sweep(jj), ...
            median(sample_fixed{jj}.ref_err, 'omitnan'), median(sample_fixed{jj}.grad_norm, 'omitnan'), sample_fixed{jj}.time / num_samples);
end
sample_adapt = cell(1, length(sample_epsilon_nl));
for jj = 1:length(sample_epsilon_nl)
    sample_adapt{jj} = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'adaptive_mts', [], beta_bar, sample_epsilon_nl(jj), sample_bt_options, sample_indices, polish_M);
    sample_adapt{jj}.ref_err = Rel_Col_Err(sample_adapt{jj}.z, sample_ref.z);
    fprintf('adaptive eps_nl = %g | dt_prior min/max %.3e/%.3e | work mean/max %.1f/%d | median ref err %.3e | median ||grad|| %.3e | %.2fs/sample | failures %d\n', ...
            sample_epsilon_nl(jj), min(sample_adapt{jj}.dt_prior), max(sample_adapt{jj}.dt_prior), ...
            mean(sample_adapt{jj}.work, 'omitnan'), max(sample_adapt{jj}.work), ...
            median(sample_adapt{jj}.ref_err, 'omitnan'), median(sample_adapt{jj}.grad_norm, 'omitnan'), ...
            sample_adapt{jj}.time / num_samples, sum(sample_adapt{jj}.failures));
end

[~, show_idx] = min(abs(log(sample_epsilon_nl) - log(show_sample_epsilon)));

results = struct();
results.settings = struct('num_data_points', num_data_points, 'num_samples', num_samples, ...
                          'fixed_mean_steps', fixed_mean_steps, 'reference_steps', reference_steps, ...
                          'sample_fixed_sweep', sample_fixed_sweep, 'sample_epsilon_nl', sample_epsilon_nl, ...
                          'epsilon_curv', epsilon_curv, 'show_idx', show_idx, 'polish_M', polish_M);
results.sample_ref = sample_ref;
results.sample_fixed = sample_fixed;  results.sample_adapt = sample_adapt;
results_file = fullfile(ci_dir, 'ci-adaptive-samples-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Adaptive_Samples;
