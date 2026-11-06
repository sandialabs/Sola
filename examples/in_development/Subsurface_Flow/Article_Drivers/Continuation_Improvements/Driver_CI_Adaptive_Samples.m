%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 2b: a priori time-step heuristic vs fixed steps for
% mean-to-sample continuation.  Uses the same fixed-step mean start as
% Experiment 1, isolating the sample heuristic.
%
% Every run finishes with Newton polishing at t = 1 that stops, per sample,
% once ||g_beta|| <= newton_gradient_tol (at most num_max_newton_steps steps).
% The endpoint gradient is therefore small by construction, and what we compare
% is how much work it takes to get there:
%   cost = continuation steps (+ backtracks) + newton_cost * Newton steps,
% together with the fraction of samples that reach the tolerance at all (a
% sample whose continuation ends outside Newton's basin does not).
% eps_curv is fixed at sample_epsilon_curv; the (eps_nl, eps_curv) grid is
% explored in Driver_CI_Curvature_Sweep.

ci_num_data_points = 1;                % this experiment uses 1 data point (other drivers use 3)
CI_Setup;
CI_Sample_Reference;
clear ci_num_data_points;

sample_fixed_sweep = [1, 2, 3, 4, 5];
num_max_newton_steps = 10;             % Newton polishing at t = 1: at most this many steps ...
newton_gradient_tol = lowfi_grad_norm; % ... stopping per sample once ||g_beta|| <= this
newton_cost = 0.5;                     % cost of a Newton step relative to a continuation step
sample_epsilon_nl = [0.05, 0.01, 0.005, 0.002];  % sweep for the sample heuristic
sample_epsilon_curv = 3;               % curvature limit (small values make eps_nl irrelevant)
show_sample_epsilon = 0.005;           % epsilon_nl used for single-run panels (nearest in sample_epsilon_nl)

sample_bt_options.epsilon_curv = sample_epsilon_curv;

fprintf('\n--- Experiment 2b: adaptive time steps, mean-to-sample ---\n');
fprintf('Newton polishing: at most %d steps, tol ||g|| <= %.3e; cost = continuation steps + %.2f * Newton steps\n', ...
        num_max_newton_steps, newton_gradient_tol, newton_cost);
fprintf('%-26s | cont mean/max | Newton mean/max | cost mean / p90 / max | converged | s/sample\n', 'setting');

sample_fixed = cell(1, length(sample_fixed_sweep));
for jj = 1:length(sample_fixed_sweep)
    sample_fixed{jj} = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', sample_fixed_sweep(jj), beta_bar, [], [], sample_indices, num_max_newton_steps, newton_gradient_tol);
    sample_fixed{jj}.ref_err = Rel_Col_Err(sample_fixed{jj}.z, sample_ref.z);
    sample_fixed{jj}.summary = Summarize_CI_Cost(sample_fixed{jj}, newton_cost, newton_gradient_tol, sprintf('fixed N_c = %d', sample_fixed_sweep(jj)));
    Print_Summary(sample_fixed{jj}.summary);
end
sample_adapt = cell(1, length(sample_epsilon_nl));
for jj = 1:length(sample_epsilon_nl)
    sample_adapt{jj} = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'adaptive_mts', [], beta_bar, sample_epsilon_nl(jj), sample_bt_options, sample_indices, num_max_newton_steps, newton_gradient_tol);
    sample_adapt{jj}.ref_err = Rel_Col_Err(sample_adapt{jj}.z, sample_ref.z);
    sample_adapt{jj}.summary = Summarize_CI_Cost(sample_adapt{jj}, newton_cost, newton_gradient_tol, sprintf('adaptive eps_nl = %g', sample_epsilon_nl(jj)));
    Print_Summary(sample_adapt{jj}.summary);
    fprintf('%-26s   dt_prior min/max %.3e/%.3e | backtracks %d | failures %d\n', '', ...
            min(sample_adapt{jj}.dt_prior), max(sample_adapt{jj}.dt_prior), sum(sample_adapt{jj}.backtracks), sum(sample_adapt{jj}.failures));
end

% Backtracking only: no a priori step size, full step first, shrink on rejection.
sample_backtrack = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'backtrack_mts', [], beta_bar, [], sample_bt_options, sample_indices, num_max_newton_steps, newton_gradient_tol);
sample_backtrack.ref_err = Rel_Col_Err(sample_backtrack.z, sample_ref.z);
sample_backtrack.summary = Summarize_CI_Cost(sample_backtrack, newton_cost, newton_gradient_tol, 'backtrack only (dt = 1 start)');
Print_Summary(sample_backtrack.summary);
fprintf('%-26s   backtracks %d | failures %d\n', '', sum(sample_backtrack.backtracks), sum(sample_backtrack.failures));

[~, show_idx] = min(abs(log(sample_epsilon_nl) - log(show_sample_epsilon)));

% Cheapest uniform choice in hindsight, for reference.
fixed_mean_costs = cellfun(@(s) s.summary.mean_cost, sample_fixed);
fixed_converged = cellfun(@(s) s.summary.converged_frac, sample_fixed);
ok_fixed = find(fixed_converged >= max(fixed_converged) - 1.e-12);
[best_cost, jb] = min(fixed_mean_costs(ok_fixed));
fprintf('\nBest uniform fixed N_c in hindsight (among fully converging): N_c = %d, mean cost %.2f\n', ...
        sample_fixed_sweep(ok_fixed(jb)), best_cost);
for jj = 1:length(sample_epsilon_nl)
    fprintf('  adaptive eps_nl = %-6g mean cost %.2f (%.2fx the best fixed), converged %.0f%%\n', sample_epsilon_nl(jj), ...
            sample_adapt{jj}.summary.mean_cost, sample_adapt{jj}.summary.mean_cost / best_cost, 100 * sample_adapt{jj}.summary.converged_frac);
end

results = struct();
results.settings = struct('num_data_points', num_data_points, 'num_samples', num_samples, ...
                          'fixed_mean_steps', fixed_mean_steps, 'reference_steps', reference_steps, ...
                          'sample_fixed_sweep', sample_fixed_sweep, 'sample_epsilon_nl', sample_epsilon_nl, ...
                          'sample_epsilon_curv', sample_epsilon_curv, 'show_idx', show_idx, ...
                          'num_max_newton_steps', num_max_newton_steps, 'newton_gradient_tol', newton_gradient_tol, ...
                          'newton_cost', newton_cost);
results.sample_ref = sample_ref;
results.sample_fixed = sample_fixed;  results.sample_adapt = sample_adapt;  results.sample_backtrack = sample_backtrack;
results_file = fullfile(ci_dir, 'ci-adaptive-samples-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Adaptive_Samples;

function Print_Summary(s)
    fprintf('%-26s |  %4.2f / %2d   |   %4.2f / %2d     |   %5.2f / %5.2f / %5.2f |   %5.1f%%  | %.2f\n', ...
            s.label, s.mean_cont, max(s.cont_steps), s.mean_newton, max(s.newton_steps), ...
            s.mean_cost, s.p90_cost, s.max_cost, 100 * s.converged_frac, s.time_per_sample);
end
