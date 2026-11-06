%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 2a: a priori time-step heuristic vs fixed steps for the
% posterior mean.  Needs only CI_Setup (no sample solves), so it is fast.

CI_Setup;

mean_fixed_sweep = [1, 2, 5];
num_max_newton_steps = 10;             % Newton polishing at t = 1: at most this many steps (each counted as work) ...
newton_gradient_tol = lowfi_grad_norm; % ... stopping once ||g_beta|| <= this
mean_epsilon_nl = [0.1, 0.05, 0.01];   % sweep for the mean heuristic

fprintf('\n--- Experiment 2a: adaptive time steps, posterior mean ---\n');
mean_fixed = [];
for jj = 1:length(mean_fixed_sweep)
    rr = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, mean_fixed_sweep(jj), [], mean_bt_options, num_max_newton_steps, newton_gradient_tol);
    rr.ref_err = Rel_Err(rr.z, mean_ref.z);
    mean_fixed = [mean_fixed, rr]; %#ok<AGROW>
    fprintf('fixed N = %2d | Newton %d | work %2d | ref err %.3e | ||grad|| %.3e | %.2fs\n', mean_fixed_sweep(jj), rr.newton_steps, ...
            rr.work, rr.ref_err, rr.grad_norm, rr.time);
end
mean_adapt = [];
for jj = 1:length(mean_epsilon_nl)
    rr = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, [], mean_epsilon_nl(jj), mean_bt_options, num_max_newton_steps, newton_gradient_tol);
    rr.ref_err = Rel_Err(rr.z, mean_ref.z);
    mean_adapt = [mean_adapt, rr]; %#ok<AGROW>
    fprintf('adaptive eps_nl = %g | dt_prior %.3e | steps %2d, backtracks %d | Newton %d | work %2d | ref err %.3e | ||grad|| %.3e | %.2fs\n', ...
            mean_epsilon_nl(jj), rr.dt_prior, rr.steps, rr.backtracks, rr.newton_steps, rr.work, rr.ref_err, rr.grad_norm, rr.time);
end

results = struct();
results.settings = struct('num_data_points', num_data_points, 'mean_reference_steps', mean_reference_steps, ...
                          'mean_fixed_sweep', mean_fixed_sweep, 'mean_epsilon_nl', mean_epsilon_nl, ...
                          'epsilon_curv', epsilon_curv, 'num_max_newton_steps', num_max_newton_steps, ...
                          'newton_gradient_tol', newton_gradient_tol);
results.mean_ref = mean_ref;  results.mean_fixed = mean_fixed;  results.mean_adapt = mean_adapt;
results_file = fullfile(ci_dir, 'ci-adaptive-mean-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Adaptive_Mean;
