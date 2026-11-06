%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 2a: is more than one continuation step helpful for the posterior mean?
%
% The posterior-mean solve is run with fixed N_c = 1, 2, ..., then Newton polishing
% that stops once ||g_beta|| <= newton_gradient_tol (at most num_max_newton_steps
% steps), and compared with backtracking-only (a full step first, shrinking only on
% rejection) and the a priori heuristic.  For each setting we report
%   - the gradient norm right after the continuation (the Newton starting point),
%   - the gradient norm after each Newton step (to see the Newton convergence rate),
%   - cost = continuation steps (+ backtracks) + newton_cost * Newton steps taken,
% so we can see whether the extra continuation steps pay for themselves by saving Newton
% steps.  Needs only CI_Setup (no sample solves), so it is fast.

CI_Setup;

mean_fixed_sweep = [1, 2, 3, 4, 5];
num_max_newton_steps = 10;             % Newton polishing at t = 1: at most this many steps (each counted as work) ...
newton_gradient_tol = lowfi_grad_norm; % ... stopping once ||g_beta|| <= this
newton_cost = 0.5;                     % cost of a Newton step relative to a continuation step
mean_epsilon_nl = [0.1, 0.05, 0.01];   % sweep for the mean heuristic

fprintf('\n--- Experiment 2a: posterior mean, is N_c > 1 helpful? ---\n');
fprintf('Newton polishing: at most %d steps, tol ||g|| <= %.3e; cost = continuation steps + %.2f * Newton steps\n', ...
        num_max_newton_steps, newton_gradient_tol, newton_cost);
fprintf('%-28s | cont | ||g|| after cont | Newton | cost  | final ||g||  | ref err    | s\n', 'setting');

mean_fixed = [];
for jj = 1:length(mean_fixed_sweep)
    rr = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, mean_fixed_sweep(jj), [], mean_bt_options, num_max_newton_steps, newton_gradient_tol);
    rr.ref_err = Rel_Err(rr.z, mean_ref.z);
    rr.cost = rr.steps + rr.backtracks + newton_cost * rr.newton_steps;
    rr.label = sprintf('fixed N_c = %d', mean_fixed_sweep(jj));
    mean_fixed = [mean_fixed, rr]; %#ok<AGROW>
    Print_Row(rr);
end

mean_backtrack = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, 'backtrack', [], mean_bt_options, num_max_newton_steps, newton_gradient_tol);
mean_backtrack.ref_err = Rel_Err(mean_backtrack.z, mean_ref.z);
mean_backtrack.cost = mean_backtrack.steps + mean_backtrack.backtracks + newton_cost * mean_backtrack.newton_steps;
mean_backtrack.label = 'backtrack only (dt = 1 start)';
Print_Row(mean_backtrack);

mean_adapt = [];
for jj = 1:length(mean_epsilon_nl)
    rr = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, [], mean_epsilon_nl(jj), mean_bt_options, num_max_newton_steps, newton_gradient_tol);
    rr.ref_err = Rel_Err(rr.z, mean_ref.z);
    rr.cost = rr.steps + rr.backtracks + newton_cost * rr.newton_steps;
    rr.label = sprintf('adaptive eps_nl = %g', mean_epsilon_nl(jj));
    mean_adapt = [mean_adapt, rr]; %#ok<AGROW>
    Print_Row(rr);
    fprintf('%-28s   dt_prior %.3e | backtracks %d\n', '', rr.dt_prior, rr.backtracks);
end

% Newton convergence trajectories (gradient after the continuation, then after each Newton step).
fprintf('\nGradient norm after the continuation and after each Newton step:\n');
all_runs = [mean_fixed, mean_backtrack, mean_adapt];
for k = 1:length(all_runs)
    fprintf('  %-28s %s\n', all_runs(k).label, sprintf('%.2e  ', all_runs(k).newton_grad_norms));
end

[best_cost, jb] = min([mean_fixed.cost]);
fprintf('\nCheapest fixed N_c: %d (cost %.2f). Costs of N_c > 1 relative to it: %s\n', mean_fixed_sweep(jb), best_cost, ...
        sprintf('%.2fx  ', [mean_fixed.cost] / best_cost));
fprintf('Newton steps saved by each extra continuation step (N_c = 1 -> 2 -> ...): %s (each costs 1.0 vs %.2f per Newton step)\n', ...
        sprintf('%.2f  ', -diff([mean_fixed.newton_steps])), newton_cost);

results = struct();
results.settings = struct('num_data_points', num_data_points, 'mean_reference_steps', mean_reference_steps, ...
                          'mean_fixed_sweep', mean_fixed_sweep, 'mean_epsilon_nl', mean_epsilon_nl, ...
                          'epsilon_curv', epsilon_curv, 'num_max_newton_steps', num_max_newton_steps, ...
                          'newton_gradient_tol', newton_gradient_tol, 'newton_cost', newton_cost);
results.mean_ref = mean_ref;  results.mean_fixed = mean_fixed;  results.mean_backtrack = mean_backtrack;
results.mean_adapt = mean_adapt;
results_file = fullfile(ci_dir, 'ci-adaptive-mean-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Adaptive_Mean;

function Print_Row(rr)
    fprintf('%-28s | %4d |   %.3e     |  %3d   | %5.2f | %.3e  | %.3e  | %.2f\n', ...
            rr.label, rr.steps + rr.backtracks, rr.grad_after_continuation, rr.newton_steps, rr.cost, rr.grad_norm, rr.ref_err, rr.time);
end
