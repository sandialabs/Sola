%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 5: a problem where one continuation step is NOT enough.
%
% On the standard problem N_c = 1 plus Newton always lands in Newton's basin,
% so the time-step heuristic has nothing to do.  Here the mean-to-sample
% increment is inflated by a factor lambda: the continuation runs from the
% posterior mean (t = 0) to mean + lambda * (sample - mean) (t = 1), see
% MD_Discrepancy_Extrapolated_Trajectory.  lambda = 1 is the ordinary problem;
% increasing lambda makes the endpoint further from the start (and more
% nonlinear), so there is a lambda* beyond which N_c = 1 followed by Newton
% fails to converge.  For each lambda we compare
%   - fixed N_c (1, 2, 4, 8, 16),
%   - backtracking only (full step first, shrink on rejection),
%   - the a priori heuristic + backtracking,
% scored by the fraction of samples that reach the Newton tolerance and by cost
%   cost = continuation steps (+ backtracks, failed trial steps) + newton_cost * Newton steps.
% The question is whether the adaptive rules converge (nearly) everywhere at a
% cost close to the cheapest fixed N_c that works for each lambda, without being
% told which N_c to use.  Newton is capped low (num_max_newton_steps) so a bad
% continuation endpoint actually shows up as non-convergence.

ci_num_data_points = 1;
CI_Setup;
clear ci_num_data_points;

lambdas = [1, 2, 4, 8, 16];        % increment inflation factors
fixed_steps = [1, 2, 4];        % fixed N_c
num_max_newton_steps = 5;              % keep modest: a bad endpoint must fail
newton_gradient_tol = max(lowfi_grad_norm, 1.e-9);   % fixed floor (see findings, item 7)
newton_cost = 0.5;
adaptive_epsilon_nl = 0.08;
adaptive_epsilon_curv = Inf;            % gentle curvature limit (section 4.3)
include_backtrack_only = false;    % also run the backtracking-only method
use_parallel_sampling = false;     % run the samples of each method in parallel (parfor; pool or Threads)
adaptive_max_steps = 60;               % the default cap (15) would mask the heuristic

bt_options = sample_bt_options;
bt_options.epsilon_curv = adaptive_epsilon_curv;
bt_options.max_steps = adaptive_max_steps;

method_labels = arrayfun(@(n) sprintf('fixed N_c = %d', n), fixed_steps, 'UniformOutput', false);
if include_backtrack_only
    method_labels{end + 1} = 'backtrack only';
end
method_labels{end + 1} = 'adaptive heuristic';
nM = length(method_labels);
nL = length(lambdas);
summaries = cell(nM, nL);
runs = cell(nM, nL);

fprintf('\n--- Experiment 5: extrapolated mean-to-sample increments ---\n');
fprintf('Newton: at most %d steps, tol ||g|| <= %.2e; cost = steps (+ backtracks) + %.2f * Newton steps\n', ...
        num_max_newton_steps, newton_gradient_tol, newton_cost);
for ll = 1:nL
    lam = lambdas(ll);
    fprintf('\nlambda = %g\n', lam);
    fprintf('%-26s | cont mean/max | Newton mean/max | cost mean / p90 / max | converged | s/sample\n', 'setting');
    for mm = 1:nM
        if mm <= length(fixed_steps)
            out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', fixed_steps(mm), beta_bar, [], [], ...
                              sample_indices, num_max_newton_steps, newton_gradient_tol, lam, use_parallel_sampling);
        elseif strcmp(method_labels{mm}, 'backtrack only')
            out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'backtrack_mts', [], beta_bar, [], bt_options, ...
                              sample_indices, num_max_newton_steps, newton_gradient_tol, lam, use_parallel_sampling);
        else
            out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'adaptive_mts', [], beta_bar, adaptive_epsilon_nl, bt_options, ...
                              sample_indices, num_max_newton_steps, newton_gradient_tol, lam, use_parallel_sampling);
        end
        s = Summarize_CI_Cost(out, newton_cost, newton_gradient_tol, method_labels{mm});
        summaries{mm, ll} = s;
        runs{mm, ll} = out;
        fprintf('%-26s |  %4.2f / %2d   |   %4.2f / %2d     |   %5.2f / %5.2f / %5.2f |   %5.1f%%  | %.2f\n', ...
                s.label, s.mean_cont, max(s.cont_steps), s.mean_newton, max(s.newton_steps), ...
                s.mean_cost, s.p90_cost, s.max_cost, 100 * s.converged_frac, s.time_per_sample);
    end
end

% Where does N_c = 1 stop working?
conv_nc1 = cellfun(@(s) s.converged_frac, summaries(1, :));
fprintf('\nFraction converged with N_c = 1 (+Newton) per lambda:\n');
fprintf('  lambda = %-6g: %.0f%%\n', [lambdas; 100 * conv_nc1]);
if all(conv_nc1 >= 1 - 1.e-12)
    fprintf('N_c = 1 converged everywhere; extend lambdas upward.\n');
end

results = struct();
results.settings = struct('num_data_points', num_data_points, 'num_samples', num_samples, ...
                          'lambdas', lambdas, 'fixed_steps', fixed_steps, 'method_labels', {method_labels}, ...
                          'num_max_newton_steps', num_max_newton_steps, 'newton_gradient_tol', newton_gradient_tol, ...
                          'newton_cost', newton_cost, 'adaptive_epsilon_nl', adaptive_epsilon_nl, ...
                          'adaptive_epsilon_curv', adaptive_epsilon_curv, 'adaptive_max_steps', adaptive_max_steps, ...
                          'include_backtrack_only', include_backtrack_only, 'use_parallel_sampling', use_parallel_sampling);
results.summaries = summaries;
results.runs = runs;
results_file = fullfile(ci_dir, 'ci-extrapolated-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Extrapolated;
