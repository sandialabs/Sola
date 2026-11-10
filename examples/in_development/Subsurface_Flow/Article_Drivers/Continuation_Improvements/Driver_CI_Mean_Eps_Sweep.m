%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 6: each limit of the a priori time-step heuristic on its own,
% for the posterior mean.  The eps_nl sweep also records the gradient norm after
% 0, 1 and 2 Newton polishing steps (one run with 2 Newton steps gives all
% three: the iterates are prefixes); the eps_curv sweep and the fixed N_c
% reference use no Newton polishing.
%
%   (a) eps_curv = Inf, sweep eps_nl   (state-map nonlinearity limit only)
%   (b) eps_nl   = Inf, sweep eps_curv (Hessian-change / curvature limit only)
%
% The endpoint gradient norm ||grad_beta J|| after the continuation alone shows
% how well each limit controls the path error.  For reference the same quantity
% is computed for fixed N_c steps (no Newton), and the low-fidelity optimizer's
% gradient norm is the usual "good enough" level.  Needs only CI_Setup (no
% sample solves), so it is fast.
%
% Note: the heuristic is evaluated once at t = 0 from the Euler direction, and
% the mean's curvature estimate eta may be exactly zero (findings, item 6); then
% the curvature limit does nothing and sweep (b) is flat.

CI_Setup;

eps_nl_values = [0.1, 0.03, 0.01, 0.003];
eps_curv_values = [1, 0.1, 0.05, 0.01];
fixed_steps = [1, 2, 3, 5, 7, 10];

mean_bt_options.max_steps = max(mean_bt_options.max_steps, 60);   % do not let the default cap mask small eps

num_polish = 3;                        % Newton polishing steps recorded for the eps_nl sweep
run_one = @(eps_nl, eps_curv, num_newton) Run_Mean_Eps(sen_op, md_hessian_analysis, z_opt, r, eps_nl, eps_curv, mean_bt_options, num_newton);

fprintf('\n--- Experiment 6: mean, heuristic limits, no Newton polishing ---\n');
fprintf('Low-fidelity optimizer ||grad|| = %.3e\n', lowfi_grad_norm);

fprintf('\n(a) eps_curv = Inf, sweep eps_nl\n');
fprintf('%-10s | steps | backtracks | dt_prior | ||grad|| after 0 / 1 / 2 Newton\n', 'eps_nl');
nl_runs = cell(size(eps_nl_values));
for ii = 1:length(eps_nl_values)
    nl_runs{ii} = run_one(eps_nl_values(ii), Inf, num_polish);
    fprintf('%-10g | %5d | %10d | %.3e | %.3e / %.3e / %.3e\n', eps_nl_values(ii), nl_runs{ii}.steps, nl_runs{ii}.backtracks, ...
            nl_runs{ii}.dt_prior, nl_runs{ii}.newton_grad_norms(1:num_polish + 1));
end

fprintf('\n(b) eps_nl = Inf, sweep eps_curv\n');
fprintf('%-10s | steps | backtracks | dt_prior | ||grad||\n', 'eps_curv');
curv_runs = cell(size(eps_curv_values));
for ii = 1:length(eps_curv_values)
    curv_runs{ii} = run_one(Inf, eps_curv_values(ii), 0);
    fprintf('%-10g | %5d | %10d | %.3e | %.3e\n', eps_curv_values(ii), curv_runs{ii}.steps, curv_runs{ii}.backtracks, ...
            curv_runs{ii}.dt_prior, curv_runs{ii}.grad_norm);
end

fprintf('\nFixed N_c reference (no Newton)\n');
fixed_runs = cell(size(fixed_steps));
for ii = 1:length(fixed_steps)
    fixed_runs{ii} = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, fixed_steps(ii), [], mean_bt_options);
    fprintf('N_c = %2d | ||grad|| %.3e\n', fixed_steps(ii), fixed_runs{ii}.grad_norm);
end

results = struct();
results.settings = struct('num_data_points', num_data_points, 'eps_nl_values', eps_nl_values, ...
                          'eps_curv_values', eps_curv_values, 'fixed_steps', fixed_steps, 'num_polish', num_polish);
results.lowfi_grad_norm = lowfi_grad_norm;
results.nl_runs = nl_runs;  results.curv_runs = curv_runs;  results.fixed_runs = fixed_runs;
results_file = fullfile(ci_dir, 'ci-mean-eps-sweep-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Mean_Eps_Sweep;

function out = Run_Mean_Eps(sen_op, hessian_analysis, z_opt, r, eps_nl, eps_curv, bt_options, num_newton)
    % Mean continuation with the a priori heuristic (+ backtracking) followed by
    % exactly num_newton Newton steps (no stopping tolerance).
    bt_options.epsilon_curv = eps_curv;
    out = Run_Mean(sen_op, hessian_analysis, z_opt, r, [], eps_nl, bt_options, num_newton, 0);
end
