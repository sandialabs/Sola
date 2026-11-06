%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shared setup for the continuation-improvements experiments (Driver_CI_*).
% Builds the posterior sampling objects, one shared sensitivity operator, the
% default backtracking options, the fixed-step posterior-mean solve that
% mean-to-sample starts from, and the posterior-mean reference solve.
%
% Results are written next to these scripts (ci_dir); the OED data is read
% from article_dir.
%
% The finished setup is cached in ci-setup.mat (next to these scripts) and
% reloaded by later calls, even from a fresh MATLAB session, as long as the
% settings below are unchanged.  Delete the file to force a new setup.
%
% All solves share one MD_Continuation_Sensitivity_Operators instance with
% materialized breve samples, so every method sees identical sample
% realizations.

% This folder is Article_Drivers/Continuation_Improvements; make
% OED_Setup (one level up) and the helpers here available from any cwd.
ci_dir = fileparts(mfilename('fullpath'));
addpath(ci_dir, fileparts(ci_dir));

% ------------------------
% Shared settings
% ------------------------
num_data_points = 3;
if exist('ci_num_data_points', 'var')  % optional per-driver override (set before calling CI_Setup)
    num_data_points = ci_num_data_points;
end
cache_tag = '';                        % non-default data counts get their own cache files
if num_data_points ~= 3
    cache_tag = sprintf('-%dpt', num_data_points);
end
num_samples = 20;
fixed_mean_steps = 7;                  % mean solve used to start mean-to-sample
mean_reference_steps = 10;             % fixed-step posterior-mean reference ...
mean_reference_newton_steps = 2;       % ... plus this many Newton polishing steps at t = 1
reference_steps = 5;                   % mean-to-sample steps for the sample reference ...
reference_newton_steps = 2;            % ... plus this many Newton polishing steps at t = 1
epsilon_curv = 100;                      % curvature tolerance used by the heuristics
sample_max_steps = 15;                 % backtracking-integrator step cap for samples

setup_file = fullfile(ci_dir, ['ci-setup', cache_tag, '.mat']);
setup_key = struct('num_data_points', num_data_points, 'num_samples', num_samples, ...
                   'fixed_mean_steps', fixed_mean_steps, 'mean_reference_steps', mean_reference_steps, ...
                   'mean_reference_newton_steps', mean_reference_newton_steps, 'reference_newton_steps', reference_newton_steps, ...
                   'epsilon_curv', epsilon_curv, 'sample_max_steps', sample_max_steps);

if isfile(setup_file) && isequal(load(setup_file, 'setup_key').setup_key, setup_key)
    % Classes must be on the path before the saved objects are loaded.
    example_dir = fileparts(fileparts(ci_dir));   % Subsurface_Flow
    addpath(genpath(example_dir));
    addpath(genpath(fullfile(example_dir, '..', '..', '..', 'src')));
    load(setup_file);
    set(0, 'DefaultAxesFontSize', 20);
    set(0, 'DefaultLineLineWidth', 3);
    set(0, 'DefaultLineMarkerSize', 20);
    fprintf('Continuation improvements setup loaded from %s.\n', setup_file);
    return;
end
vars_before = who;
OED_Setup;
oed_results = load(fullfile(article_dir, 'oed-results.mat'));

data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

fprintf('\n============================================================\n');
fprintf('Continuation improvements setup (%d data points, %d samples)\n', num_data_points, num_samples);
fprintf('============================================================\n');

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

if isempty(md_hessian_analysis.evals)
    r = length(data_interface.z_opt);
else
    r = length(md_hessian_analysis.evals);
end
z_opt = data_interface.z_opt;

% Gradient norm of the standard low-fidelity reduced objective at the
% low-fidelity solution (reference level for the endpoint gradient plots).
[~, lowfi_grad] = opt_lofi.Jhat(z_opt);
lowfi_grad_norm = norm(lowfi_grad);
fprintf('Low-fidelity optimizer: ||grad Jhat_lofi|| = %.3e\n', lowfi_grad_norm);

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

fprintf('\nMean reference (%d fixed steps + %d Newton steps)...\n', mean_reference_steps, mean_reference_newton_steps);
mean_ref = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, mean_reference_steps, [], mean_bt_options, mean_reference_newton_steps);

mean_start = Run_Mean(sen_op, md_hessian_analysis, z_opt, r, fixed_mean_steps, [], mean_bt_options);
beta_bar = mean_start.beta;
fprintf('Mean solve (%d steps): %.2fs, ||grad|| = %.3e\n', fixed_mean_steps, mean_start.time, mean_start.grad_norm);

vars_to_save = setdiff(who, [vars_before; {'vars_before'}]);
save(setup_file, vars_to_save{:}, 'setup_key', '-v7.3');
fprintf('Saved setup to %s.\n', setup_file);
