%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Diagnostics for posterior-mean and posterior-sample continuation with the
% a priori time-step heuristic and pseudo-time backtracking safeguard.
%
% The mean diagnostic compares fixed-step and heuristic mean continuation.
% The sample diagnostic then checks sample-specific mean-to-sample continuation
% using each sample's realized centered discrepancy.

OED_Setup;

oed_results = load('oed-results.mat');
num_data_points = 1;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

num_samples = 5;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis);

% Newton polishing steps at t = 1 after the continuation steps (fixed and heuristic).
num_newton_steps = 2;

%% Posterior-mean time-step heuristic diagnostic

fixed_mean_options = struct();
fixed_mean_options.num_steps = 7;
fixed_mean_options.num_newton_steps = num_newton_steps;
[u_fixed, z_fixed, beta_fixed] = md_cont_update.Posterior_Update_Mean(fixed_mean_options);

heuristic_mean_options = struct();
heuristic_mean_options.num_newton_steps = num_newton_steps;
heuristic_mean_options.time_step.use = true;
heuristic_mean_options.time_step.epsilon_nl = 0.01;
[u_heuristic, z_heuristic, beta_heuristic] = md_cont_update.Posterior_Update_Mean(heuristic_mean_options);

mean_prior_diagnostics = md_cont_update.last_mean_time_step_prior_diagnostics;
mean_time_grid = md_cont_update.last_mean_time_grid;
mean_dt_values = diff(mean_time_grid);
num_mean_heuristic_steps = length(mean_dt_values);

J_fixed = opt_hifi.Jhat(z_fixed);
J_heuristic = opt_hifi.Jhat(z_heuristic);
J_lofi = opt_hifi.Jhat(z_lofi);
[J_fixed_param, grad_fixed_param] = md_cont_update.Parameterized_RS_Objective_beta(beta_fixed, 0);
[J_heuristic_param, grad_heuristic_param] = md_cont_update.Parameterized_RS_Objective_beta(beta_heuristic, 0);

fprintf('\nMean time-step heuristic diagnostics (subsurface flow)\n');
fprintf('----------------------------------------------------\n');
fprintf('Fixed continuation steps:       %d\n', fixed_mean_options.num_steps);
fprintf('Heuristic continuation steps:   %d\n', num_mean_heuristic_steps);
fprintf('A priori dt:                    %.4e\n', mean_prior_diagnostics.dt);
fprintf('Accepted dt range:              [%.4e, %.4e]\n', min(mean_dt_values), max(mean_dt_values));
fprintf('A priori A, eta:                %.4e, %.4e\n', mean_prior_diagnostics.A, mean_prior_diagnostics.eta);
fprintf('A priori dt_nl, dt_curv:        %.4e, %.4e\n', mean_prior_diagnostics.dt_nonlinearity, mean_prior_diagnostics.dt_curvature);
fprintf('HiFi objective at LoFi opt:     %.6e\n', J_lofi);
fprintf('HiFi objective fixed mean:      %.6e\n', J_fixed);
fprintf('HiFi objective heuristic mean:  %.6e\n', J_heuristic);
fprintf('Param objective fixed mean:     %.6e  (||g_beta|| = %.4e)\n', J_fixed_param, norm(grad_fixed_param));
fprintf('Param objective heuristic mean: %.6e  (||g_beta|| = %.4e)\n', J_heuristic_param, norm(grad_heuristic_param));
fprintf('||z_heuristic - z_fixed||/||z_fixed|| = %.4e\n', ...
        norm(z_heuristic - z_fixed) / max(norm(z_fixed), eps));

%% Posterior-sample mean-to-sample time-step heuristic diagnostic

sample_options = struct();
% sample_options.num_mean_to_sample_steps = 4;
sample_options.materialize_breve_samples = true;
sample_options.num_newton_steps = num_newton_steps;
sample_options.time_step.use = true;
sample_options.time_step.epsilon_nl = 0.01;
sample_options.time_step.options.epsilon_curv = 1;

% Reuse the fixed-step mean as the mean-to-sample starting point so this block
% isolates the sample-specific heuristic.
u_bar = u_fixed;
z_bar = z_fixed;
beta_bar = beta_fixed;
[u_samples, z_samples, beta_samples] = md_cont_update.Posterior_Update_Samples_From_Mean(u_bar, z_bar, beta_bar, sample_options);
sample_options = md_cont_update.options.samples;

sample_prior_diags = md_cont_update.last_sample_time_step_prior_diagnostics;
sample_time_grids = md_cont_update.last_sample_time_grids;
sample_step_diags = md_cont_update.last_sample_time_step_diagnostics;

sample_dt_prior = zeros(num_samples, 1);
sample_A = zeros(num_samples, 1);
sample_eta = zeros(num_samples, 1);
sample_dt_nl = zeros(num_samples, 1);
sample_dt_curv = zeros(num_samples, 1);
sample_num_steps = zeros(num_samples, 1);
sample_total_backtracks = zeros(num_samples, 1);
sample_min_dt = zeros(num_samples, 1);
sample_max_dt = zeros(num_samples, 1);

for sample_idx = 1:num_samples
    d = sample_prior_diags{sample_idx};
    tg = sample_time_grids{sample_idx};
    sd = sample_step_diags{sample_idx};

    assert(isfield(d, 'dt'), sprintf('Missing prior dt for sample %d.', sample_idx));
    assert(d.dt > 0 && d.dt <= 1, sprintf('Invalid prior dt for sample %d.', sample_idx));
    assert(abs(tg(1)) < 1.e-14, sprintf('Sample %d time grid does not start at zero.', sample_idx));
    assert(abs(tg(end) - 1.0) < 1.e-12, sprintf('Sample %d time grid does not end at one.', sample_idx));
    assert(all(diff(tg) > 0), sprintf('Sample %d time grid is not strictly increasing.', sample_idx));
    assert(numel(diff(tg)) <= sample_options.time_step.options.max_steps, sprintf('Sample %d exceeded max sample steps.', sample_idx));

    sample_dt_prior(sample_idx) = d.dt;
    sample_A(sample_idx) = d.A;
    sample_eta(sample_idx) = d.eta;
    sample_dt_nl(sample_idx) = d.dt_nonlinearity;
    sample_dt_curv(sample_idx) = d.dt_curvature;
    sample_num_steps(sample_idx) = numel(diff(tg));
    sample_min_dt(sample_idx) = min(diff(tg));
    sample_max_dt(sample_idx) = max(diff(tg));

    for k = 1:numel(sd)
        sample_total_backtracks(sample_idx) = sample_total_backtracks(sample_idx) + sd{k}.num_backtracks;
        assert(sd{k}.accepted, sprintf('Sample %d step %d was not accepted.', sample_idx, k));
        assert(sd{k}.grad_corr_norm <= max(sample_options.time_step.options.gradient_abs_tol, sample_options.time_step.options.gradient_reduction_factor * sd{k}.grad_pred_norm) || sd{k}.objective_accept, ...
               sprintf('Sample %d step %d did not satisfy acceptance diagnostics.', sample_idx, k));
    end
end

if max(sample_dt_prior) - min(sample_dt_prior) <= 1.e-12
    warning('SampleTimeStepDiagnostic:NoVariation', ...
            'All sample a priori dt values are numerically identical for this run.');
end

sample_grad_norm = md_cont_update.last_sample_grad_norms;   % ||g_beta|| at t = 1, after Newton polishing
J_bar = opt_hifi.Jhat(z_bar);

fprintf('\nSample time-step heuristic diagnostics (subsurface flow)\n');
fprintf('------------------------------------------------------\n');
fprintf('Number of samples:             %d\n', num_samples);
fprintf('Mean-to-sample base steps:      %d\n', sample_options.num_mean_to_sample_steps);
fprintf('LoFi HiFi objective:           %.6e\n', J_lofi);
fprintf('Posterior mean HiFi objective: %.6e\n', J_bar);
fprintf('Posterior mean param objective: %.6e  (||g_beta|| = %.4e)\n', J_fixed_param, norm(grad_fixed_param));
fprintf('\nPer-sample a priori and accepted-step diagnostics:\n');
fprintf(' sample    dt_prior        A          eta       dt_nl     dt_curv   steps  backtracks   accepted dt range\n');
for sample_idx = 1:num_samples
    fprintf(' %4d   %10.4e  %9.3e  %9.3e  %9.3e  %9.3e  %5d     %5d     [%9.3e, %9.3e]\n', ...
            sample_idx, sample_dt_prior(sample_idx), sample_A(sample_idx), sample_eta(sample_idx), ...
            sample_dt_nl(sample_idx), sample_dt_curv(sample_idx), sample_num_steps(sample_idx), ...
            sample_total_backtracks(sample_idx), sample_min_dt(sample_idx), sample_max_dt(sample_idx));
end
fprintf('\nPer-sample endpoint (t = 1, after %d Newton polishing steps) gradient norm:\n', num_newton_steps);
fprintf(' sample    ||g_beta||\n');
for sample_idx = 1:num_samples
    fprintf(' %4d     %10.4e\n', sample_idx, sample_grad_norm(sample_idx));
end
fprintf('\nSample ||g_beta|| range:        [%.4e, %.4e]\n', min(sample_grad_norm), max(sample_grad_norm));
fprintf('Sample dt_prior range:          [%.4e, %.4e]\n', min(sample_dt_prior), max(sample_dt_prior));
fprintf('Sample eta range:               [%.4e, %.4e]\n', min(sample_eta), max(sample_eta));

output_file = 'time_step_heuristic_subsurface_diagnostics.mat';
save(output_file);
fprintf('\nSaved diagnostics to %s\n', output_file);

