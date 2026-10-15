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
num_data_points = 3;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

num_samples = 5;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

%% Posterior-mean time-step heuristic diagnostic

num_fixed_steps = 7;
fixed_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_fixed_steps);
[u_fixed, z_fixed, beta_fixed] = fixed_update.Posterior_Update_Mean();

mean_continuation_options = struct();
mean_continuation_options.discard_cache = true;
mean_continuation_options.materialize_breve_samples = false;
mean_continuation_options.use_parallel_sampling = false;
% mean_continuation_options.num_mean_to_sample_steps = 5;
mean_continuation_options.mean_time_step.use = true;
mean_continuation_options.mean_time_step.epsilon_nl = 0.01;

heuristic_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_fixed_steps, ...
                                          mean_continuation_options);
[u_heuristic, z_heuristic, beta_heuristic] = heuristic_update.Posterior_Update_Mean();

mean_prior_diagnostics = heuristic_update.last_mean_time_step_prior_diagnostics;
mean_time_grid = heuristic_update.last_mean_time_grid;
mean_step_diagnostics = heuristic_update.last_mean_time_step_diagnostics;
mean_dt_values = diff(mean_time_grid);
num_mean_heuristic_steps = length(mean_dt_values);

J_fixed = opt_hifi.Jhat(z_fixed);
J_heuristic = opt_hifi.Jhat(z_heuristic);
J_lofi = opt_hifi.Jhat(z_lofi);

fprintf('\nMean time-step heuristic diagnostics (subsurface flow)\n');
fprintf('----------------------------------------------------\n');
fprintf('Fixed continuation steps:       %d\n', num_fixed_steps);
fprintf('Heuristic continuation steps:   %d\n', num_mean_heuristic_steps);
fprintf('A priori dt:                    %.4e\n', mean_prior_diagnostics.dt);
fprintf('Accepted dt range:              [%.4e, %.4e]\n', min(mean_dt_values), max(mean_dt_values));
fprintf('A priori A, eta:                %.4e, %.4e\n', mean_prior_diagnostics.A, mean_prior_diagnostics.eta);
fprintf('A priori dt_nl, dt_curv:        %.4e, %.4e\n', mean_prior_diagnostics.dt_nonlinearity, mean_prior_diagnostics.dt_curvature);
fprintf('HiFi objective at LoFi opt:     %.6e\n', J_lofi);
fprintf('HiFi objective fixed mean:      %.6e\n', J_fixed);
fprintf('HiFi objective heuristic mean:  %.6e\n', J_heuristic);
fprintf('||z_heuristic - z_fixed||/||z_fixed|| = %.4e\n', ...
        norm(z_heuristic - z_fixed) / max(norm(z_fixed), eps));

if ~isempty(mean_step_diagnostics)
    fprintf('\nMean backtracking diagnostics:\n');
    fprintf('  step       t_start        dt      backtracks   ||g_pred||   ||g_corr||\n');
    for k = 1:length(mean_step_diagnostics)
        d = mean_step_diagnostics{k};
        fprintf('  %4d    %10.4e  %10.4e      %4d    %10.4e  %10.4e\n', ...
                k, d.t_start, d.dt_used, d.num_backtracks, ...
                d.grad_pred_norm, d.grad_corr_norm);
    end
end

%% Posterior-sample mean-to-sample time-step heuristic diagnostic

sample_continuation_options = struct();
sample_continuation_options.discard_cache = true;
sample_continuation_options.materialize_breve_samples = true;
sample_continuation_options.use_parallel_sampling = false;
% sample_continuation_options.num_mean_to_sample_steps = 4;
sample_continuation_options.sample_time_step.use = true;
sample_continuation_options.sample_time_step.epsilon_nl = 0.1;

md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_fixed_steps, ...
                                        sample_continuation_options);

% Reuse the fixed-step mean as the mean-to-sample starting point so this block
% isolates the sample-specific heuristic.
u_bar = u_fixed;
z_bar = z_fixed;
beta_bar = beta_fixed;
[u_samples, z_samples, beta_samples] = md_cont_update.Posterior_Update_Samples_From_Mean(u_bar, z_bar, beta_bar);

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
    assert(numel(diff(tg)) <= sample_max_steps, sprintf('Sample %d exceeded max sample steps.', sample_idx));

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
        assert(sd{k}.grad_corr_norm <= max(sample_heuristic_options.gradient_abs_tol, sample_heuristic_options.gradient_reduction_factor * sd{k}.grad_pred_norm) || sd{k}.objective_accept, ...
               sprintf('Sample %d step %d did not satisfy acceptance diagnostics.', sample_idx, k));
    end
end

if max(sample_dt_prior) - min(sample_dt_prior) <= 1.e-12
    warning('SampleTimeStepDiagnostic:NoVariation', ...
            'All sample a priori dt values are numerically identical for this run.');
end

J_samples = zeros(num_samples, 1);
for sample_idx = 1:num_samples
    J_samples(sample_idx) = opt_hifi.Jhat(z_samples(:, sample_idx));
end
J_bar = opt_hifi.Jhat(z_bar);

fprintf('\nSample time-step heuristic diagnostics (subsurface flow)\n');
fprintf('------------------------------------------------------\n');
fprintf('Number of samples:             %d\n', num_samples);
fprintf('Mean-to-sample base steps:      %d\n', num_mean_to_sample_steps);
fprintf('LoFi HiFi objective:           %.6e\n', J_lofi);
fprintf('Posterior mean HiFi objective: %.6e\n', J_bar);
fprintf('\nPer-sample a priori and accepted-step diagnostics:\n');
fprintf(' sample    dt_prior        A          eta       dt_nl     dt_curv   steps  backtracks   accepted dt range\n');
for sample_idx = 1:num_samples
    fprintf(' %4d   %10.4e  %9.3e  %9.3e  %9.3e  %9.3e  %5d     %5d     [%9.3e, %9.3e]\n', ...
            sample_idx, sample_dt_prior(sample_idx), sample_A(sample_idx), sample_eta(sample_idx), ...
            sample_dt_nl(sample_idx), sample_dt_curv(sample_idx), sample_num_steps(sample_idx), ...
            sample_total_backtracks(sample_idx), sample_min_dt(sample_idx), sample_max_dt(sample_idx));
end
fprintf('\nSample HiFi objective range:    [%.6e, %.6e]\n', min(J_samples), max(J_samples));
fprintf('Sample dt_prior range:          [%.4e, %.4e]\n', min(sample_dt_prior), max(sample_dt_prior));
fprintf('Sample eta range:               [%.4e, %.4e]\n', min(sample_eta), max(sample_eta));

save('time_step_heuristic_subsurface_diagnostics.mat');
fprintf('\nSaved diagnostics to %s\n', output_file);

