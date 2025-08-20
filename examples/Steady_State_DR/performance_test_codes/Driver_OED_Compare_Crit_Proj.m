% Get the OED Setup ready
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

load("oed-results.mat");

% Perform Offline OED Computations - USES data_interface
md_oed = MD_OED_DeltaCov(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

%% Perform Random Data Point Comparison
N = 5;
num_samples = 5;
crit_rand = zeros(N, num_samples);
crit_DC_oed = zeros(N, 1);
crit_same_pt = zeros(N, 1);
betas_oed = md_hessian_analysis.evecs \ (Z_oed - z_lofi);
beta_bars = md_hessian_analysis.evecs \ (z_bars - z_lofi);
beta_bar_hifi = md_hessian_analysis.evecs \ (z_hifi - z_lofi);
rng(0);

wb = waitbar(0, 'Starting Random');
for p = 2:N
    % fprintf('\nStep %d:\n-------------\n', p);
    % Generate random data points
    if p == 2
        covar_coeff = W_z_norm(z_bars(:, p - 1) - z_lofi)^2 / n;
    else
        covar_coeff = W_z_norm(z_bars(:, p - 1) - z_bars(:, p - 2))^2 / n;
    end
    md_oed.Set_Covariance_Coefficient(covar_coeff);

    % Also get optimizer's criterion:
    crit_fn = @(beta_input) md_oed.Evaluate_OED_Objective(reshape([betas_oed(:, 2:p - 1) beta_input], [], 1), alpha_d, 0, beta_bars(:, p - 1));
    crit_DC_oed(p) = crit_fn(betas_oed(:, p));
    crit_same_pt(p) = crit_fn(zeros(num_evals, 1));

    % Obtains random points in a neighborhood of zbar from a projected subspace, where distance is dictated by covar_coeff
    perturb = sqrt(covar_coeff * (n / num_evals)) / W_z_norm(z_lofi) * md_hessian_analysis.evecs * randn(num_evals, num_samples);
    z_ps = z_bars(:, p - 1) + sqrt(covar_coeff * (n / num_evals)) / W_z_norm(z_lofi) * md_hessian_analysis.evecs * randn(num_evals, num_samples);

    for i = 1:num_samples
        % disp(i)
        waitbar((num_samples * (p - 2) + i) / (num_samples * (N - 1)), wb, sprintf('Random Progress: %d %%', floor((num_samples * (p - 2) + i) / (num_samples * (N - 1)) * 100)));
        % Obtain random data point
        z_p = z_ps(:, i);
        beta_p = md_hessian_analysis.evecs \ (z_p - z_lofi);
        % Z_rand = [Z_oed(:, 1:p - 1) z_p];

        % Obtain Optimal Solution Update via Continuation
        crit_rand(p, i) = crit_fn(beta_p);
    end
    % waitbar((p-1) / (N-1), wb, sprintf('Random Progress: %d %%', floor((p-1) / (N-1) * 100)));

end
close(wb);

figure;
plot(2:N, (crit_rand(2:end, :)' - crit_DC_oed(2:end)') ./ (crit_same_pt(2:end)' - crit_DC_oed(2:end)'), ".", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "HandleVisibility", "off");
% hold on;
% xlim([2 N]);
% % plot(2:N, crit_DC_oed, ".-", "Color", "#BAB86C", "DisplayName", "OED");
% xlabel("Evaluations ($N$)", "Interpreter", "latex");
% ylabel("Criterion", "Interpreter", "latex");
% legend("location", "east", "Interpreter", "latex");
% title("Optimization Objective over Evals");
% saveas(gcf, "oed-rand-proj.eps", "epsc");
