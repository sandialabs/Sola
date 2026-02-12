% Get the OED Setup ready
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

load("oed-results-con2-new.mat");

% Perform Offline OED Computations - USES data_interface
md_oed = MD_OED_DeltaCov(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();
V = md_hessian_analysis.evecs;
alpha_k_denom = trace(z_prior_interface.Apply_W_z_Inverse(z_prior_interface.M));
alpha_k_denom_proj = trace(V' * z_prior_interface.M * V) / W_z_norm(z_lofi)^2;

% Perform Random Data Point Comparison
N = 5;
num_samples = 100;
crit_rand = zeros(num_samples, 1);
crit_rand_hifi = zeros(num_samples, 1);
betas_oed = V \ (Z_oed - z_lofi);
beta_bars = V \ (z_bars - z_lofi);
beta_bar_hifi = V \ (z_hifi - z_lofi);
rng(0);

p = 3;
alpha_k_num = M_z_norm(z_bars(:, p - 1) - z_bars(:, p - 2))^2;
alpha_k = alpha_k_num / alpha_k_denom;
alpha_k_proj = alpha_k_num / alpha_k_denom_proj;
md_oed.Set_Covariance_Coefficient(alpha_k);

% Also get optimizer's criterion:
crit_fn = @(beta_input) md_oed.Evaluate_OED_Objective(reshape([betas_oed(:, 2:p - 1) beta_input], [], 1), alpha_d, 0, beta_bars(:, p - 1));
crit_hifi = @(beta_input) md_oed.Evaluate_OED_Objective(reshape([betas_oed(:, 2:p - 1) beta_input], [], 1), alpha_d, 0, beta_bar_hifi);

crit_DC_oed = crit_fn(betas_oed(:, p));
crit_DC_oed_hifi = crit_hifi(betas_oed(:, p));
dist_DC_oed_hifi = W_z_norm(z_lofi + V * betas_oed(:, p) - z_hifi) / W_z_norm(z_hifi);
dist_DC_oed_bar = W_z_norm(V * (betas_oed(:, p) - beta_bars(:, p - 1))) / W_z_norm(V * beta_bars(:, p - 1));
crit_same_pt = crit_fn(zeros(num_evals, 1));
crit_same_pt_hifi = crit_hifi(zeros(num_evals, 1));

% M_z_norm(z_lofi + V * betas_oed(:, p) - z_bars(:, p - 1)))^2
% Obtains random points in a neighborhood of zbar from a projected subspace, where distance is dictated by covar_coeff
% z_ps = z_bars(:, p - 1) + 1*sqrt(covar_coeff * (n / num_evals)) / W_z_norm(z_lofi) * V * randn(num_evals, num_samples);
% (V' *  z_prior_interface.Apply_W_z(V))/W_z_norm(z_lofi)^2 = I
z_ps = z_bars(:, p - 1) + sqrt(alpha_k_proj) / W_z_norm(z_lofi) * V * randn(num_evals, 2 * num_samples);
% z_ps = z_ps(:, diag(M_z_norm(z_ps - z_bars(:, p - 1))).^2 < 100*alpha_k_num); % NOTE: no rejection sampling
% disp("Rejection Sampling, accepted pct: " + size(z_ps, 2) / (2 * num_samples));
z_ps = z_ps(:, 1:num_samples);

for i = 1:num_samples
    z_p = z_ps(:, i);
    beta_p = V \ (z_p - z_lofi);
    crit_rand(i) = crit_fn(beta_p);
    crit_rand_hifi(i) = crit_hifi(beta_p);
    dist_rand_hifi(i) = W_z_norm(z_lofi + V * beta_p - z_hifi) / W_z_norm(z_hifi);
    dist_rand_bar(i) = W_z_norm(V * (beta_p - beta_bars(:, p - 1))) / W_z_norm(V * beta_bars(:, p - 1));
end

% md_oed.Set_Covariance_Coefficient(0);
% crit_hifi = @(beta_input) md_oed.Evaluate_OED_Objective(reshape([betas_oed(:, 2:p - 1) beta_input], [], 1), alpha_d, 0, beta_bar_hifi);
% crit_DC_oed_hifi = crit_hifi(betas_oed(:, p));
% crit_same_pt_hifi = crit_hifi(zeros(num_evals, 1));

% for i = 1:num_samples
%     z_p = z_ps(:, i);
%     beta_p = V \ (z_p - z_lofi);
%     crit_rand_hifi(i) = crit_hifi(beta_p);
% end

% figure;
% plot(crit_rand - crit_same_pt,  dist_rand_hifi, ".", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "HandleVisibility", "off");
% hold on;
% plot(crit_DC_oed - crit_same_pt, dist_DC_oed_hifi, "r*", "DisplayName", "Optimal Data Point");
% xlabel("OED Criterion"); ylabel("Distance to z_{hifi}"); legend;

% figure;
% plot(crit_rand - crit_same_pt,  dist_rand_bar, ".", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "HandleVisibility", "off");
% hold on;
% plot(crit_DC_oed - crit_same_pt, dist_DC_oed_bar, "r*", "DisplayName", "Optimal Data Point");
% xlabel("OED Criterion"); ylabel("Distance to z_{bar}"); legend;

% figure;
% plot((crit_rand - crit_same_pt) / median(crit_rand - crit_same_pt),  (crit_rand_hifi - crit_same_pt_hifi) / median(crit_rand_hifi - crit_same_pt_hifi), ".", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "HandleVisibility", "off");
% hold on;
% plot((crit_DC_oed - crit_same_pt) / median(crit_rand - crit_same_pt), (crit_DC_oed_hifi - crit_same_pt_hifi) / median(crit_rand_hifi - crit_same_pt_hifi), "r*", "DisplayName", "Optimal Data Point");
% xlabel("Uncertainty reduction near $\bar{z}_k$", "Interpreter", "latex");
% ylabel("Uncertainty reduction near $z^{\star}$", "Interpreter", "latex");
% legend('Location', 'northwest');

figure;
plot((crit_rand - crit_same_pt) / (crit_DC_oed - crit_same_pt),  (crit_rand_hifi - crit_same_pt_hifi) / (crit_DC_oed_hifi - crit_same_pt_hifi), ".", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "HandleVisibility", "off");
hold on;
plot((crit_DC_oed - crit_same_pt) / (crit_DC_oed - crit_same_pt), (crit_DC_oed_hifi - crit_same_pt_hifi) / (crit_DC_oed_hifi - crit_same_pt_hifi), "r*", "DisplayName", "Optimal Data Point");
xlabel("Uncertainty reduction near $\bar{z}_k$", "Interpreter", "latex");
ylabel("Uncertainty reduction near $z^{\star}$", "Interpreter", "latex");
legend('Location', 'northwest');
ylim([0 inf]);
% saveas(gcf, 'ScatterPlot2', 'epsc');

% figure;
% scatter3(crit_rand - crit_same_pt,  crit_rand_hifi - crit_same_pt_hifi,  dist_rand_hifi, ".", "Color", [0.7 0.8 0.9], "HandleVisibility", "off");
% hold on;
% scatter3(crit_DC_oed - crit_same_pt, crit_DC_oed_hifi - crit_same_pt_hifi, dist_DC_oed_hifi, "r*", "DisplayName", "Optimal Data Point");
% xlabel("OED Criterion"); ylabel("Expected Uncertainty near z_{hifi}"); zlabel("Distance to z_{hifi}"); legend;

% figure;
% scatter3(crit_rand - crit_same_pt,  crit_rand_hifi - crit_same_pt_hifi,  dist_rand_bar, ".", "Color", [0.7 0.8 0.9], "HandleVisibility", "off");
% hold on;
% scatter3(crit_DC_oed - crit_same_pt, crit_DC_oed_hifi - crit_same_pt_hifi, dist_DC_oed_bar, "r*", "DisplayName", "Optimal Data Point");
% xlabel("OED Criterion"); ylabel("Expected Uncertainty near z_{hifi}"); zlabel("Distance to z_{bar}"); %legend;
