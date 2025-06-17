clear;
close all;
clc;
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

load Truth_Results.mat;
load Std_OED_Results.mat;
load Seq_OED_Results.mat;

N = length(Jhat_std_oed);

std_theta_error = zeros(N, 1);
std_theta_norm_error = zeros(N, 1);
std_theta_error_z_hifi = zeros(N, 1);
std_theta_error_z_data = zeros(N, 1);
std_z_error = zeros(N, 1);

seq_theta_error = zeros(N, 1);
seq_theta_norm_error = zeros(N, 1);
seq_theta_error_z_hifi = zeros(N, 1);
seq_theta_error_z_data = zeros(N, 1);
seq_z_error = zeros(N, 1);

znorm = @(z) sqrt(sum(z .* z_prior_interface.Apply_M_z(z), 1))';
unorm = @(u) sqrt(sum(u .* u_prior_interface.Apply_M_u(u), 1))';
batch_unorm = @(u) norm(unorm(u)) / sqrt(size(u, 2));

% Generate samples from a neighborhood of z_hifi
% Normalize with magnitude 2*(z_hifi-z_lofi)
p = 100;
z_samples = z_prior_interface.Sample_with_Covariance_W_z_Inverse(p);
normalize = znorm(z_samples);
r = znorm(z_hifi - z_lofi);
for k = 1:p
    z_samples(:, k) = z_hifi + 2 * r * z_samples(:, k) / normalize(k); % z_lofi
end

m = length(z_lofi);
Im = eye(m);
Mz = z_prior_interface.M;
delta_eval = @(z, theta) [Im kron(Im, z' * Mz)] * theta;
true_delta = zeros(m, p);
for k = 1:p
    true_delta(:, k) = delta_eval(z_samples(:, k), best_theta);
end
true_delta_z_hifi = delta_eval(z_hifi, best_theta);

delta_z_hifi_normalization = unorm(true_delta_z_hifi);
delta_normalization = batch_unorm(true_delta);
z_normalization = znorm(z_hifi);

for i = 1:N
    disp(i);
    % Obtain error in delta at samples (standard OED)
    std_delta = zeros(m, p);
    for k = 1:p
        std_delta(:, k) = delta_eval(z_samples(:, k), std_oed_mean_theta{i});
    end
    std_theta_error(i) = batch_unorm(std_delta - true_delta) / delta_normalization;

    % Obtain error in delta at data points (standard OED)
    diff = zeros(m, i);
    for k = 1:i
        diff(:, k) = delta_eval(std_oed_Z{i}(:, k), std_oed_mean_theta{i}) - delta_eval(std_oed_Z{i}(:, k), best_theta);
    end
    std_theta_error_z_data(i) = batch_unorm(diff) / delta_normalization;

    % Obtain error in delta at z_hifi (standard OED)
    diff = delta_eval(z_hifi, best_theta) - delta_eval(z_hifi, std_oed_mean_theta{i});
    std_theta_error_z_hifi(i) = unorm(diff) / delta_z_hifi_normalization;

    % Obtain error in z_update (standard OED)
    diff = std_oed_mean_z{i} - z_hifi;
    std_z_error(i) = znorm(diff) / z_normalization;

    % Obtain error in delta at samples (sequential OED)
    seq_delta = zeros(m, p);
    for k = 1:p
        seq_delta(:, k) = delta_eval(z_samples(:, k), seq_oed_mean_theta{i});
    end
    diff = seq_delta - true_delta;
    seq_theta_error(i) = batch_unorm(diff) / delta_normalization;

    % Obtain error in delta at data points (sequential OED)
    diff = zeros(m, i);
    for k = 1:i
        diff(:, k) = delta_eval(seq_oed_Z{i}(:, k), seq_oed_mean_theta{i}) - delta_eval(seq_oed_Z{i}(:, k), best_theta);
    end
    seq_theta_error_z_data(i) = batch_unorm(diff) / delta_normalization;

    % Obtain error in delta at z_hifi (sequential OED)
    diff = delta_eval(z_hifi, best_theta) - delta_eval(z_hifi, seq_oed_mean_theta{i});
    seq_theta_error_z_hifi(i) = unorm(diff) / delta_z_hifi_normalization;

    % Obtain error in z_update (sequential OED)
    diff = seq_oed_mean_z{i} - z_hifi;
    seq_z_error(i) = znorm(diff) / z_normalization;

end

best_z_error = znorm(best_z - z_hifi) / z_normalization;

delta_jac = @(theta)  reshape(theta((m + 1):end), m, m)' * Mz;

u_lofi = opt_prob_interface.State_Solve(z_lofi);
Im = eye(m);
B1 = @(x) opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(opt_prob_interface.Apply_Misfit_Hessian([Im kron(Im, z_lofi' * Mz)] * x, u_lofi, z_lofi), z_lofi);
B2 = @(x) [zeros(m, m) kron(opt_prob_interface.Misfit_Gradient(u_lofi, z_lofi)', Mz)] * x;
B = @(x) B1(x) + B2(x);
PHinvB = @(x) md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(B(x));

% disp(norm(std_theta_error)-1.2896)
% disp(norm(std_theta_error_z_hifi)-1.2730)
% disp(norm(std_theta_error_z_data)-0.0021)
% disp(norm(std_z_error)-0.1896)
% disp(norm(seq_theta_error)-1.2896)
% disp(norm(seq_theta_error_z_hifi)-1.2730)
% disp(norm(seq_theta_error_z_data)-0.0021)
% disp(norm(seq_z_error)-0.1896)

x = linspace(0, 1, m)';
for k = 1:6
    figure(1);
    plot(x, std_oed_Z{k});
    title('Std OED');
    figure(2);
    plot(x, seq_oed_Z{k});
    title('Seq OED');
    pause();
end
