clear;
close all;
clc;

load Truth_Results.mat;
load Std_OED_Results.mat;
load Seq_OED_Results.mat;

N = size(seq_oed_mean_z, 2);
std_theta_error = zeros(N, 1);
std_theta_error_z_hifi = zeros(N, 1);
std_theta_error_z_data = zeros(N, 1);
std_z_error = zeros(N, 1);
seq_theta_error = zeros(N, 1);
seq_theta_error_z_hifi = zeros(N, 1);
seq_theta_error_z_data = zeros(N, 1);
seq_z_error = zeros(N, 1);

% Generate samples from a neighborhood of z_hifi
% Normalize with magnitude 2*(z_hifi-z_lofi)
p = 100;
z_samples = z_prior_interface.Sample_with_Covariance_W_z_Inverse(p);
normalize = sqrt(diag(z_samples' * z_prior_interface.Apply_M_z(z_samples)));
r = sqrt((z_hifi - z_lofi)' * z_prior_interface.Apply_M_z(z_hifi - z_lofi));
for k = 1:p
    % z_samples(:,k) = z_lofi + 2*r*z_samples(:,k)/normalize(k);
    z_samples(:, k) = z_hifi + 2 * r * z_samples(:, k) / normalize(k);
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

delta_z_hifi_normalization = sqrt(true_delta_z_hifi' * u_prior_interface.Apply_M_u(true_delta_z_hifi));
delta_normalization = sqrt(mean(diag(true_delta' * u_prior_interface.Apply_M_u(true_delta))));
z_normalization = sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

for i = 1:N

    std_delta = zeros(m, p);
    for k = 1:p
        std_delta(:, k) = delta_eval(z_samples(:, k), std_oed_mean_theta(:, i));
    end
    diff = std_delta - true_delta;
    std_theta_error(i) = sqrt(mean(diag(diff' * u_prior_interface.Apply_M_u(diff)))) / delta_normalization;

    diff = zeros(m, i);
    for k = 1:i
        diff(:, k) = delta_eval(std_oed_Z{i}(:, k), std_oed_mean_theta(:, i)) - delta_eval(std_oed_Z{i}(:, k), best_theta);
    end
    std_theta_error_z_data(i) = sqrt(mean(diag(diff' * u_prior_interface.Apply_M_u(diff)))) / delta_normalization;

    diff = delta_eval(z_hifi, best_theta) - delta_eval(z_hifi, std_oed_mean_theta(:, i));
    std_theta_error_z_hifi(i) = sqrt(diff' * u_prior_interface.Apply_M_u(diff)) / delta_z_hifi_normalization;

    diff = std_oed_mean_z(:, i) - z_hifi;
    std_z_error(i) = sqrt(diff' * z_prior_interface.Apply_M_z(diff)) / z_normalization;

    seq_delta = zeros(m, p);
    for k = 1:p
        seq_delta(:, k) = delta_eval(z_samples(:, k), seq_oed_mean_theta(:, i));
    end
    diff = seq_delta - true_delta;
    seq_theta_error(i) = sqrt(mean(diag(diff' * u_prior_interface.Apply_M_u(diff)))) / delta_normalization;

    diff = zeros(m, i);
    for k = 1:i
        diff(:, k) = delta_eval(seq_oed_Z{i}(:, k), seq_oed_mean_theta(:, i)) - delta_eval(seq_oed_Z{i}(:, k), best_theta);
    end
    seq_theta_error_z_data(i) = sqrt(mean(diag(diff' * u_prior_interface.Apply_M_u(diff)))) / delta_normalization;

    diff = delta_eval(z_hifi, best_theta) - delta_eval(z_hifi, seq_oed_mean_theta(:, i));
    seq_theta_error_z_hifi(i) = sqrt(diff' * u_prior_interface.Apply_M_u(diff)) / delta_z_hifi_normalization;

    diff = seq_oed_mean_z(:, i) - z_hifi;
    seq_z_error(i) = sqrt(diff' * z_prior_interface.Apply_M_z(diff)) / z_normalization;

end

diff = best_z - z_hifi;
best_z_error = sqrt(diff' * z_prior_interface.Apply_M_z(diff)) / z_normalization;

delta_jac = @(theta)  reshape(theta((m + 1):end), m, m)' * Mz;

u_lofi = opt_prob_interface.State_Solve(z_lofi);
Im = eye(m);
B1 = @(x) opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(opt_prob_interface.Apply_Misfit_Hessian([Im kron(Im, z_lofi' * Mz)] * x, u_lofi, z_lofi), z_lofi);
B2 = @(x) [zeros(m, m) kron(opt_prob_interface.Misfit_Gradient(u_lofi, z_lofi)', Mz)] * x;
B = @(x) B1(x) + B2(x);
PHinvB = @(x) md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(B(x));

% x = linspace(0,1,m)';
% for k = 1:6
%     figure(1)
%     plot(x,std_oed_Z{k})
%     title('Std OED')
%     figure(2)
%     plot(x,seq_oed_Z{k})
%     title('Seq OED')
%     pause()
% end
