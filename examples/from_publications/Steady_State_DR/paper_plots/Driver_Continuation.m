% Setup OED
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

% NOTE: THIS USES RANDOM DATA!

% % GENERATE DATA
generate_data = false;
test_approximation = false;

if generate_data
    I = eye(n);
    h = 1.e-4;
    hI = h * I;
    z_pert = z_hifi + hI;
    disp("Started...");
    tic;
    D_pert = Evaluate_Discrepancy(con_hifi, con_lofi, z_pert); % Takes ~25s
    toc;
    disp("Done...");
    D_lin = (1 / h) * (D_pert - d_hifi);
    d_hifi = Evaluate_Discrepancy(con_hifi, con_lofi, z_hifi);
    D_aff = d_hifi - D_lin * z_hifi;
    save("discrep_evals.mat", "D_lin", "D_aff", "D_pert", "z_pert", "d_hifi", "z_hifi");
end

load discrep_evals;

local_discrep_approx = @(z) D_lin * (z - z_hifi) + d_hifi;

if test_approximation
    Im = eye(m);
    Mz = z_prior_interface.M;
    D_lin_Mzinv = D_lin / Mz;
    z_test = z_hifi + norm(z_hifi - z_lofi) / sqrt(n) * randn(n, 1);
    d_test = Evaluate_Discrepancy(con_hifi, con_lofi, z_test);
    d_approx = local_discrep_approx(z_test);
    disp("Distance from z_hifi:      " + norm(z_test - z_hifi) / norm(z_test));
    disp("Distance from d_true:      " + norm(d_test - d_approx) / norm(d_test));
    best_theta = [D_aff; reshape(D_lin_Mzinv', m * n, 1)];
    eval_discrep_theta = @(z, theta) [Im kron(Im, z' * Mz)] * theta;
    discrep_est = eval_discrep_theta(z_test, best_theta);
    disp(norm(d_test - discrep_est) / norm(d_test));
end

% Perform Calibration
% Z = z_pert(:, 1:20:200);
% D = D_pert(:, 1:20:200);
% P = md_hessian_analysis.evecs * z_prior_interface.Apply_W_z(md_hessian_analysis.evecs)' / (z_prior_interface.Apply_W_z(z_lofi)' * z_lofi);
% Z = [z_lofi z_pert(:, 1:5:200)];
% Z = z_lofi + P*(z_hifi + 0.1*randn(n, 10) - z_lofi) ;
Z = [z_lofi z_hifi + 0.1 * norm(z_hifi) / sqrt(n) * randn(n, 9)];
D = local_discrep_approx(Z);
disp("Number of data points: " + size(Z, 2));
data_interface.Set_Z_and_D(Z, D);
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-3;
md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

max_cont_steps = 5;
betas = cell(max_cont_steps, 1);
zs = cell(max_cont_steps, 1);
Jhat_post = zeros(max_cont_steps, 1);
Jhat_grad_norm = zeros(max_cont_steps, 1);
MV_norm = @(beta) M_z_norm(md_hessian_analysis.evecs * beta);

% Obtain Optimal Solution Update via Continuation
for num_continuation_steps = 1:max_cont_steps
    disp(num_continuation_steps);
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    [~, z_cont, beta_cont] = md_cont_update.Posterior_Update_Mean();
    betas{num_continuation_steps} = beta_cont(:, end);
    zs{num_continuation_steps} = z_cont(:, end);
    [Jhat_post(num_continuation_steps), grad_beta] = md_cont_update.Jhat_Posterior_beta(beta_cont(:, end));
    Jhat_grad_norm(num_continuation_steps) = MV_norm(grad_beta);
end

% Try to get the true solution

options = optimoptions(@fminunc, ...
                       'Display', 'iter-detailed', ...
                       'SpecifyObjectiveGradient', true);

beta0 = beta_cont(:, end);
[beta_opt, Jhat_post_opt, ~, ~, Jhat_grad_opt, ~] = fminunc(@(beta) md_cont_update.Jhat_Posterior_beta(beta), beta0, options);
z_opt = z_lofi + md_hessian_analysis.evecs * beta_opt;
Jhat_grad_norm_opt = MV_norm(Jhat_grad_opt);

for i = 1:max_cont_steps
    fprintf("\nContinuation Steps = " + i + " \n");
    rel_z_dist(i) = M_z_norm(zs{i} - z_opt) / M_z_norm(z_opt);
    Jhat_err(i) = (Jhat_post(i) - Jhat_post_opt) / Jhat_post_opt;
    disp("Rel. Distance Error:  " + 100 * rel_z_dist(i) + "%");
    disp("Rel. Objective Diff:  " + 100 * Jhat_err(i) + "%");
    disp("Grad Norm:            " + Jhat_grad_norm(i));
end

rel_z_dist_lofi = M_z_norm(z_lofi - z_opt);
[Jhat_post_lofi, grad_beta_lofi] = md_cont_update.Jhat_Posterior_beta(0 * beta_cont(:, end));
Jhat_err_lofi = (Jhat_post_lofi - Jhat_post_opt) / Jhat_post_opt;
figure;
semilogy(0:max_cont_steps, [Jhat_err_lofi Jhat_err], ".-");
title("Jhat err.");
xlabel("Number of Continuation Steps");
ylabel("Relative Error in Objective");

% figure;
% semilogy(0:max_cont_steps, [rel_z_dist_lofi rel_z_dist], "*")
% title("Dist. to z-opt")
% figure;
% semilogy(0:max_cont_steps, [MV_norm(grad_beta_lofi); Jhat_grad_norm], "*")
% title("Grad Norm")
