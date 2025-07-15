% Setup OED
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

% % GENERATE DATA
% I = eye(n);
% h = 1.e-4;
% hI = h*I;
% z_pert = z_hifi + hI;
% disp("Started...")
% tic
% D_pert = Evaluate_Discrepancy(con_hifi, con_lofi, z_pert); % Takes ~25s
% toc
% disp("Done...")
% D_lin = (1/h) * (D_pert - d_hifi)
% d_hifi = Evaluate_Discrepancy(con_hifi, con_lofi, z_hifi);
% D_aff = d_hifi - D_lin*z_hifi;
% save("discrep_evals.mat", "D_lin", "D_aff", "D_pert", "z_pert", "d_hifi", "z_hifi")

load discrep_evals;

% % Check if discrepancy approximation locally is good
% Im = eye(m);
% Mz = z_prior_interface.M;
% D_lin_Mzinv = D_lin / Mz;
% z_test = z_hifi + norm(z_hifi-z_lofi)/sqrt(n)*randn(n, 1);
% d_test = Evaluate_Discrepancy(con_hifi, con_lofi, z_test);
% local_discrep_approx = @(z) D_lin * (z - z_hifi) + d_hifi;
% d_approx = local_discrep_approx(z_test);
% disp("Distance from z_hifi:      " + norm(z_test - z_hifi)/norm(z_test))
% disp("Distance from d_true:      " + norm(d_test - d_approx)/norm(d_test))
% best_theta = [D_aff; reshape(D_lin_Mzinv', m * n, 1)];
% eval_discrep_theta = @(z, theta) [Im kron(Im, z' * Mz)] * theta;
% discrep_est = eval_discrep_theta(z_test, best_theta);
% disp(norm(d_test-discrep_est)/norm(d_test))

% Perform Calibration
Z = z_pert(:, 1:10:200);
D = D_pert(:, 1:10:200);
disp("Number of data points: " + size(Z, 2));
data_interface.Set_Z_and_D(Z, D);
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-6;
md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

% Obtain Optimal Solution Update via Continuation
num_continuation_steps = 1;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[~, z_cont, ~] = md_cont_update.Posterior_Update_Mean_beta();
z_bar = z_cont(:, end);

fprintf("\nContinuation Steps = " + num_continuation_steps + " \n");
disp("Rel. Distance Error:  " + 100 * norm(z_bar - z_hifi) / norm(z_hifi) + "%");
disp("Rel. Objective Diff:  " + 100 * (opt_hifi.Jhat(z_bar) - Jhat_hifi) / Jhat_hifi + "%");

% Obtain Optimal Solution Update via Continuation
num_continuation_steps = 2;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[~, z_cont, ~] = md_cont_update.Posterior_Update_Mean_beta();
z_bar = z_cont(:, end);

fprintf("\nContinuation Steps = " + num_continuation_steps + " \n");
disp("Rel. Distance Error:  " + 100 * norm(z_bar - z_hifi) / norm(z_hifi) + "%");
disp("Rel. Objective Diff:  " + 100 * (opt_hifi.Jhat(z_bar) - Jhat_hifi) / Jhat_hifi + "%");
