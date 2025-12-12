% Set Default Font Axes and Line Width
OED_Setup;

% Perform Offline OED Computations
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();
md_oed.verbosity = false;

% Set Parameters for OED
N = 5;
rng(0);
beta_0 = randn(num_evals * (N - 1), 1);
reg_coeff = 1.e-6;
% [betas, Z] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);

% Generate Design (Generate_Random_Design(N), Generate_Random_Design_from_Subspace(N), Generate_Optimal_Design(...))
Z = md_oed.Generate_Random_Design_from_Subspace(N);
% D = con_hifi.State_Solve(Z) - con_lofi.State_Solve(Z);
D = Evaluate_Discrepancy(con_hifi, con_lofi, Z);
data_interface.Set_Z_and_D(Z, D);

% % Sample from Posterior (i.e., solve problem) - USES DATA INTERFACE
num_post_samples = 1;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
% [delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z);

% Methods to compute H^{-1}B at lofi solution
Im = eye(m);
Mz = z_prior_interface.M;
B1_lofi = @(x) opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(opt_prob_interface.Apply_Misfit_Hessian([Im kron(Im, z_lofi' * Mz)] * x, u_lofi, z_lofi), z_lofi);
B2_lofi = @(x) [zeros(n, n) kron(opt_prob_interface.Misfit_Gradient(u_lofi, z_lofi)', Mz)] * x;
B_lofi = @(x) B1_lofi(x) + B2_lofi(x);
PHinvB_lofi = @(x) md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(B_lofi(x));
HinvB_lofi = @(theta, z) md_hessian_analysis.Apply_RS_Hessian_Inverse_CG(B_lofi(theta), z);

% Extract posterior mean of theta (and its HDSA propagation)
theta_post = Extract_mean_theta(md_post_sampling.post_data);
z_post = z_lofi - PHinvB_lofi(theta_post);

% Implement H^{-1}B for general theta
% Subfunctions
theta_zero = zeros(size(theta_post));
discrep_eval = @(z, theta) [Im kron(Im, z' * Mz)] * theta;
discrep_z = @(theta) reshape(theta(m + 1:end)', n, m)' * Mz;
discrep_zt = @(theta_in) discrep_z(theta_in);
u_solve = @(z, theta) con_lofi.State_Solve(z) + discrep_eval(z, theta);
J_u = @(u, z) J_u_fn(obj, u, z);
J_uu = @(u_in, u, z) opt_prob_interface.Apply_Misfit_Hessian(u_in, u, z);
J_z = @(u, z) J_z_fn(obj, u, z);

% Compute action of B
B1_p1 = @(u_in, z, u, theta) opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_in, z) + discrep_z(theta)' * u_in;
% B1 = @(theta_in, z, u, theta) B1_p1(obj.J_uu_Apply(discrep_eval(z, theta_in), u, z), z, u, theta); % B1(theta_zero, z_lofi, u_lofi, theta_zero)
B1 = @(theta_in, z, u, theta) B1_p1(J_uu(discrep_eval(z, theta_in), u, z), z, u, theta); % B1(theta_zero, z_lofi, u_lofi, theta_zero)
B2 = @(theta_in, z, u, theta) discrep_zt(theta_in)' * J_u(u, z);
B = @(theta_in, z, u, theta) B1(theta_in, z, u, theta) + B2(theta_in, z, u, theta);

% Compute action of H
H1 = @(z_in, z, u, theta) opt_prob_interface.Apply_RS_Hessian(z_in, z); % NOTE: this computes J_zz + Sz' * J_uu * Sz
H2_a = @(u_in, z, u, theta) B1_p1(u_in, z, u, theta);
H2_b = @(z_in, z, u, theta) J_uu(discrep_z(theta) * z_in, u, z);
H2 = @(z_in, z, u, theta) H2_a(H2_b(z_in, z, u, theta), z, u, theta);
H3 = @(z_in, z, u, theta) discrep_z(theta)' * J_uu(opt_prob_interface.Apply_Solution_Operator_z_Jacobian(z_in, z), u, z);
H = @(z_in, z, u, theta) H1(z_in, z, u, theta) + H2(z_in, z, u, theta) + H3(z_in, z, u, theta);
H_inv = @(z_in, z, u, theta) Apply_CG_Inverse(H, z_in, z, u, theta);
HinvB = @(theta_in, z, u, theta) H_inv(B(theta_in, z, u, theta), z, u, theta);

% Substitute u = S(z) + delta(z, theta)
B_sub = @(theta_in, z, theta) B(theta_in, z, u_solve(z, theta), theta);
H_sub = @(z_in, z, theta) H(z_in, z, u_solve(z, theta), theta);
Hinv_sub = @(z_in, z, theta) H_inv(z_in, z, u_solve(z, theta), theta);
HinvB_sub = @(theta_in, z, theta) HinvB(theta_in, z, u_solve(z, theta), theta);

% Compute gradient
JJz_tmp = @(u_in, z, u, theta) discrep_z(theta)' * u_in + opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_in, z);
Jz = @(z, u, theta) J_z(u, z) + JJz_tmp(J_u(u, z), z, u, theta);
Jz_sub = @(z, theta) Jz(z, u_solve(z, theta), theta);

% Compute Projected Computations
V = md_hessian_analysis.evecs;
get_z = @(beta) z_lofi + V * beta;
VtB = @(theta_in, beta, u, theta) V' * B(theta_in, get_z(beta), u, theta);
VHVt = @(beta_in, beta, u, theta) V' * H(V * beta_in, get_z(beta), u, theta);
VHVt_inv = @(beta_in, beta, u, theta) Apply_CG_Inverse(VHVt, beta_in, beta, u, theta);
VHVt_inv_VtB = @(theta_in, beta, u, theta) VHVt_inv(VtB(theta_in, beta, u, theta), beta, u, theta);
VHVt_sub = @(beta_in, beta, theta) VHVt(beta_in, beta, u_solve(get_z(beta), theta), theta);
VHVt_inv_sub = @(beta_in, beta, theta) VHVt_inv(beta_in, beta, u_solve(get_z(beta), theta), theta);
VHVt_inv_VtB_sub = @(theta_in, beta, theta) VHVt_inv_VtB(theta_in, beta, u_solve(get_z(beta), theta), theta);

% Need to compare above methods w/ MD_Continuation_Update
num_continuation_steps = 8;
dt = 1 / num_continuation_steps;
ts = dt * (0:num_continuation_steps);
thetas = theta_post * ts;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);

% Diagonostics
disp("Lofi Objective: " + opt_hifi.Jhat(z_lofi));
disp("Hifi Objective: " + opt_hifi.Jhat(z_hifi));
fprintf("\n\n");

% Continuation
cont_method = "fwd-euler-proj";
disp("Continuation (method: " + cont_method + ", num_steps: " + num_continuation_steps + "):");

switch cont_method
    case "fwd-euler"
        zs = {z_lofi};
        for n = 1:num_continuation_steps
            zs{n + 1} = zs{n} - dt * HinvB_sub(theta_post, zs{n}, thetas(:, n));
        end
        z_cont_update = zs{end};
    case "pred-corr"
        zs = {z_lofi};
        for n = 1:num_continuation_steps
            z_pred = zs{n} - dt * HinvB_sub(theta_post, zs{n}, thetas(:, n));
            zs{n + 1} = z_pred - Hinv_sub(Jz_sub(z_pred, thetas(:, n + 1)), z_pred, thetas(:, n + 1));
        end
        z_cont_update = zs{end};
    case "fwd-euler-proj"
        bs = {zeros(size(V, 2), 1)};
        for n = 1:num_continuation_steps
            % bs{n + 1} = bs{n} - dt * VHVt_inv_VtB_sub(theta_post, bs{n}, thetas(:, n));
            u_n = u_solve(get_z(bs{n}), thetas(:, n));
            % disp(norm(u_n))
            t_in = VtB(theta_post, bs{n}, u_n, thetas(:, n)); % difference here...
            theta_in = theta_post;
            z = get_z(bs{n});
            u = u_n;
            theta = thetas(:, n);
            % t_in_1 = B1(theta_post, get_z(bs{n}), u_n, thetas(:, n));
            u_tmp2 = J_uu(discrep_eval(z, theta_in), u, z);
            % disp(norm(u_tmp2))
            % disp(norm(z)-2.9111e3)
            % disp(norm(opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z)))
            % disp(norm(t_in))
            bs_crr = VHVt_inv(t_in, bs{n}, u_n, thetas(:, n));
            % disp(norm(bs_crr))
            bs{n + 1} = bs{n} - dt * bs_crr;
            % disp(norm(bs{n + 1}))
        end
        z_cont_update = get_z(bs{end});
    case "pred-corr-proj"
        bs = {zeros(size(V, 2), 1)};
        for n = 1:num_continuation_steps
            b_pred = bs{n} - dt * VHVt_inv_VtB_sub(theta_post, bs{n}, thetas(:, n));
            bs{n + 1} = b_pred - VHVt_inv_sub(V' * Jz_sub(get_z(b_pred), thetas(:, n + 1)), b_pred, thetas(:, n + 1));
        end
        z_cont_update = get_z(bs{end});
    otherwise
        error("Unknown Continuation Method: " + cont_method);
end

disp("Post-continuation Objective: " + opt_hifi.Jhat(z_cont_update));

[u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_beta();
disp("Post-continuation Objective via Interface: " + opt_hifi.Jhat(z_cont(:, end)));

% md_update = MD_Update(md_post_sampling, md_hessian_analysis);
% z_update_mean = md_update.Posterior_Update_Mean();
% disp("Post-continuation Objective via HDSA: " + opt_hifi.Jhat(z_update_mean));

%% Functions
function out = J_u_fn(obj, u, z)
    [~, out, ~] = obj.J(u, z);
end

function out = J_z_fn(obj, u, z)
    [~, ~, out] = obj.J(u, z);
end

function [z_out] = Apply_CG_Inverse(fn, z_in, z, u, theta)
    z_out = 0 * z_in;
    for k = 1:size(z_in, 2)
        tol = 1.e-7;
        max_iter = length(z) + 10;
        [z_out(:, k), flag, ~, ~, ~] = pcg(@(x)fn(x, z, u, theta), z_in(:, k), tol, max_iter);
        if flag ~= 0
            disp('CG did not converge');
            disp(flag);
        end
    end
end
