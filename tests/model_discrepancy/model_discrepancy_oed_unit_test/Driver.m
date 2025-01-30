%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

print_output = true; %false;

m = 200;
diff_coeff = 1;
vel_coeff = 1 / 2;
robin_coeff = 2;
reg_coeff = 10;
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
data_interface = MD_Data_Interface_PDE_Test_Problem();
data_interface.Load_Data();

alpha_u = 1 / (2^2);
alpha_z = (199^2) / (600^2);
u_prior_interface = MD_Elliptic_u_Prior_Interface_PDE_Test_Problem(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_PDE_Test_Problem(alpha_z, opt_lofi);

%%
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

oed_interface = MD_OED_Interface_Diff(data_interface, con_lofi);

md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.verbosity = print_output;
md_oed.Offline_Computation();

%%
alpha_d = 1.e-5;
N = 5;
r = length(md_hessian_analysis.evals);
beta = randn((N - 1) * r, 1);
val = md_oed.Evaluate_Posterior_Cov_Trace(beta, alpha_d);

Z = zeros(m, N);
Z(:, 1) = data_interface.z_opt;
for i = 2:N
    Z(:, i) = data_interface.z_opt + md_hessian_analysis.evecs * beta(((i - 2) * r + 1):((i - 1) * r));
end

Ztilde = Z - data_interface.z_opt * ones(1, N);
M_Ztilde = z_prior_interface.Apply_M_z(Ztilde);
G = ones(N, N) + M_Ztilde' * z_prior_interface.Apply_W_z_Inverse(M_Ztilde);
[g, mu] = eig(G, 'vector');

Wu_inv = u_prior_interface.Apply_W_u_Inverse(eye(m));
Wd = u_prior_interface.Apply_M_u(eye(m));
[X, lambda] = eig(Wu_inv * Wd, 'vector');
for j = 1:m
    X(:, j) = X(:, j) / sqrt(X(:, j)' * Wd * X(:, j));
end
lambda = 1 ./ lambda;

y = zeros(m, N);
for i = 1:N
    y(:, i) = Z * g(:, i) - sum(g(:, i)) * data_interface.z_opt;
end

w_ij = cell(N, 1);
tmp1 = opt_prob_interface.Apply_Misfit_Hessian(X, data_interface.u_opt, data_interface.z_opt);
tmp2 = opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp1, data_interface.z_opt);
Ju = opt_prob_interface.Misfit_Gradient(data_interface.u_opt, data_interface.z_opt);
for i = 1:N
    w_ij{i} = sum(g(:, i)) * tmp2;
    Mz_yi = z_prior_interface.Apply_M_z(y(:, i));
    tmp3 = z_prior_interface.Apply_M_z(z_prior_interface.Apply_W_z_Inverse(Mz_yi));
    for j = 1:m
        w_ij{i}(:, j) =  w_ij{i}(:, j) + (Ju' * X(:, j)) * tmp3;
    end
end

f = 0;
for i = 1:N
    for j = 1:m
        tmp = md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(w_ij{i}(:, j));
        f = f + (1 / (mu(i) + alpha_d * lambda(j))) * (1 / lambda(j)) * (tmp' * z_prior_interface.Apply_W_z(tmp));
    end
end
f = f / (data_interface.z_opt' * z_prior_interface.Apply_W_z(data_interface.z_opt));

error = abs(f - val) / abs(f);
if print_output
    disp(['error = ', num2str(error)]);
end

if error > 1.e-14
    disp('Error in model discrepancy OED unit test');
end

%%
beta = randn((N - 1) * r, 1);
[g, mu, Mg, g_jac, mu_jac, Mg_jac] = md_oed.G_eigs(beta);

i = 3;
dbeta = randn((N - 1) * r, 1);
h = 10.^(-2:-1:-6)';
p = length(h);
g_fd_error = zeros(p, 1);
mu_fd_error = zeros(p, 1);
Mg_fd_error = zeros(p, 1);
for s = 1:p
    [g_s, mu_s, Mg_s] = md_oed.G_eigs(beta + h(s) * dbeta);
    sign_normalization = sign(g(1,i))*sign(g_s(1,i));
    g_s = sign_normalization * g_s;
    Mg_s = sign_normalization * Mg_s;

    fd_approx = (g_s(:, i) - g(:, i)) / h(s);
    g_fd_error(s) = norm(fd_approx - g_jac{i} * dbeta);

    fd_approx = (mu_s(i) - mu(i)) / h(s);
    mu_fd_error(s) = norm(fd_approx - mu_jac{i}' * dbeta);

    fd_approx = (Mg_s(:, i) - Mg(:, i)) / h(s);
    Mg_fd_error(s) = norm(fd_approx - Mg_jac{i} * dbeta);
end

if print_output
    disp('G eigenvector finite difference error:');
    for s = 1:p
        disp(['Stepsize = ', num2str(h(s)), ' and error = ', num2str(g_fd_error(s))]);
    end

    disp('G eigenvalue finite difference error:');
    for s = 1:p
        disp(['Stepsize = ', num2str(h(s)), ' and error = ', num2str(mu_fd_error(s))]);
    end

    disp('M*gi finite difference error:');
    for s = 1:p
        disp(['Stepsize = ', num2str(h(s)), ' and error = ', num2str(Mg_fd_error(s))]);
    end
end

%%
beta = randn((N - 1) * r, 1);
alpha_d = rand;
[val, grad] = md_oed.Evaluate_Posterior_Cov_Trace(beta, alpha_d);

beta_d = randn((N - 1) * r, 1);
h = 10.^(-2:-1:-6);
p = length(h);
post_cov_fd_error = zeros(p, 1);
for s = 1:p
    beta_s = beta + h(s) * beta_d;
    val_s = md_oed.Evaluate_Posterior_Cov_Trace(beta_s, alpha_d);
    fd_approx = (val_s - val) / h(s);
    post_cov_fd_error(s) = abs(fd_approx - grad' * beta_d);
end

if print_output
    disp('Posterior covariance finite difference error:');
    for s = 1:p
        disp(['Stepsize = ', num2str(h(s)), ' and error = ', num2str(post_cov_fd_error(s))]);
    end
end

%%
beta = randn((N - 1) * r, 1);
alpha_d = rand;
reg_coeff = rand;
[val, grad] = md_oed.Evaluate_OED_Objective(beta, alpha_d, reg_coeff);

beta_d = randn((N - 1) * r, 1);
h = 10.^(-2:-1:-6);
p = length(h);
obj_fd_error = zeros(p, 1);
for s = 1:p
    beta_s = beta + h(s) * beta_d;
    val_s = md_oed.Evaluate_OED_Objective(beta_s, alpha_d, reg_coeff);
    fd_approx = (val_s - val) / h(s);
    obj_fd_error(s) = abs(fd_approx - grad' * beta_d);
end

if print_output
    disp('OED objective finite difference error:');
    for s = 1:p
        disp(['Stepsize = ', num2str(h(s)), ' and error = ', num2str(obj_fd_error(s))]);
    end
end

%%
beta_0 = randn((N - 1) * r, 1);
alpha_d = rand;
reg_coeff = (1.e-6) * rand;
[beta, Z] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);

%%
% save('Reference_Solution.mat','g_fd_error','mu_fd_error','Mg_fd_error','post_cov_fd_error','obj_fd_error','beta','Z')
ref = load('Reference_Solution.mat');

error = norm(ref.Mg_fd_error - Mg_fd_error);
error = max(error, norm(ref.mu_fd_error - mu_fd_error));
error = max(error, norm(ref.g_fd_error - g_fd_error));
error = max(error, norm(ref.post_cov_fd_error - post_cov_fd_error));
error = max(error, norm(ref.obj_fd_error - obj_fd_error));
error = max(error, norm(ref.beta - beta));
error = max(error, norm(ref.Z - Z));

if error > 1.e-14
    disp('Error in model discrepancy OED unit test');
end
