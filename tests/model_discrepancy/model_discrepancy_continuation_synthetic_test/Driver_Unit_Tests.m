%%
clear;
close all;
clc;
addpath(genpath('../../../src'));
rng(121234);

print_output = true;
max_error = 0;

con = Synthetic_Test_Constraint();
con_hifi = Synthetic_Test_Hifi_Constraint();
obj = Synthetic_Test_Objective(con);
opt = Reduced_Space_Optimization(obj, con);

data_interface = MD_Data_Interface_continuation_synthetic_test();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);
alpha_u = 1.e1;
alpha_z = 1.e-2;
u_prior_interface = MD_Elliptic_u_Prior_Interface_continuation_synthetic_test(alpha_u, con.m);
z_prior_interface = MD_Elliptic_z_Prior_Interface_continuation_synthetic_test(alpha_z, con.n);

md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);

num_continuation_steps = 5;
md_continuation_update = MD_Continuation_Update(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, num_continuation_steps);

alpha_d = 1.e-4;
num_post_samples = 100;
md_continuation_update.Compute_Posterior_Data(alpha_d, num_post_samples);
md_update.Compute_Posterior_Data(alpha_d, num_post_samples);

%%
delta_true = @(z) -(1 / 2) * linsolve(con.A, con.B) * z;
z_test = randn(con.n, 5);
error = norm(delta_true(z_test) - (con_hifi.State_Solve(z_test) - con.State_Solve(z_test)));
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
I_true = zeros(con.m, 1);
K_true = delta_true(eye(9)) * inv(z_prior_interface.M_z);
tmp = K_true';
theta_true = [I_true; tmp(:)];
z_test = randn(con.n, 1);
error = norm(delta_true(z_test) - [eye(con.m), kron(eye(con.m), z_test' * z_prior_interface.M_z)] * theta_true);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
m = con.m;
n = con.n;
p = m * (n + 1);
N = md_continuation_update.post_data.N;
theta_est = zeros(p, 1);
for ell = 1:N
    coeff = md_continuation_update.post_data.a_ell(ell);
    u = md_continuation_update.post_data.u_ell(:, ell);
    z_tmp = md_continuation_update.post_data.Z(:, ell) - data_interface.z_opt;
    z = linsolve(z_prior_interface.M_z, z_prior_interface.Apply_W_z_Inverse(z_tmp));
    tmp = [coeff * u; kron(u, z)];
    theta_est = theta_est + tmp;

    for i = 1:N
        coeff = md_continuation_update.si(i);
        u = md_continuation_update.post_data.u_i_ell{i}(:, ell);
        z_tmp = md_continuation_update.W_z_inv_yi(:, i);
        z = linsolve(z_prior_interface.M_z, z_tmp);
        tmp = [coeff * u; kron(u, z)];
        theta_est = theta_est - md_continuation_update.post_data.b_i_ell(i, ell) * tmp;
    end
end
theta_est = (1 / alpha_d) * theta_est;

I_est = theta_est(1:m);
K_est = reshape(theta_est((m + 1):end), n, m)';
delta_est = @(z) I_est + K_est * z_prior_interface.M_z * z;

z_test = randn(n, 1);
tmp = md_update.Posterior_Discrepancy_Samples(z_test);
error = norm(tmp{1} - delta_est(z_test));
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
t_n = rand;
delta = md_continuation_update.Discrepancy_Evaluation(z_test, t_n);
error = norm(t_n * tmp{1} - delta);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
t_n = rand;
delta_z_test = md_continuation_update.Apply_Discrepancy_z_Jacobian(z_test, t_n);
error = norm(delta_z_test - t_n * K_est * z_prior_interface.M_z * z_test);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
t_n = rand;
u_test = randn(m, 1);
delta_u_test = md_continuation_update.Apply_Discrepancy_z_Jacobian_transpose(u_test, t_n);
error = norm(delta_u_test - t_n * z_prior_interface.M_z * K_est' * u_test);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
z_test = randn(n, 1);
delta_theta_jac = [eye(m), kron(eye(m), z_test' * z_prior_interface.M_z)];
u_out = md_continuation_update.Apply_Discrepancy_theta_Jacobian(z_test);
error = norm(u_out - delta_theta_jac * theta_est);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
u_test = randn(m, 1);
z_out = md_continuation_update.Apply_Discrepancy_z_theta_Hessian(u_test);

z_out_test = zeros(n, 1);
nom = [u_test', kron(u_test', z_test' * z_prior_interface.M_z)] * theta_est;
for k = 1:n
    e = zeros(n, 1);
    e(k) = 1;
    z_out_test(k) = [u_test', kron(u_test', (z_test + e)' * z_prior_interface.M_z)] * theta_est - nom;
end
error = norm(z_out - z_out_test);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
z_n = randn(n, 1);
t_n = rand;
u_n = con.State_Solve(z_n);
F = linsolve(con.A, con.B) + t_n * K_est * z_prior_interface.M_z;
d = obj.uT - t_n * I_est;
Jhat = @(z) (1 / 2) * (F * z - d)' * obj.M * (F * z - d) + (1 / 2) * z' * obj.R * z;
Jhat_grad = @(z)F' * obj.M * (F * z - d) + obj.R * z;
Jhat_hess = F' * obj.M * F + obj.R;
z_out = md_continuation_update.Apply_Parameterized_RS_Hessian(z_test, u_n, z_n, t_n);
error = norm(z_out - Jhat_hess * z_test);
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
z_out = md_continuation_update.Apply_Parameterized_RS_Hessian_Inverse(z_test, u_n, z_n, t_n);
error = norm(z_out - linsolve(Jhat_hess, z_test));
if print_output
    disp(['error = ', num2str(error)]);
end
max_error = max(max_error, error);

%%
Btheta = md_continuation_update.Apply_B(u_n, z_n, t_n);
h = 10.^(-2:-1:-6);
p = length(h);
fd_error = zeros(p, 1);
for k = 1:p
    F_pert = linsolve(con.A, con.B) + (t_n + h(k)) * K_est * z_prior_interface.M_z;
    d_pert = obj.uT - (t_n + h(k)) * I_est;
    Jhat_grad_pert = @(z)F_pert' * obj.M * (F_pert * z - d_pert) + obj.R * z;
    Btheta_fd = (Jhat_grad_pert(z_n) - Jhat_grad(z_n)) / h(k);
    fd_error(k) = norm(Btheta_fd - Btheta);
end

disp('Finite difference test for B');
for k = 1:p
    disp(['Finite difference error = ', num2str(fd_error(k)), ' for step size ', num2str(h(k))]);
end
