%%
clear;
close all;
clc;
addpath(genpath('../../../src'));
rng(121234);

con = Synthetic_Test_Constraint();
con_hifi = Synthetic_Test_Hifi_Constraint();
obj = Synthetic_Test_Objective(con);
opt = Reduced_Space_Optimization(obj, con);

alpha_u = 1.e1;
alpha_z = 1.e-2;
md_interface = Synthetic_Test_MD_Interface_Elliptic_Prior(opt, alpha_u, alpha_z, con, con_hifi, opt);
md_update = HDSA_MD_Update(md_interface);

md_continuation_interface = HDSA_Sabl_MD_Continuation_Interface(md_interface, con);

num_continuation_steps = 100;
md_continuation_update = HDSA_MD_Continuation_Update(md_continuation_interface, num_continuation_steps);

alpha_d = 1.e-4;
num_post_samples = 100;
md_continuation_update.Compute_Posterior_Data(alpha_d, num_post_samples);
md_update.Compute_Posterior_Data(alpha_d, num_post_samples);

%%
[u_update, z_update] = md_continuation_update.Posterior_Update_Mean();

%%
m = con.m;
n = con.n;
N = md_update.post_data.N;
p = m * (n + 1);
theta_est = zeros(p, 1);
for ell = 1:N
    coeff = md_continuation_update.post_data.a_ell(ell);
    u = md_continuation_update.post_data.u_ell(:, ell);
    z_tmp = md_continuation_update.post_data.Z(:, ell) - md_interface.Load_Optimal_z();
    z = linsolve(md_interface.M_z, md_interface.Apply_W_z_Inverse(z_tmp));
    tmp = [coeff * u; kron(u, z)];
    theta_est = theta_est + tmp;

    for i = 1:N
        coeff = md_continuation_update.si(i);
        u = md_continuation_update.post_data.u_i_ell{i}(:, ell);
        z_tmp = md_continuation_update.W_z_inv_yi(:, i);
        z = linsolve(md_interface.M_z, z_tmp);
        tmp = [coeff * u; kron(u, z)];
        theta_est = theta_est - md_continuation_update.post_data.b_i_ell(i, ell) * tmp;
    end
end
theta_est = (1 / alpha_d) * theta_est;

I_est = theta_est(1:m);
K_est = reshape(theta_est((m + 1):end), n, m)';
F = linsolve(con.A, con.B) + K_est * md_interface.M_z;
d = obj.uT - I_est;
Jhat = @(z) (1 / 2) * (F * z - d)' * obj.M * (F * z - d) + (1 / 2) * z' * obj.R * z;
Jhat_grad = @(z)F' * obj.M * (F * z - d) + obj.R * z;
Jhat_hess = F' * obj.M * F + obj.R;

g = Jhat_grad(zeros(n, 1));
z_pert_star = linsolve(Jhat_hess, -g);

%%
z = load('Opt_Data.mat').z;
z_hifi = load('Opt_Data.mat').z_hifi;

figure;
hold on;
plot(z, 'LineWidth', 3);
plot(z_pert_star, 'LineWidth', 3);
plot(z_hifi, 'LineWidth', 3);
plot(z_update(:, end), '--', 'LineWidth', 3);
legend({'Low-fidelity solution', 'Estimated high-fidelity solution', 'High-fidelity solution', 'Mean update'});
