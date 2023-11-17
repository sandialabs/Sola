%%
clear;
close all;
clc;
addpath(genpath('../../../src'));
rng(1234423);

m = 200;
diff_coeff = 1;
vel_coeff = 5;
robin_coeff = 2;
reg_coeff = 10;
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(obj, con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_hifi.x;

%%
data_interface = MD_Data_Interface_Adv_Diff();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
alpha_u = 1 / (2^2);
alpha_z = 1 / (600^2);
u_prior_interface = MD_Elliptic_u_Prior_Interface_Adv_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Adv_Diff(alpha_z, opt_lofi);

md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);

%%
continuation_step_sweep = [1, 3, 7, 12, 20];
p = length(continuation_step_sweep);
u = cell(p, 1);
z = cell(p, 1);

error = zeros(p, 1);
z_hifi = load('z_hifi.mat').z_hifi;
for k = 1:p
    num_continuation_steps = continuation_step_sweep(k);
    md_update = MD_Continuation_Update(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, num_continuation_steps);
    alpha_d = 1.e-5;
    num_post_samples = 100;
    md_update.Compute_Posterior_Data(alpha_d, num_post_samples);
    [u{k}, z{k}] = md_update.Posterior_Update_Mean();
    error(k) = sqrt((z{k}(:, end) - z_hifi)' * con_lofi.M * (z{k}(:, end) - z_hifi)) / sqrt(z_hifi' * con_lofi.M * z_hifi);
end

%%
figure;
hold on;
plot(x, data_interface.z_opt, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
lg = cell(p + 2, 1);
lg{1} = 'Low-fidelity solution';
lg{2} = 'High-fidelity solution';
for k = 1:p
    plot(x, z{k}(:, end), '--', 'LineWidth', 3);
    lg{k + 2} = [num2str(continuation_step_sweep(k)), ' step posterior mean update'];
end
legend(lg);
