%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
rng(1234423);

random_numbers = randn(4 * 10^5, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(1234423);

m = 200;
diff_coeff = 1;
vel_coeff = 1 / 2;
robin_coeff = 2;
reg_coeff = 10;
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(obj, con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_hifi.x;

%%
data_interface = MD_Data_Interface_PDE_Test_Problem();
data_interface.Load_Data();

u_lofi = con_lofi.State_Solve(data_interface.Z);

writematrix(data_interface.u_opt, 'u_opt.txt');
writematrix(data_interface.z_opt, 'z_opt.txt');
writematrix(data_interface.D(:, 1) + u_lofi(:, 1), 'U_Hifi_1.txt');
writematrix(u_lofi(:, 1), 'U_Lofi_1.txt');
writematrix(data_interface.Z(:, 1), 'Z_1.txt');
writematrix(data_interface.D(:, 2) + u_lofi(:, 2), 'U_Hifi_2.txt');
writematrix(u_lofi(:, 2), 'U_Lofi_2.txt');
writematrix(data_interface.Z(:, 2), 'Z_2.txt');

alpha_u = 1 / (2^2);
alpha_z = 1 / (3^2);
u_prior_interface = MD_Elliptic_u_Prior_Interface_PDE_Test_Problem(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_PDE_Test_Problem(alpha_z, opt_lofi);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

prior_delta_z_opt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

%%
z = zeros(m, 3);
z(:, 1) = 1 + sin(2 * pi * x);
z(:, 2) = 1 + cos(2 * pi * x);
z(:, 3) = 1 + sin(20 * pi * x);
prior_delta = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);

%%
md_post_samples = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
num_post_samples = 100;
md_post_samples.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(m, 3);
Z_test(:, 1:2) = md_post_samples.post_data.Z;
Z_test(:, 3) = 1.5 * ones(m, 1);
[post_delta_mean, post_delta_samples] = md_post_samples.Posterior_Discrepancy_Samples(Z_test);

%%
opt_prob_interface = MD_Opt_Prob_Interface_Sola(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_samples, md_hessian_analysis);

[post_z_mean, post_z_samples] = md_update.Posterior_Update_Samples();

%%
post_delta_mean = reshape(cell2mat(post_delta_mean), 200, 3);
save('Sola_Output.mat', 'prior_delta', 'prior_delta_z_opt', 'post_delta_mean', 'post_delta_samples', 'post_z_mean', 'post_z_samples');
