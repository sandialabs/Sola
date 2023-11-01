%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
rng(2451423);

x = adv_diff.pde_meshing.x;
y = adv_diff.pde_meshing.y;
M = adv_diff.pde_meshing.M;
m = length(x);
obj = Adv_Diff_Objective(adv_diff, reg_coeff);
con = Adv_Diff_Constraint(adv_diff);
opt = Reduced_Space_Optimization(obj, con);

alpha_u = 2^2;
alpha_z = 1.e-8;
md_interface = Adv_Diff_HDSA(opt, alpha_u, alpha_z);

%%
num_prior_samples = 10;
md_prior_sampling = HDSA_MD_Prior_Sampling(md_interface);

delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
for k = 1:10
    name = ['Prior discrepancy sample ', num2str(k)];
    adv_diff.pde_meshing.Plot_Field(delta_samples(:, k), name);
end

%%
md_update = HDSA_MD_Update(md_interface);

alpha_d = 1.e-2;
num_post_samples = 100;
md_update.Compute_Posterior_Data(alpha_d, num_post_samples);
n = length(z_lofi);
Z_test = zeros(n, 2);
Z_test(:, 1) = z_lofi;
Z_test(:, 2) = mean(z_lofi) * ones(n, 1);
[delta_mean, delta_samples] = md_update.Posterior_Discrepancy_Samples(Z_test);

%%
name = 'D_1';
adv_diff.pde_meshing.Plot_Field(D(:, 1), name);

name = 'D_1 discrepancy mean';
adv_diff.pde_meshing.Plot_Field(delta_mean{1}, name);

for k = 1:5
    name = 'D_1 discrepancy sample';
    adv_diff.pde_meshing.Plot_Field(delta_samples{1}(:, k), name);
end

diff = D(:, 1) - delta_mean{1};
normalize = sqrt(D(:, 1)' * md_interface.Apply_M_u(D(:, 1)));
mean_diff = sqrt(diff' * md_interface.Apply_M_u(diff)) / normalize;
sample_diff = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    diff = delta_mean{1} - delta_samples{1}(:, k);
    sample_diff(k) = sqrt(diff' * md_interface.Apply_M_u(diff)) / normalize;
end
figure;
hold on;
histogram(sample_diff);
title(['Mean discrepancy error = ', num2str(mean_diff)]);

%%
for k = 1:5
    name = 'Discrepancy at z sample';
    adv_diff.pde_meshing.Plot_Field(delta_samples{2}(:, k), name);
end

normalize = sqrt(delta_mean{2}' * md_interface.Apply_M_u(delta_mean{2}));
sample_diff = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    diff = delta_mean{2} - delta_samples{2}(:, k);
    sample_diff(k) = sqrt(diff' * md_interface.Apply_M_u(diff)) / normalize;
end
figure;
hold on;
histogram(sample_diff);

%%
num_evals = 1;
oversampling = 10;
md_update.Compute_Hessian_GEVP(num_evals, oversampling);

%%
[z_update_mean_1, z_update_samples_1] = md_update.Posterior_Update_Samples();
I = find(x > opt.obj.control_xlim(1));
I = intersect(I, find(x < opt.obj.control_xlim(2)));
I = intersect(I, find(y > opt.obj.control_ylim(1)));
I = intersect(I, find(y < opt.obj.control_ylim(2)));
tmp = opt.obj.Map_z_vec_to_mesh(z_lofi);
zmin = min(tmp(I));
zmax = max(tmp(I));
tmp = opt.obj.Map_z_vec_to_mesh(z_update_mean_1);
zmin = min([zmin; tmp(I)]);
zmax = max([zmax; tmp(I)]);

name = 'Low-fidelity control';
adv_diff.pde_meshing.Plot_Field(opt.obj.Map_z_vec_to_mesh(z_lofi), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
caxis([zmin, zmax]);

name = 'Updated control mean';
adv_diff.pde_meshing.Plot_Field(opt.obj.Map_z_vec_to_mesh(z_update_mean_1), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
caxis([zmin, zmax]);

name = 'Updated control standard deviation';
adv_diff.pde_meshing.Plot_Field(std(opt.obj.Map_z_vec_to_mesh(z_update_samples_1), [], 2), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);

%%
u_lofi = nonlinear_adv_diff.State_Solve(opt.obj.Map_z_vec_to_mesh(z_lofi));
u_update_mean_1 = nonlinear_adv_diff.State_Solve(opt.obj.Map_z_vec_to_mesh(z_update_mean_1));

val_lofi = opt.obj.J(u_lofi, z_lofi);
val_update_1 = opt.obj.J(u_update_mean_1, z_update_mean_1);

u_update_samples_1 = zeros(m, num_post_samples);
val_update_samples_1 = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    u_update_samples_1(:, k) = nonlinear_adv_diff.State_Solve(opt.obj.Map_z_vec_to_mesh(z_update_samples_1(:, k)));
    val_update_samples_1(k) = opt.obj.J(u_update_samples_1(:, k), z_update_samples_1(:, k));
end

disp(['Objective at low-fidelity solution = ', num2str(val_lofi)]);
disp(['Objective at mean update solution = ', num2str(val_update_1)]);

%%
num_evals = 2;
oversampling = 10;
md_update.Compute_Hessian_GEVP(num_evals, oversampling);

%%
[z_update_mean_2, z_update_samples_2] = md_update.Posterior_Update_Samples();
I = find(x > opt.obj.control_xlim(1));
I = intersect(I, find(x < opt.obj.control_xlim(2)));
I = intersect(I, find(y > opt.obj.control_ylim(1)));
I = intersect(I, find(y < opt.obj.control_ylim(2)));
tmp = opt.obj.Map_z_vec_to_mesh(z_lofi);
zmin = min(tmp(I));
zmax = max(tmp(I));
tmp = opt.obj.Map_z_vec_to_mesh(z_update_mean_2);
zmin = min([zmin; tmp(I)]);
zmax = max([zmax; tmp(I)]);

name = 'Low-fidelity control';
adv_diff.pde_meshing.Plot_Field(opt.obj.Map_z_vec_to_mesh(z_lofi), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
caxis([zmin, zmax]);

name = 'Updated control mean';
adv_diff.pde_meshing.Plot_Field(opt.obj.Map_z_vec_to_mesh(z_update_mean_2), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
caxis([zmin, zmax]);

name = 'Updated control standard deviation';
adv_diff.pde_meshing.Plot_Field(std(opt.obj.Map_z_vec_to_mesh(z_update_samples_2), [], 2), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);

%%
u_lofi = nonlinear_adv_diff.State_Solve(opt.obj.Map_z_vec_to_mesh(z_lofi));
u_update_mean_2 = nonlinear_adv_diff.State_Solve(opt.obj.Map_z_vec_to_mesh(z_update_mean_2));

val_lofi = opt.obj.J(u_lofi, z_lofi);
val_update_2 = opt.obj.J(u_update_mean_2, z_update_mean_2);

u_update_samples_2 = zeros(m, num_post_samples);
val_update_samples_2 = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    u_update_samples_2(:, k) = nonlinear_adv_diff.State_Solve(opt.obj.Map_z_vec_to_mesh(z_update_samples_2(:, k)));
    val_update_samples_2(k) = opt.obj.J(u_update_samples_2(:, k), z_update_samples_2(:, k));
end

disp(['Objective at low-fidelity solution = ', num2str(val_lofi)]);
disp(['Objective at mean update solution = ', num2str(val_update_2)]);

%%
save('HDSA_Results.mat');
