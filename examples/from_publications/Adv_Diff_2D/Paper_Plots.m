%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

addpath(genpath('../../src'));
load('HDSA_Results.mat');

working_path = pwd;
write_path = '~/Desktop/output_figures/adv_diff/';

name = 'Low-fidelity state';
u_lofi = adv_diff.State_Solve(opt.con.Map_z_vec_to_mesh(z_lofi));
adv_diff.pde_meshing.Plot_Field(u_lofi, name);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$\tilde{S}(\tilde{z})$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'lofi_state', 'epsc');
cd(working_path);

name = 'Low-fidelity control';
adv_diff.pde_meshing.Plot_Field(opt.con.Map_z_vec_to_mesh(z_lofi), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$f(\tilde{z})$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'lofi_source', 'epsc');
cd(working_path);

name = 'High-fidelity state';
u_hifi = nonlinear_adv_diff.State_Solve(opt.con.Map_z_vec_to_mesh(z_lofi));
adv_diff.pde_meshing.Plot_Field(u_hifi, name);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$S(\tilde{z})$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'hifi_state', 'epsc');
cd(working_path);

name = 'Updated control mean';
adv_diff.pde_meshing.Plot_Field(opt.con.Map_z_vec_to_mesh(z_update_mean_1), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
caxis([zmin, zmax]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$f(\overline{z})$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'mean_source_update_1', 'epsc');
cd(working_path);

name = 'Updated control standard deviation';
adv_diff.pde_meshing.Plot_Field(std(opt.con.Map_z_vec_to_mesh(z_update_samples_1), [], 2), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$f(z)$ Posterior Pointwise Standard Deviation', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'std_source_update_1', 'epsc');
cd(working_path);

name = 'Updated control mean';
adv_diff.pde_meshing.Plot_Field(opt.con.Map_z_vec_to_mesh(z_update_mean_2), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
caxis([zmin, zmax]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$f(\overline{z})$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'mean_source_update_2', 'epsc');
cd(working_path);

name = 'Updated control standard deviation';
adv_diff.pde_meshing.Plot_Field(std(opt.con.Map_z_vec_to_mesh(z_update_samples_2), [], 2), name);
xlim(opt.obj.control_xlim);
ylim(opt.obj.control_ylim);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
title('$f(z)$ Posterior Pointwise Standard Deviation', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'std_source_update_2', 'epsc');
cd(working_path);

edges = linspace(0.5, 4, 10)' * (10^-3);
figure;
hold on;
histogram(val_update_samples_1, edges);
plot([val_update_1, val_update_1], [0, 35], '--', 'LineWidth', 3, 'Color', 'red');
plot([val_lofi, val_lofi], [0, 35], 'LineWidth', 3, 'Color', 'black');
xlim([0.5, 4] * (10^-3));
ylim([0, 40]);
xlabel('High-fidelity objective function value');
% yticks([])
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_obj_fun_vals_rank_1', 'epsc');
cd(working_path);

figure;
hold on;
histogram(val_update_samples_2, edges);
plot([val_update_2, val_update_2], [0, 35], '--', 'LineWidth', 3, 'Color', 'red');
plot([val_lofi, val_lofi], [0, 35], 'LineWidth', 3, 'Color', 'black');
xlim([0.5, 4] * (10^-3));
ylim([0, 40]);
xlabel('High-fidelity objective function value');
% yticks([])
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_obj_fun_vals_rank_2', 'epsc');
cd(working_path);
