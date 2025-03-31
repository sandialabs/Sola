%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
load Assembled_Operators.mat;
rng(2451423);

x = pde_meshing.x;
y = pde_meshing.y;
m = length(x);

con = Diff_Constraint(pde_meshing, diff_coeff);
obj = Diff_Objective(con, reg_coeff);
opt = Reduced_Space_Optimization(obj, con);

data_interface = MD_Data_Interface_Diff();
data_centering = true;

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x, y, data_centering);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(pde_meshing.S, pde_meshing.M, data_interface, u_hyperparam_interface);

num_state_solves = 100;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(num_state_solves, x, y, con);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(pde_meshing.S, pde_meshing.M, data_interface, z_hyperparam_interface, u_prior_interface);

% z_hyperparam_interface.beta_z = (5) * z_hyperparam_interface.beta_z;
% z_prior_interface.Set_beta_z(z_hyperparam_interface.beta_z);
%
% z_hyperparam_interface.discrepancy_percent_z_variation = 100 * z_hyperparam_interface.discrepancy_percent_z_variation;
% z_hyperparam_interface.Determine_alpha_z(z_prior_interface);
% z_prior_interface.Set_alpha_z(z_hyperparam_interface.alpha_z);

md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

md_prior_sampling.Visualization_for_Prior_Discrepancy_at_z_opt(1);

% num_prior_samples = 100;
% num_perts_init = round((18 / pi^2) / z_prior_interface.beta_z);
% [delta_samples_z_opt, delta_samples_z_pert, z_pert] = md_prior_sampling.Prior_Discrepancy_Samples_for_Visualization(num_prior_samples, num_perts_init);
% num_perts = size(z_pert, 2);
%
% z_corr_len = zeros(num_perts, 1);
% delta_pert_mag = zeros(num_prior_samples, num_perts);
% initial_guess = 0;
% for k = 1:num_perts
%     z_corr_len(k) = computeCorrelationLength_2D(x, y, z_pert(:, k), initial_guess);
%     initial_guess = z_corr_len(k);
%
%     for i = 1:num_prior_samples
%         delta_pert_mag(i, k) = max(abs(delta_samples_z_pert{k}(:, i)));
%     end
% end
%
% figure;
% scatter(ones(num_prior_samples, 1) * z_corr_len', delta_pert_mag);
%
% pert_id = 1;
% [~, sample_id] = max(delta_pert_mag(:, pert_id));
%
% name = 'Discrepancy sample at z_{opt}';
% pde_meshing.Plot_Field(delta_samples_z_opt(:, sample_id), name);
%
% name = 'Discrepancy sample at pertubed z 1';
% pde_meshing.Plot_Field(delta_samples_z_pert{pert_id}(:, sample_id), name);
%
% name = 'Perturbed z 1';
% pde_meshing.Plot_Field(z_pert(:, pert_id), name);
%
% pert_id = num_perts;
% [~, sample_id] = max(delta_pert_mag(:, pert_id));
%
% name = 'Discrepancy sample at pertubed z end';
% pde_meshing.Plot_Field(delta_samples_z_pert{pert_id}(:, sample_id), name);
%
% name = 'Perturbed z end';
% pde_meshing.Plot_Field(z_pert(:, pert_id), name);
