clear;
close all;
addpath(genpath('../../../src'));
rng(132253);

n_y = 100;
n_t = 51;
T = 1;
n_z = n_y;
obj = Adv_Diff_Objective(n_y, n_z, T, n_t);
con_hifi = Adv_Diff_Constraint(n_y, n_z, T, n_t);
con = Diff_Constraint(n_y, n_z, T, n_t);
opt = Reduced_Space_Optimization(obj, con);

data_interface = MD_Data_Interface_Transient_Test_Problem();
data_interface.Load_Data();

u_hyperparams = MD_u_Hyperparameters_Transient_Test_Problem(data_interface, n_y);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(u_hyperparams, T, n_t, n_y);
spatial_u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(con_hifi.S, con_hifi.M, u_hyperparams);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(spatial_u_prior_interface, transient_prior_cov);

num_state_solves = 100;
z_hyperparams = MD_z_Hyperparameters_Transient_Test_Problem(data_interface, u_prior_interface, num_state_solves, con, n_y);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(con_hifi.S, con_hifi.M, z_hyperparams);

% transient_prior_cov.Compute_Time_Covariance_GEVP_test();

num_samples = 10;
u_samples = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samples);

t = u_hyperparams.Load_Time_Node_Data();
x = u_hyperparams.Load_Spatial_Node_Data();
b = max(abs(u_samples(:)));

% for k = 1:num_samples
%     tmp = reshape(u_samples(:, k), n_y, n_t);
%     figure;
%     surf(t, x, tmp);
%     xlabel('Time');
%     ylabel('Space');
%     title(['Sample number ', num2str(k)]);
%     view(2);
%     colorbar();
%     caxis([-b, b]);
%     set(gca, 'fontsize', 24);
% end
