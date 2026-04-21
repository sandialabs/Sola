%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
addpath(genpath('../../src'));
rng(132253);

working_path = pwd;
write_path = '~/Documents/dasco/papers/Model_Discrepancy_Hyperparameters/figures/';

n_y = 100;
n_t = 51;
T = 1;
n_z = n_y;
obj = Adv_Diff_Objective(n_y, n_z, T, n_t);
con_hifi = Adv_Diff_Constraint(n_y, n_z, T, n_t);
con = Diff_Constraint(n_y, n_z, T, n_t);
opt = Reduced_Space_Optimization(obj, con);

data_interface = MD_Data_Interface_Transient_Test_Problem();

num_samples = 500;

%%
adapt_time_variance = true;
u_hyperparam_interface = MD_u_Hyperparameter_Interface_Transient_Test_Problem(n_t, n_y, adapt_time_variance);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface, T, n_t, n_y);
spatial_u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(con_hifi.S, con_hifi.M, data_interface, u_hyperparam_interface);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_u_prior_interface, transient_prior_cov);

u_samples_adapt = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samples);

%%
adapt_time_variance = false;
u_hyperparam_interface = MD_u_Hyperparameter_Interface_Transient_Test_Problem(n_t, n_y, adapt_time_variance);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface, T, n_t, n_y);
spatial_u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(con_hifi.S, con_hifi.M, data_interface, u_hyperparam_interface);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_u_prior_interface, transient_prior_cov);

u_samples_no_adapt = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samples);

%%
ts_adapt = zeros(n_t, num_samples);
ts_no_adapt = zeros(n_t, num_samples);
for k = 1:num_samples
    tmp = reshape(u_samples_adapt(:, k), n_y, n_t);
    ts_adapt(:, k) = real(sqrt(diag(tmp' * spatial_u_prior_interface.Apply_M_u(tmp))));

    tmp = reshape(u_samples_no_adapt(:, k), n_y, n_t);
    ts_no_adapt(:, k) = real(sqrt(diag(tmp' * spatial_u_prior_interface.Apply_M_u(tmp))));
end

tmp = reshape(data_interface.D(:, 1), n_y, n_t);
d1_ts = sqrt(diag(tmp' * spatial_u_prior_interface.Apply_M_u(tmp)));

t = u_hyperparam_interface.Load_Time_Node_Data();
b = 1.1 * max([ts_adapt(:); ts_no_adapt(:)]);

figure;
hold on;
plot(t, ts_adapt, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(t, d1_ts, '--', 'LineWidth', 3, 'color', 'red');
xlabel('Time');
ylabel('Discrepancy Norm');
ylim([0, b]);
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'adapted_time_variance', 'epsc');
cd(working_path);

figure;
hold on;
plot(t, ts_no_adapt, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(t, d1_ts, '--', 'LineWidth', 3, 'color', 'red');
xlabel('Time');
ylabel('Discrepancy Norm');
ylim([0, b]);
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'nonadapted_time_variance', 'epsc');
cd(working_path);
