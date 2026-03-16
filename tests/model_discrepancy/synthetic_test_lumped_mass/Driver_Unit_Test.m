%%
clear;
close all;
rng(121234);

m = 101;
x = linspace(0, 1, m)';

[M, S] = Assemble_Mass_and_Stiffness(m);

data_interface = MD_Data_Interface_synthetic_test_lumped_mass(m);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_lumped_mass(m);
u_hyperparam_interface.alpha_u = 0.048969233204560;
u_hyperparam_interface.beta_u = 0.007702351792463;
u_prior_interface = MD_Lumped_Mass_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_synthetic_test_lumped_mass(m);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

M_u = M;
M_lumped = diag(M * ones(m, 1));
M_lumped_inv = linsolve(M_lumped, eye(m));
E_u = u_hyperparam_interface.beta_u * S + M;
M_u_inv = linsolve(M, eye(m));
W_u = (1 / u_hyperparam_interface.alpha_u) * E_u' * M_lumped_inv * E_u;
E_u_inv = linsolve(E_u, eye(m));
W_u_inv = u_hyperparam_interface.alpha_u * E_u_inv * M_lumped * E_u_inv;

%%
diff = [];
lumped_mass_error = [];

u = linsolve(E_u, randn(m, 1));
tmp1 = M * u;
tmp2 = u_prior_interface.Apply_M_u(u);
local_diff = norm(tmp1 - tmp2) / norm(tmp1);
diff = [diff; local_diff];

u = linsolve(E_u, randn(m, 1));
tmp1 = W_u_inv * u;
tmp2 = u_prior_interface.Apply_W_u_Inverse(u);
local_diff = norm(tmp1 - tmp2) / norm(tmp1);
diff = [diff; local_diff];

tmp3 = u_hyperparam_interface.alpha_u * E_u_inv * M_u * E_u_inv * u;
lumped_mass_error_1 = norm(tmp3 - tmp2) / norm(tmp3);
lumped_mass_error = [lumped_mass_error; lumped_mass_error_1];

u = linsolve(E_u, randn(m, 1));
scalar = 5e4;
tmp1 = linsolve(W_u + scalar * M_u, u);
tmp2 = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u, scalar);
local_diff = norm(tmp1 - tmp2) / norm(tmp1);
diff = [diff; local_diff];

tmp3 = linsolve((1 / u_hyperparam_interface.alpha_u) * E_u * M_u_inv * E_u + scalar * M_u, u);
lumped_mass_error_2 = norm(tmp3 - tmp2) / norm(tmp3);
lumped_mass_error = [lumped_mass_error; lumped_mass_error_2];

sampling_diff = [];
num_samps = 10000;
R = chol(W_u);
tmp1 = linsolve(R, randn(m, num_samps));
test1 = cov(tmp1');
tmp2 = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
test2 = cov(tmp2');
local_diff = norm(test1 - test2, 'fro') / norm(test1, 'fro');
sampling_diff = [sampling_diff; local_diff];

R = chol((1 / u_hyperparam_interface.alpha_u) * E_u' * M_u_inv * E_u);
tmp3 = linsolve(R, randn(m, num_samps));
test3 = cov(tmp3');
lumped_mass_error_3 = norm(test3 - test2, 'fro') / norm(test3, 'fro');
lumped_mass_error = [lumped_mass_error; lumped_mass_error_3];

R = chol(W_u + scalar * M_u);
tmp1 = linsolve(R, randn(m, num_samps));
test1 = cov(tmp1');
tmp2 = u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samps, scalar);
test2 = cov(tmp2');
local_diff = norm(test1 - test2, 'fro') / norm(test1, 'fro');
sampling_diff = [sampling_diff; local_diff];

R = chol((1 / u_hyperparam_interface.alpha_u) * E_u' * M_u_inv * E_u + scalar * M_u);
tmp3 = linsolve(R, randn(m, num_samps));
test3 = cov(tmp3');
lumped_mass_error_4 = norm(test3 - test2, 'fro') / norm(test3, 'fro');
lumped_mass_error = [lumped_mass_error; lumped_mass_error_4];

disp('diff:');
disp(diff');
disp('sampling_diff');
disp(sampling_diff');
disp('lumped_mass_error');
disp(lumped_mass_error');
