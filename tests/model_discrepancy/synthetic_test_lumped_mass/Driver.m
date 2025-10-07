%%
clear;
close all;
addpath(genpath('../../../src'));
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
E_u = u_hyperparam_interface.beta_u * S + M;
M_u_inv = linsolve(M_u,eye(m));
W_u = (1/u_hyperparam_interface.alpha_u) * E_u' * M_u_inv * E_u;
E_u_inv = linsolve(E_u,eye(m));
W_u_inv = u_hyperparam_interface.alpha_u * E_u_inv * M_u * E_u_inv;


%%
diff = [];
u = randn(m,1);
tmp1 = M * u;
tmp2 = u_prior_interface.Apply_M_u(u);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff = [diff ; local_diff];

u = randn(m,1);
tmp1 = W_u_inv * u;
tmp2 = u_prior_interface.Apply_W_u_Inverse(u);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff = [diff ; local_diff];

u = randn(m,1);
scalar = rand;
tmp1 = linsolve(W_u + scalar*M_u,u);
tmp2 = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u,scalar);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff = [diff ; local_diff];

sampling_diff = [];
num_samps = 100000;
R = chol(W_u);
tmp1 = linsolve(R,randn(m,num_samps));                 
test1 = cov(tmp1');
tmp2 = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
test2 = cov(tmp2');
local_diff = norm(test1-test2,'fro')/norm(test1,'fro');
sampling_diff = [sampling_diff ; local_diff];

scalar = rand;
R = chol(W_u + scalar*M_u);
tmp1 = linsolve(R,randn(m,num_samps));  
test1 = cov(tmp1');
tmp2 = u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samps,scalar);
test2 = cov(tmp2');
local_diff = norm(test1-test2,'fro')/norm(test1,'fro');
sampling_diff = [sampling_diff ; local_diff];

if max(diff) > 2.e-3
    disp('Error in model_discrepancy/synthetic_test_lumped_mass');
end

if max(sampling_diff) > 2.e-2
    disp('Error in model_discrepancy/synthetic_test_lumped_mass');
end