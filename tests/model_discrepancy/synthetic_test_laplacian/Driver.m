%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

m = 501;

[M, S] = Assemble_Mass_and_Stiffness(m);

data_interface = MD_Data_Interface_synthetic_test_laplacian(m);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_laplacian(m);
u_hyperparam_interface.alpha_u = 0.048969233204560;
u_hyperparam_interface.beta_u = 0.007702351792463;

%%
u_prior_interface = MD_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

M_u = M;
E_u = u_hyperparam_interface.beta_u * S + M;
M_u_inv = linsolve(M_u,eye(m));
W_u = (1/u_hyperparam_interface.alpha_u) * E_u' * M_u_inv * E_u;
E_u_inv = linsolve(E_u,eye(m));
W_u_inv = u_hyperparam_interface.alpha_u * E_u_inv * M_u * E_u_inv;


%%
diff1 = [];
u = randn(m,1);
tmp1 = M * u;
tmp2 = u_prior_interface.Apply_M_u(u);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff1 = [diff1 ; local_diff];

u = randn(m,1);
tmp1 = W_u_inv * u;
tmp2 = u_prior_interface.Apply_W_u_Inverse(u);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff1 = [diff1 ; local_diff];

u = randn(m,1);
scalar = (1.e-2) * rand;
tmp1 = linsolve(W_u + scalar*M_u,u);
tmp2 = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u,scalar);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff1 = [diff1 ; local_diff];

if max(diff1) > 1.e-8
    disp('Error in model_discrepancy/synthetic_test_laplacian');
end

%%
u_prior_interface = MD_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface, true);

%%
diff2 = [];
u = randn(m,1);
tmp1 = M * u;
tmp2 = u_prior_interface.Apply_M_u(u);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff2 = [diff2 ; local_diff];

u = randn(m,1);
tmp1 = W_u_inv * u;
tmp2 = u_prior_interface.Apply_W_u_Inverse(u);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff2 = [diff2 ; local_diff];

u = randn(m,1);
scalar = (1.e-2) * rand;
tmp1 = linsolve(W_u + scalar*M_u,u);
tmp2 = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u,scalar);
local_diff = norm(tmp1 - tmp2)/norm(tmp1);
diff2 = [diff2 ; local_diff];

if max(diff2) > 1.e-4
    disp('Error in model_discrepancy/synthetic_test_laplacian');
end

%%
m = 51;

[M, S] = Assemble_Mass_and_Stiffness(m);

data_interface = MD_Data_Interface_synthetic_test_laplacian(m);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_laplacian(m);
u_hyperparam_interface.alpha_u = 0.048969233204560;
u_hyperparam_interface.beta_u = 0.007702351792463;

%%
u_prior_interface = MD_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

M_u = M;
E_u = u_hyperparam_interface.beta_u * S + M;
M_u_inv = linsolve(M_u,eye(m));
W_u = (1/u_hyperparam_interface.alpha_u) * E_u' * M_u_inv * E_u;
E_u_inv = linsolve(E_u,eye(m));
W_u_inv = u_hyperparam_interface.alpha_u * E_u_inv * M_u * E_u_inv;

%%

sampling_diff1 = [];
num_samps = 10000;
R = chol(W_u);
tmp1 = linsolve(R,randn(m,num_samps));                 
test1 = cov(tmp1');
tmp2 = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
test2 = cov(tmp2');
local_diff = norm(test1-test2,'fro')/norm(test1,'fro');
sampling_diff1 = [sampling_diff1 ; local_diff];

scalar = (1.e-2) * rand;
R = chol(W_u + scalar*M_u);
tmp1 = linsolve(R,randn(m,num_samps));  
test1 = cov(tmp1');
tmp2 = u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samps,scalar);
test2 = cov(tmp2');
local_diff = norm(test1-test2,'fro')/norm(test1,'fro');
sampling_diff1 = [sampling_diff1 ; local_diff];

if max(sampling_diff1) > 5.e-2
    disp('Error in model_discrepancy/synthetic_test_laplacian');
end

%%
u_prior_interface = MD_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface, true);

%%

sampling_diff2 = [];
num_samps = 10000;
R = chol(W_u);
tmp1 = linsolve(R,randn(m,num_samps));                 
test1 = cov(tmp1');
tmp2 = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
test2 = cov(tmp2');
local_diff = norm(test1-test2,'fro')/norm(test1,'fro');
sampling_diff2 = [sampling_diff2 ; local_diff];

scalar = (1.e-2) * rand;
R = chol(W_u + scalar*M_u);
tmp1 = linsolve(R,randn(m,num_samps));  
test1 = cov(tmp1');
tmp2 = u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samps,scalar);
test2 = cov(tmp2');
local_diff = norm(test1-test2,'fro')/norm(test1,'fro');
sampling_diff2 = [sampling_diff2 ; local_diff];

if max(sampling_diff2) > 5.e-2
    disp('Error in model_discrepancy/synthetic_test_laplacian');
end